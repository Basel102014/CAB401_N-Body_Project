// nbody_parallel.c  — high-performance cond-var + SoA + cache blocking
#define _POSIX_C_SOURCE 200809L
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <sched.h>

#include "nbody_parallel.h"  // brings Body, Snapshots, G, EPS2

/* ------------------------ Timing ------------------------ */
static inline double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* --------------------- Tuning knobs --------------------- */
#ifndef NBLOCK
#define NBLOCK 128              // try 64, 128, 256 for your CPU
#endif

/* ---------------- Optional core pinning ----------------- */
static inline void maybe_pin_to_core(int want, int core_id) {
#ifdef __linux__
    if (!want) return;
    if (core_id < 0) core_id = 0;
    cpu_set_t set; CPU_ZERO(&set);
    CPU_SET((unsigned)core_id, &set);
    (void)pthread_setaffinity_np(pthread_self(), sizeof(set), &set);
#else
    (void)want; (void)core_id;
#endif
}

/* ------------------- Read-only snapshot ----------------- */
typedef struct {
    const double *x;
    const double *y;
    const double *z;
    const double *m;
} ROSnap;

/* --------------------- Worker context ------------------- */
typedef struct {
    Body   *restrict b;     // shared bodies (write-only: fx/fy/fz in our i-range)
    size_t n;
    size_t start, end;      // i-range [start, end)
    size_t steps;           // total steps (for loop bound)
    double dt;

    const double *restrict rx;  // per-step SoA snapshot (read-only)
    const double *restrict ry;
    const double *restrict rz;
    const double *restrict rm;

    // sync
    struct Sync* sync;
    size_t tid, T;

#ifdef __linux__
    int want_affinity;
    int core_id;
#endif
} TD;

/* ------------------- Step synchronisation ----------------
   Single broadcast from main per step; workers track an epoch.
-----------------------------------------------------------*/
typedef struct Sync {
    pthread_mutex_t m;
    pthread_cond_t  cv_start;   // broadcast to start a new step
    pthread_cond_t  cv_done;    // signal when remaining reaches 0
    unsigned        epoch;      // increments each step
    size_t          remaining;  // workers yet to finish this step
    int             stop;       // shutdown flag
} Sync;

/* ----------------------- Worker ------------------------- */
static void *worker(void *arg) {
    TD *td = (TD*)arg;
    Body *restrict b = td->b;
    const size_t n   = td->n;
    const size_t is  = td->start;
    const size_t ie  = td->end;

#ifdef __linux__
    maybe_pin_to_core(1, (int)td->tid);
#endif

    unsigned seen_epoch = 0;

    for (;;) {
        // Wait for next epoch
        pthread_mutex_lock(&td->sync->m);
        while (td->sync->epoch == seen_epoch && !td->sync->stop)
            pthread_cond_wait(&td->sync->cv_start, &td->sync->m);
        if (td->sync->stop) {
            pthread_mutex_unlock(&td->sync->m);
            break;
        }
        seen_epoch = td->sync->epoch;

        // Snapshot pointers for this epoch (published by main)
        const double *restrict X = td->rx;
        const double *restrict Y = td->ry;
        const double *restrict Z = td->rz;
        const double *restrict M = td->rm;
        pthread_mutex_unlock(&td->sync->m);

        // ---- Force phase (cache-blocked over j) ----
        for (size_t i = is; i < ie; ++i) {
            const double xi = X[i], yi = Y[i], zi = Z[i];
            const double mi = M[i];

            double fx = 0.0, fy = 0.0, fz = 0.0;

            for (size_t jb = 0; jb < n; jb += (size_t)NBLOCK) {
                const size_t jend = (jb + (size_t)NBLOCK < n) ? (jb + (size_t)NBLOCK) : n;

                for (size_t j = jb; j < jend; ++j) {
                    if (j == i) continue;

                    const double dx = X[j] - xi;
                    const double dy = Y[j] - yi;
                    const double dz = Z[j] - zi;

                    const double r2   = dx*dx + dy*dy + dz*dz + EPS2;
                    const double invr = 1.0 / sqrt(r2);
                    const double invr3 = invr * invr * invr;

                    const double sF = G * mi * M[j] * invr3;
                    fx += sF * dx;
                    fy += sF * dy;
                    fz += sF * dz;
                }
            }

            // Write results directly to bodies (no extra copy pass)
            b[i].fx = fx;
            b[i].fy = fy;
            b[i].fz = fz;
        }

        // Done with this step
        pthread_mutex_lock(&td->sync->m);
        if (--td->sync->remaining == 0)
            pthread_cond_signal(&td->sync->cv_done);
        pthread_mutex_unlock(&td->sync->m);
    }

    return NULL;
}

/* ----------------- Public Entry Point (API) -------------- */
int run_nbody_parallel(Body *bodies,
                       size_t n,
                       size_t steps,
                       double dt,
                       size_t num_threads,
                       double *avg_force_time,
                       double *avg_update_time,
                       size_t stride,
                       Snapshots *snaps,
                       double *total_time_out)
{
    if (!bodies || n == 0 || steps == 0 || num_threads == 0) return -1;

    // Allocate SoA snapshot once (aligned for cache/TLB friendliness)
    double *rx = NULL, *ry = NULL, *rz = NULL, *rm = NULL;
    if (posix_memalign((void**)&rx, 64, n*sizeof(double)) ||
        posix_memalign((void**)&ry, 64, n*sizeof(double)) ||
        posix_memalign((void**)&rz, 64, n*sizeof(double)) ||
        posix_memalign((void**)&rm, 64, n*sizeof(double))) {
        perror("posix_memalign");
        free(rx); free(ry); free(rz); free(rm);
        return -1;
    }

    // Seed snapshot for step 0
    for (size_t i = 0; i < n; ++i) {
        rx[i] = bodies[i].x; ry[i] = bodies[i].y; rz[i] = bodies[i].z; rm[i] = bodies[i].mass;
    }

    // Prepare worker contexts
    TD *td = (TD*)calloc(num_threads, sizeof(*td));
    pthread_t *ths = (pthread_t*)calloc(num_threads, sizeof(*ths));
    if (!td || !ths) {
        perror("calloc");
        free(ths); free(td);
        free(rm); free(rz); free(ry); free(rx);
        return -1;
    }

    Sync sync;
    pthread_mutex_init(&sync.m, NULL);
    pthread_cond_init(&sync.cv_start, NULL);
    pthread_cond_init(&sync.cv_done, NULL);
    sync.epoch = 0U;
    sync.remaining = 0;
    sync.stop = 0;

    // Partition bodies
    const size_t T = num_threads;
    const size_t chunk = (n + T - 1) / T;

    // Launch workers
    for (size_t t = 0; t < T; ++t) {
        const size_t start = t * chunk;
        const size_t end   = (start < n) ? ((start + chunk > n) ? n : start + chunk) : n;

        td[t].b = bodies;
        td[t].n = n;
        td[t].start = start;
        td[t].end   = end;
        td[t].steps = steps;
        td[t].dt    = dt;
        td[t].rx = rx; td[t].ry = ry; td[t].rz = rz; td[t].rm = rm;
        td[t].sync = &sync;
        td[t].tid = t; td[t].T = T;
#ifdef __linux__
        td[t].want_affinity = 1;
        td[t].core_id = (int)t;
#endif
        pthread_create(&ths[t], NULL, worker, &td[t]);
    }

    // Main loop
    double total_force_time  = 0.0;
    double total_update_time = 0.0;

    const double t0 = now_sec();

    for (size_t s = 0; s < steps; ++s) {
        // Publish the current snapshot pointers (arrays are stable; data updated below)
        for (size_t t = 0; t < T; ++t) {
            td[t].rx = rx; td[t].ry = ry; td[t].rz = rz; td[t].rm = rm;
        }

        // Start step: set epoch, set remaining, broadcast once
        const double tf0 = now_sec();
        pthread_mutex_lock(&sync.m);
        sync.remaining = T;
        ++sync.epoch;                   // workers will observe new epoch
        pthread_cond_broadcast(&sync.cv_start);
        pthread_mutex_unlock(&sync.m);

        // Wait until all workers finish forces
        pthread_mutex_lock(&sync.m);
        while (sync.remaining > 0)
            pthread_cond_wait(&sync.cv_done, &sync.m);
        pthread_mutex_unlock(&sync.m);

        total_force_time += (now_sec() - tf0);

        // Serial Euler update
        const double tu0 = now_sec();
#pragma GCC ivdep
        for (size_t i = 0; i < n; ++i) {
            const double inv_m = 1.0 / bodies[i].mass;
            const double vx = bodies[i].vx + bodies[i].fx * inv_m * dt;
            const double vy = bodies[i].vy + bodies[i].fy * inv_m * dt;
            const double vz = bodies[i].vz + bodies[i].fz * inv_m * dt;

            bodies[i].x += vx * dt; bodies[i].y += vy * dt; bodies[i].z += vz * dt;
            bodies[i].vx = vx;       bodies[i].vy = vy;       bodies[i].vz = vz;
        }
        total_update_time += (now_sec() - tu0);

        // Snapshot output (optional)
        if (snaps && stride && (s % stride) == 0) {
            const size_t frame = s / stride;
            float *dst = snaps->xyz + frame * n * 3u;
            for (size_t i = 0; i < n; ++i) {
                dst[i*3+0] = (float)bodies[i].x;
                dst[i*3+1] = (float)bodies[i].y;
                dst[i*3+2] = (float)bodies[i].z;
            }
        }

        // Build next step's read-only SoA from updated bodies
        if (s + 1 < steps) {
            for (size_t i = 0; i < n; ++i) {
                rx[i] = bodies[i].x; ry[i] = bodies[i].y; rz[i] = bodies[i].z; rm[i] = bodies[i].mass;
            }
        }
    }

    const double total_time = now_sec() - t0;

    // Clean shutdown
    pthread_mutex_lock(&sync.m);
    sync.stop = 1;
    ++sync.epoch; // ensure any waiters wake
    pthread_cond_broadcast(&sync.cv_start);
    pthread_mutex_unlock(&sync.m);

    for (size_t t = 0; t < T; ++t) pthread_join(ths[t], NULL);

    // Timing outs
    if (avg_force_time)  *avg_force_time  = total_force_time  / (double)steps;
    if (avg_update_time) *avg_update_time = total_update_time / (double)steps;
    if (total_time_out)  *total_time_out  = total_time;

    // Cleanup
    pthread_cond_destroy(&sync.cv_done);
    pthread_cond_destroy(&sync.cv_start);
    pthread_mutex_destroy(&sync.m);

    free(ths);
    free(td);

    free(rm);
    free(rz);
    free(ry);
    free(rx);

    return 0;
}
