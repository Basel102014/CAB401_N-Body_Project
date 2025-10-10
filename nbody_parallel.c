#define _POSIX_C_SOURCE 200809L
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#ifdef __linux__
#include <sched.h>
#endif
#include "nbody_parallel.h"

/* High-precision time */
static inline double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Optional: pin a thread to a specific core to improve cache locality */
static inline void maybe_pin_to_core(int want, int core_id) {
#ifdef __linux__
    if (!want) return;
    cpu_set_t set;
    CPU_ZERO(&set);
    if (core_id < 0) core_id = 0;
    CPU_SET((unsigned)core_id, &set);
    pthread_setaffinity_np(pthread_self(), sizeof(set), &set);
#else
    (void)want; (void)core_id;
#endif
}

/* A little cache-blocking on the inner j-loop helps a bit on mid-size N */
#ifndef NBLOCK
#define NBLOCK 64
#endif

/* Worker thread:
   - For each step:
     * Compute forces for i in [start,end), reading positions from b[0..n)
     * Write directly into b[i].fx/fy/fz (only this thread touches its i-range)
     * barrier #1: signal "forces done" and wait for main to update
     * barrier #2: wait until main finishes Euler update (positions ready for next step)
*/
static void *worker(void *arg)
{
    ThreadData *td = (ThreadData*)arg;

    Body *restrict b = td->b;
    const size_t n   = td->n;
    const size_t is  = td->start;
    const size_t ie  = td->end;
    const size_t steps = td->steps;
    pthread_barrier_t *bar = td->bar;

    maybe_pin_to_core(
#ifdef __linux__
        td->want_affinity
#else
        0
#endif
        ,
#ifdef __linux__
        td->core_id
#else
        0
#endif
    );

    for (size_t s = 0; s < steps; ++s) {

        /* Compute forces for my i-range (with simple j-blocking) */
        for (size_t i = is; i < ie; ++i) {
            const double mi = b[i].mass;
            const double xi = b[i].x, yi = b[i].y, zi = b[i].z;

            double fx = 0.0, fy = 0.0, fz = 0.0;

            for (size_t jb = 0; jb < n; jb += NBLOCK) {
                const size_t jend = (jb + NBLOCK < n) ? (jb + NBLOCK) : n;
                for (size_t j = jb; j < jend; ++j) {
                    if (j == i) continue;
                    const double dx = b[j].x - xi;
                    const double dy = b[j].y - yi;
                    const double dz = b[j].z - zi;
                    const double r2 = dx*dx + dy*dy + dz*dz + EPS2;
                    const double invr  = 1.0 / sqrt(r2);
                    const double invr3 = invr * invr * invr;
                    const double sF = G * mi * b[j].mass * invr3;
                    fx += sF * dx;
                    fy += sF * dy;
                    fz += sF * dz;
                }
            }

            /* Write into my own slot—no races */
            b[i].fx = fx; b[i].fy = fy; b[i].fz = fz;
        }

        /* Barrier #1: signal that all forces are ready; wait for main */
        pthread_barrier_wait(bar);

        /* Barrier #2: wait for main to complete the Euler update */
        pthread_barrier_wait(bar);
    }

    return NULL;
}

/* Public API: forces parallel, Euler update in main, barriers for sync */
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

    pthread_t *ths     = (pthread_t*)calloc(num_threads, sizeof(*ths));
    ThreadData *td     = (ThreadData*)calloc(num_threads, sizeof(*td));
    pthread_barrier_t bar;

    if (!ths || !td) { perror("calloc"); free(ths); free(td); return -1; }

    /* Main participates in the barrier, so count = workers + 1 */
    if (pthread_barrier_init(&bar, NULL, (unsigned)(num_threads + 1)) != 0) {
        perror("pthread_barrier_init");
        free(ths); free(td);
        return -1;
    }

    /* Partition i-range */
    const size_t chunk = (n + num_threads - 1) / num_threads;

    /* Launch workers */
    for (size_t t = 0; t < num_threads; ++t) {
        const size_t start = t * chunk;
        const size_t end   = (start < n) ? ((start + chunk > n) ? n : start + chunk) : n;

        td[t].b = bodies;
        td[t].n = n;
        td[t].start = start;
        td[t].end   = end;
        td[t].steps = steps;
        td[t].num_threads = num_threads;
        td[t].dt = dt;
        td[t].bar = &bar;
#ifdef __linux__
        td[t].want_affinity = 1;
        td[t].core_id = (int)t;   /* simple 1:1 pin */
#endif
        pthread_create(&ths[t], NULL, worker, &td[t]);
    }

    /* Main loop: barrier with workers, do Euler, optionally snapshot */
    double total_force_time  = 0.0;
    double total_update_time = 0.0;

    const double t_sim0 = now_sec();

    for (size_t s = 0; s < steps; ++s) {
        /* Wait for workers to finish forces; time this interval */
        const double tf0 = now_sec();
        pthread_barrier_wait(&bar);    /* matches workers' barrier #1 */
        total_force_time += (now_sec() - tf0);

        /* Sequential Euler update (tiny vs forces) */
        const double tu0 = now_sec();
        for (size_t i = 0; i < n; ++i) {
            const double inv_m = 1.0 / bodies[i].mass;
            bodies[i].vx += bodies[i].fx * inv_m * dt;
            bodies[i].vy += bodies[i].fy * inv_m * dt;
            bodies[i].vz += bodies[i].fz * inv_m * dt;
            bodies[i].x  += bodies[i].vx * dt;
            bodies[i].y  += bodies[i].vy * dt;
            bodies[i].z  += bodies[i].vz * dt;
        }
        total_update_time += (now_sec() - tu0);

        /* Save snapshot after the update (positions at t+dt) */
        if (snaps && (s % stride) == 0) {
            const size_t frame_index = s / stride;
            float *dst = snaps->xyz + frame_index * n * 3u;
            for (size_t i = 0; i < n; ++i) {
                dst[i*3 + 0] = (float)bodies[i].x;
                dst[i*3 + 1] = (float)bodies[i].y;
                dst[i*3 + 2] = (float)bodies[i].z;
            }
        }

        /* Let workers proceed to next step (they will read new positions) */
        pthread_barrier_wait(&bar);    /* matches workers' barrier #2 */
    }

    const double total_time = now_sec() - t_sim0;

    /* Join & destroy */
    for (size_t t = 0; t < num_threads; ++t) pthread_join(ths[t], NULL);
    pthread_barrier_destroy(&bar);
    free(td);
    free(ths);

    /* Return timings (no printing; main handles output/CSV) */
    if (avg_force_time)  *avg_force_time  = total_force_time  / (double)steps;
    if (avg_update_time) *avg_update_time = total_update_time / (double)steps;
    if (total_time_out)  *total_time_out  = total_time;

    return 0;
}
