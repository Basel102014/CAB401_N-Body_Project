#define _POSIX_C_SOURCE 200809L
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include "nbody_parallel.h"

/* ------------------------ Timing ------------------------ */

static inline double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* --------------------- Helpers (local) ------------------ */

static inline void zero_acc_range(Force *acc, size_t start, size_t end) {
    for (size_t i = start; i < end; ++i)
        acc[i].fx = acc[i].fy = acc[i].fz = 0.0;
}

/* Race-free force accumulation:
   Each thread owns i in [start,end) and computes sum over j for acc[i].
   Shared b[j] reads are fine; each thread writes only its own acc[i]. */
static inline void compute_forces_range(const Body *b, size_t n,
                                        size_t start, size_t end,
                                        Force *acc)
{
    for (size_t i = start; i < end; ++i) {
        const double mi = b[i].mass;
        const double xi = b[i].x, yi = b[i].y, zi = b[i].z;

        /* Full j loop (skip self); softening avoids singularities. */
        for (size_t j = 0; j < n; ++j) {
            if (j == i) continue;
            const double dx = b[j].x - xi;
            const double dy = b[j].y - yi;
            const double dz = b[j].z - zi;
            const double r2   = dx*dx + dy*dy + dz*dz + EPS2;
            const double invr = 1.0 / sqrt(r2);
            const double invr3 = invr * invr * invr;
            const double s = G * mi * b[j].mass * invr3;

            acc[i].fx += s * dx;
            acc[i].fy += s * dy;
            acc[i].fz += s * dz;
        }
    }
}

/* --------------------- Worker Thread -------------------- */

static void *worker(void *arg)
{
    ThreadData *td = (ThreadData*)arg;
    Body  *b   = td->b;
    Force *acc = td->acc;
    const size_t n = td->n;
    const size_t start = td->start, end = td->end;
    const double dt = td->dt;

    for (size_t s = 0; s < td->steps; ++s) {
        /* Phase A: first half-kick + drift */
        if (td->force_t_out && start == 0) (void)now_sec(); // placeholder read to avoid unused warnings if you add per-phase timing later
        for (size_t i = start; i < end; ++i) {
            const double inv_m = 1.0 / b[i].mass;
            b[i].vx += 0.5 * dt * b[i].fx * inv_m;
            b[i].vy += 0.5 * dt * b[i].fy * inv_m;
            b[i].vz += 0.5 * dt * b[i].fz * inv_m;

            b[i].x  += b[i].vx * dt;
            b[i].y  += b[i].vy * dt;
            b[i].z  += b[i].vz * dt;
        }

        /* Barrier #1: ensure all positions are updated before force reads */
        pthread_barrier_wait(td->barrier);

        /* Phase B: compute forces at x(t+dt) */
        double tF0 = 0.0, tF1 = 0.0;
        if (td->force_t_out && start == 0) tF0 = now_sec();

        zero_acc_range(acc, start, end);
        compute_forces_range(b, n, start, end, acc);

        /* Barrier #2: ensure all acc[i] totals are complete */
        pthread_barrier_wait(td->barrier);

        if (td->force_t_out && start == 0) {
            tF1 = now_sec();
            td->force_t_out[s] = (tF1 - tF0);
        }

        /* Phase C: second half-kick (using new forces) */
        for (size_t i = start; i < end; ++i) {
            const double inv_m = 1.0 / b[i].mass;
            b[i].vx += 0.5 * dt * acc[i].fx * inv_m;
            b[i].vy += 0.5 * dt * acc[i].fy * inv_m;
            b[i].vz += 0.5 * dt * acc[i].fz * inv_m;

            /* Keep forces on Body for the next step’s first half-kick (viewer/debug) */
            b[i].fx = acc[i].fx;
            b[i].fy = acc[i].fy;
            b[i].fz = acc[i].fz;
        }

        /* Barrier #3: ensure all velocities are final before next step */
        pthread_barrier_wait(td->barrier);
    }

    return NULL;
}

/* ----------------- Public Entry Point (API) -------------- */

int run_nbody_parallel(Body *bodies,
                       size_t n, size_t steps, double dt,
                       size_t num_threads,
                       double *avg_force_time,   /* out, nullable */
                       double *avg_update_time)  /* out, nullable */
{
    int rc = 0;

    if (!bodies || n == 0 || steps == 0 || num_threads == 0) return -1;

    pthread_t   *ths = NULL;
    ThreadData  *td  = NULL;
    Force       *acc = NULL;
    double      *force_t_series = NULL;

    pthread_barrier_t barrier;
    int barrier_inited = 0;

    ths = (pthread_t*)  calloc(num_threads, sizeof(*ths));
    td  = (ThreadData*) calloc(num_threads, sizeof(*td));
    acc = (Force*)      calloc(n, sizeof(*acc));
    force_t_series = (double*)calloc(steps, sizeof(double));  /* thread 0 writes */

    if (!ths || !td || !acc || !force_t_series) { rc = -1; goto cleanup; }

    /* Initial force computation (so first half-kick uses correct a(t)) */
    zero_acc_range(acc, 0, n);
    compute_forces_range(bodies, n, 0, n, acc);
    for (size_t i = 0; i < n; ++i) {
        bodies[i].fx = acc[i].fx;
        bodies[i].fy = acc[i].fy;
        bodies[i].fz = acc[i].fz;
    }

    if (pthread_barrier_init(&barrier, NULL, (unsigned)num_threads) != 0) {
        perror("pthread_barrier_init");
        rc = -1; goto cleanup;
    }
    barrier_inited = 1;

    /* Partition */
    const size_t chunk = (n + num_threads - 1) / num_threads;

    /* Launch workers */
    const double t_sim_start = now_sec();
    for (size_t t = 0; t < num_threads; ++t) {
        const size_t start = t * chunk;
        const size_t end   = (start < n) ? ((start + chunk > n) ? n : start + chunk) : n;

        td[t].b = bodies;
        td[t].acc = acc;
        td[t].n = n;
        td[t].start = start;
        td[t].end = end;
        td[t].dt = dt;
        td[t].steps = steps;
        td[t].barrier = &barrier;
        td[t].force_t_out = (t == 0 ? force_t_series : NULL);

        if (pthread_create(&ths[t], NULL, worker, &td[t]) != 0) {
            perror("pthread_create");
            rc = -1;
            /* Try to join any already created threads */
            for (size_t j = 0; j < t; ++j) pthread_join(ths[j], NULL);
            goto cleanup;
        }
    }

    for (size_t t = 0; t < num_threads; ++t) pthread_join(ths[t], NULL);
    const double t_sim_end = now_sec();

    /* Averages */
    const double total_sim  = t_sim_end - t_sim_start;
    const double avg_total  = total_sim / (double)steps;

    double sum_force = 0.0;
    for (size_t s = 0; s < steps; ++s) sum_force += force_t_series[s];
    const double avg_force  = sum_force / (double)steps;
    const double avg_update = avg_total - avg_force;

    if (avg_force_time)  *avg_force_time  = avg_force;
    if (avg_update_time) *avg_update_time = avg_update;

    /* Optional: summary print (comment out if you want a quiet library call) */
    printf("\n--- Parallel Timing Summary ---\n");
    printf("Threads: %zu\n", num_threads);
    printf("Bodies : %zu\nSteps: %zu\n", n, steps);
    printf("Avg Force Time : %.8f s\n", avg_force);
    printf("Avg Update Time: %.8f s\n", avg_update);
    printf("Avg Total Step : %.8f s\n", avg_total);
    printf("Total Time     : %.8f s\n", total_sim);

cleanup:
    if (barrier_inited) pthread_barrier_destroy(&barrier);
    free(force_t_series);
    free(acc);
    free(td);
    free(ths);
    return rc;
}
