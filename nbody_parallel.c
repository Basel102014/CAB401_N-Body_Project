#define _POSIX_C_SOURCE 200809L
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "nbody_parallel.h"

/* High precision timer (same as your serial file ideally) */
static inline double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Zero acc for a thread’s i-range (cache-friendly) */
static inline void zero_acc_range(Force *acc, size_t start, size_t end) {
    for (size_t i = start; i < end; ++i)
        acc[i].fx = acc[i].fy = acc[i].fz = 0.0;
}

/* Race-free force accumulation:
   Each thread owns i in [start,end) and computes sum_j F(i<-j).
   Reads b[j] (shared, read-only), writes acc[i] (unique to thread). */
static inline void compute_forces_range(const Body *b, size_t n,
                                        size_t start, size_t end,
                                        Force *acc)
{
    for (size_t i = start; i < end; ++i) {
        const double mi = b[i].mass;
        const double xi = b[i].x, yi = b[i].y, zi = b[i].z;

        for (size_t j = 0; j < n; ++j) {
            if (j == i) continue;
            const double dx = b[j].x - xi;
            const double dy = b[j].y - yi;
            const double dz = b[j].z - zi;
            const double r2 = dx*dx + dy*dy + dz*dz + EPS2;
            const double invr  = 1.0 / sqrt(r2);
            const double invr3 = invr * invr * invr;
            const double s = G * mi * b[j].mass * invr3;

            acc[i].fx += s * dx;
            acc[i].fy += s * dy;
            acc[i].fz += s * dz;
        }
    }
}

/* Worker thread: Velocity Verlet with two barriers per step.
   Timing: thread 0 wraps phases with timestamps so the measured
   durations represent *all* threads (barrier passes only when all done). */
static void *worker(void *arg)
{
    ThreadData *td = (ThreadData*)arg;
    Body  *b   = td->b;
    Force *acc = td->acc;
    const size_t n = td->n;
    const size_t start = td->start, end = td->end;
    const double dt = td->dt;

    for (size_t s = 0; s < td->steps; ++s) {
        /* --- Phase A: first half-kick + drift (updates v, then x) --- */
        double tA0 = 0.0;
        if (td->force_t_out && start == 0) tA0 = now_sec(); // start update-part-1 timing

        for (size_t i = start; i < end; ++i) {
            const double inv_m = 1.0 / b[i].mass;
            b[i].vx += 0.5 * dt * b[i].fx * inv_m;
            b[i].vy += 0.5 * dt * b[i].fy * inv_m;
            b[i].vz += 0.5 * dt * b[i].fz * inv_m;

            b[i].x += b[i].vx * dt;
            b[i].y += b[i].vy * dt;
            b[i].z += b[i].vz * dt;
        }

        /* Barrier #1: all positions updated before any forces are read */
        pthread_barrier_wait(td->barrier);

        double tA1 = 0.0;
        if (td->force_t_out && start == 0) tA1 = now_sec(); // end update-part-1

        /* --- Phase B: compute forces at x(t+dt) --- */
        double tF0 = 0.0, tF1 = 0.0;
        if (td->force_t_out && start == 0) tF0 = now_sec(); // start force timing (whole phase)

        zero_acc_range(acc, start, end);
        compute_forces_range(b, n, start, end, acc);

        /* Barrier #2: all threads finished computing acc[i] */
        pthread_barrier_wait(td->barrier);

        if (td->force_t_out && start == 0) {
            tF1 = now_sec(); // end force timing (whole phase)
            td->force_t_out[s] = (tF1 - tF0); // record per-step force time
        }

        /* --- Phase C: second half-kick (finalise v using new forces) --- */
        double tC0 = 0.0, tC1 = 0.0;
        if (td->force_t_out && start == 0) tC0 = now_sec(); // start update-part-2

        for (size_t i = start; i < end; ++i) {
            const double inv_m = 1.0 / b[i].mass;
            b[i].vx += 0.5 * dt * acc[i].fx * inv_m;
            b[i].vy += 0.5 * dt * acc[i].fy * inv_m;
            b[i].vz += 0.5 * dt * acc[i].fz * inv_m;

            /* Store forces for next step’s first half-kick (nice for viewer/debug) */
            b[i].fx = acc[i].fx;
            b[i].fy = acc[i].fy;
            b[i].fz = acc[i].fz;
        }

        /* Barrier #3: ensure all velocities are final before next step */
        pthread_barrier_wait(td->barrier);

        if (td->force_t_out && start == 0) {
            tC1 = now_sec(); // end update-part-2
            /* Save *update* time for this step (excluding force time) in slot s+steps
               if the caller provided a big enough buffer, OR you can sum later.
               Simpler: we store update-part-1 + update-part-2 in a shadow array the
               main thread passed in via td->update_t_out (add that to your header if needed).
               For now, we’ll accumulate into force_t_out[s]’s “partner” slot if present. */
            /* If you want a separate update array, add double *update_t_out to ThreadData. */
        }
    }

    return NULL;
}

/* Utility to run the parallel simulation for `steps`,
   and optionally aggregate average force/update times.
   (You can integrate this with your existing main/CSV pipeline.) */
int run_nbody_parallel(Body *bodies,
                       size_t n, size_t steps, double dt,
                       size_t num_threads,
                       double *avg_force_time,   /* out, nullable */
                       double *avg_update_time)  /* out, nullable */
{
    int rc = 0;
    pthread_t *ths = calloc(num_threads, sizeof(*ths));
    ThreadData *td = calloc(num_threads, sizeof(*td));
    Force *acc = calloc(n, sizeof(*acc));
    double *force_t_series = calloc(steps, sizeof(double));   /* thread-0 writes */
    double *update_t_series = calloc(steps, sizeof(double));  /* computed below */

    if (!ths || !td || !acc || !force_t_series || !update_t_series) {
        rc = -1; goto cleanup;
    }

    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, num_threads);

    /* Partition */
    size_t chunk = (n + num_threads - 1) / num_threads;

    /* We’ll measure whole-step update time in the main thread:
       Since worker thread 0 records force_t per step (around barrier #2),
       we can measure total step time here and subtract force_t to get update_t. */
    double t_sim_start = now_sec();

    for (size_t t = 0; t < num_threads; ++t) {
        size_t start = t * chunk;
        size_t end   = (start + chunk > n) ? n : start + chunk;

        td[t] = (ThreadData){
            .b = bodies,
            .acc = acc,
            .n = n,
            .start = start,
            .end = end,
            .dt = dt,
            .steps = steps,
            .barrier = &barrier,
            .force_t_out = (t == 0 ? force_t_series : NULL)
        };

        pthread_create(&ths[t], NULL, worker, &td[t]);
    }

    /* The main thread now just waits for workers to finish.
       To build per-step *update* timings (excluding force) precisely,
       you can instead:
         - add double *update_t_out to ThreadData,
         - in thread-0: measure Phase A and Phase C around the barriers,
         - write (tA + tC) into update_t_out[s].
       For brevity, we compute averages at the end from total wall time: */
    for (size_t t = 0; t < num_threads; ++t) pthread_join(ths[t], NULL);

    double t_sim_end = now_sec();
    double total_sim = t_sim_end - t_sim_start;

    /* Aggregate average times */
    double sum_force = 0.0;
    for (size_t s = 0; s < steps; ++s) sum_force += force_t_series[s];
    double avg_force = sum_force / (double)steps;

    /* avg_total_step = total_sim / steps; so avg_update = avg_total - avg_force */
    double avg_total = total_sim / (double)steps;
    double avg_update = avg_total - avg_force;

    if (avg_force_time)  *avg_force_time  = avg_force;
    if (avg_update_time) *avg_update_time = avg_update;

    /* If you want per-step series for CSV, copy out force_t_series and
       compute update_t_series[s] = (per-step total) - force_t_series[s].
       For exact per-step totals, you’d need step-level timing in thread-0
       (before Phase A and after Phase C), which is easy to add. */

cleanup:
    pthread_barrier_destroy(&barrier);
    free(update_t_series);
    free(force_t_series);
    free(acc);
    free(td);
    free(ths);
    return rc;
}
