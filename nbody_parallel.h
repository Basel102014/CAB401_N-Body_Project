#ifndef NBODY_PARALLEL_H
#define NBODY_PARALLEL_H

#include <pthread.h>
#include "nbody.h"
#include "viewer.h"

/* Per-body force accumulator (kept for compatibility in case you still use it
 * elsewhere, but the parallel kernel writes directly into Body.fx/fy/fz) */
typedef struct {
    double fx, fy, fz;
} Force;

/* Per-thread context (persistent) */
typedef struct {
    Body *b;                  // shared bodies (read positions, write own fx/fy/fz)
    size_t n;                 // total bodies
    size_t start;             // inclusive i-start (owned by this thread)
    size_t end;               // exclusive i-end
    size_t steps;             // total steps
    size_t num_threads;       // worker count
    double dt;                // time step
    pthread_barrier_t *bar;   // barrier shared by all workers + main (count = T+1)
#ifdef __linux__
    int want_affinity;        // pin to core if nonzero
    int core_id;              // core index to pin to
#endif
} ThreadData;

/* Parallel N-Body (forces-only parallelism; sequential Euler update)
 *
 * - Writes timing outputs via avg_force_time / avg_update_time / total_time_out
 * - Saves snapshots every `stride` steps into `snaps` if provided
 * - No printing: keep output uniform with your sequential path
 */
int run_nbody_parallel(Body *bodies,
                       size_t n,
                       size_t steps,
                       double dt,
                       size_t num_threads,
                       double *avg_force_time,
                       double *avg_update_time,
                       size_t stride,
                       Snapshots *snaps,
                       double *total_time_out);

#endif /* NBODY_PARALLEL_H */
