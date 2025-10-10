#ifndef NBODY_PARALLEL_H
#define NBODY_PARALLEL_H

#include <pthread.h>
#include "nbody.h"
#include "viewer.h"

/* -------------------------------------------------------------------------
 *  Parallel N-body simulation (forces-only parallelism)
 *  Using persistent worker threads with condition-variable sync.
 * ------------------------------------------------------------------------- */

/* Per-body force accumulator */
typedef struct
{
    double fx, fy, fz;
} Force;

/* Synchronisation structure shared between main and workers */
typedef struct
{
    pthread_mutex_t mutex;
    pthread_cond_t cond_start;
    pthread_cond_t cond_done;
    int step_ready;       // 1 when a new step is ready to process
    size_t workers_done;  // number of threads that have finished current step
    int stop;             // signal to stop all workers
} WorkSync;

/* Per-thread context (persistent across all steps) */
typedef struct
{
    Body *b;             // shared array of bodies
    Force *acc;          // shared array of force accumulators
    size_t n;            // total number of bodies
    size_t start;        // inclusive start index for this thread
    size_t end;          // exclusive end index for this thread
    size_t num_threads;  // total number of worker threads
    WorkSync *sync;      // shared sync object between threads
} ThreadData;

/* -------------------------------------------------------------------------
 *  Public API
 * ------------------------------------------------------------------------- */

/**
 * @brief Run the parallel N-body simulation (forces only, persistent threads).
 *
 * This version parallelises only the force computation phase, keeping
 * velocity/position updates sequential for minimal synchronisation overhead.
 * Threads are persistent and synchronised each step using condition variables.
 *
 * @param bodies           Shared array of bodies.
 * @param n                Number of bodies.
 * @param steps            Number of simulation steps.
 * @param dt               Time step size.
 * @param num_threads      Number of worker threads.
 * @param avg_force_time   Output: average per-step force computation time (nullable).
 * @param avg_update_time  Output: average per-step update time (nullable).
 * @param stride           Step interval between saved frames (e.g., every Nth step).
 * @param snaps            Optional pointer to Snapshots structure for viewer output.
 *                         If non-null, frames will be written every `stride` steps.
 *
 * @return 0 on success, non-zero on failure.
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
