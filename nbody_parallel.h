#ifndef NBODY_PARALLEL_H
#define NBODY_PARALLEL_H

#include <pthread.h>
#include "nbody.h"

/* -------------------------------------------------------------------------
 *  Parallel N-body simulation types and function prototypes
 * ------------------------------------------------------------------------- */

/* Per-body force accumulator */
typedef struct {
    double fx, fy, fz;
} Force;

/* Per-thread data context for parallel N-body simulation */
typedef struct {
    Body   *b;          // Shared array of bodies
    Force  *acc;        // Shared per-body accumulators (fx, fy, fz)
    size_t  n;          // Total number of bodies
    size_t  start;      // Inclusive start index for this thread
    size_t  end;        // Exclusive end index for this thread
    double  dt;         // Time step (Δt)
    size_t  steps;      // Total simulation steps

    pthread_barrier_t *barrier; // Synchronisation barrier shared by all threads

    // Optional: per-step timing (only thread 0 writes)
    double *force_t_out;   // Array length = steps
    double *update_t_out;  // Array length = steps (optional)
    double *total_t_out;   // Array length = steps (optional)
} ThreadData;

/* -------------------------------------------------------------------------
 *  Function prototypes
 * ------------------------------------------------------------------------- */

/**
 * @brief Run the parallel N-body simulation using PThreads.
 *
 * @param bodies          Shared array of bodies.
 * @param n               Number of bodies.
 * @param steps           Number of simulation steps.
 * @param dt              Time step size.
 * @param num_threads     Number of worker threads.
 * @param avg_force_time  Output pointer for average force computation time (nullable).
 * @param avg_update_time Output pointer for average update time (nullable).
 *
 * @return 0 on success, non-zero on failure.
 */
int run_nbody_parallel(Body *bodies,
                       size_t n,
                       size_t steps,
                       double dt,
                       size_t num_threads,
                       double *avg_force_time,
                       double *avg_update_time);

/**
 * @brief Worker thread entry function (used internally).
 * @param arg Pointer to ThreadData for this thread.
 */
void *worker(void *arg);

/**
 * @brief Compute forces for a range of bodies (thread-safe if ranges don’t overlap).
 */
void compute_forces_range(const Body *b, size_t n,
                          size_t start, size_t end,
                          Force *acc);

/**
 * @brief Zero-out force accumulators in a range.
 */
void zero_acc_range(Force *acc, size_t start, size_t end);

#endif /* NBODY_PARALLEL_H */
