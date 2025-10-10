#ifndef NBODY_SEQ_H
#define NBODY_SEQ_H

#include <stddef.h>
#include "nbody.h"
#include "viewer.h"
#include "nbody_tests.h"
#include "cli_helpers.h"

/* -------------------------------------------------------------------------
 *  Sequential N-body Simulation (Reference Implementation)
 * -------------------------------------------------------------------------
 *  Provides single-threaded physics integration using the Velocity Verlet
 *  method for performance benchmarking and correctness comparison against
 *  the parallel version.
 * ------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
 *  High-precision timer
 * ------------------------------------------------------------------------- */

/**
 * @brief Get current high-precision time in seconds.
 *
 * Uses POSIX clock_gettime() with CLOCK_MONOTONIC for stable timing.
 *
 * @return Current time in seconds as a double.
 */
double now_sec(void);

/* -------------------------------------------------------------------------
 *  Core simulation routines
 * ------------------------------------------------------------------------- */

/**
 * @brief Compute pairwise gravitational forces on all bodies.
 *
 * @param bodies Array of Body structs (positions, masses, etc.).
 * @param n      Number of bodies.
 * @return Time taken for the force calculation in seconds.
 */
double compute_forces(Body *bodies, size_t n);

/**
 * @brief Integrate positions and velocities using Velocity Verlet.
 *
 * Performs:
 *   - First half-kick (update velocity by ½Δt)
 *   - Drift (update position by Δt)
 *   - Recalculate forces
 *   - Second half-kick (update velocity by ½Δt)
 *
 * @param b               Array of bodies.
 * @param n               Number of bodies.
 * @param dt              Time step (Δt).
 * @param force_time_out  Optional output for mid-step force computation time.
 * @return Time taken for the update phase (excluding force time).
 */
double update_bodies(Body *b, size_t n, double dt, double *force_time_out);

/**
 * @brief Execute one complete simulation step (update + force).
 *
 * @param b             Array of bodies.
 * @param n             Number of bodies.
 * @param dt            Time step (Δt).
 * @param force_time    Output: time spent computing forces.
 * @param update_time   Output: time spent updating velocities/positions.
 * @return Total time for the step (force + update).
 */
double sim_step(Body *b, size_t n, double dt,
                double *force_time, double *update_time);

#endif /* NBODY_SEQ_H */
