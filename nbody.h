#ifndef NBODY_H
#define NBODY_H

#include <stddef.h>

/* -------------------------------------------------------------------------
 *  Core N-body simulation types and constants
 * ------------------------------------------------------------------------- */

/* Gravitational constants and simulation parameters */
#define G       1e-3   // Gravitational constant (scaled)
#define EPS2    1e-2   // Softening factor to avoid singularities
#define DT      1e-3   // Default time step (Δt per substep)
#define M_PI    3.14159265358979323846 // Pi constant (for geometry)

/* -------------------------------------------------------------------------
 *  Data structures
 * ------------------------------------------------------------------------- */

/**
 * @brief Represents a single body in the simulation.
 *
 * Each body stores its mass, position, velocity, and force components.
 */
typedef struct
{
    double mass;           // Mass of the body
    double x, y, z;        // Position
    double vx, vy, vz;     // Velocity
    double fx, fy, fz;     // Force accumulator
} Body;

/**
 * @brief Function pointer type for user-defined setup routines.
 *
 * Custom initialisation functions can be passed to init_bodies()
 * to generate reproducible configurations (e.g., binary star, cluster, etc.).
 */
typedef void (*BodySetupFn)(Body *bodies, size_t n);

/* -------------------------------------------------------------------------
 *  Shared API
 * ------------------------------------------------------------------------- */

/**
 * @brief Initialise bodies with default or user-defined setup.
 *
 * @param bodies      Array of Body structs to initialise.
 * @param n           Number of bodies.
 * @param custom_init Optional user-provided setup function.
 */
void init_bodies(Body *bodies, size_t n, BodySetupFn custom_init);

#endif /* NBODY_H */
