#ifndef NBODY_TESTS_H
#define NBODY_TESTS_H

#include "nbody.h"

/* -------------------------------------------------------------------------
 *  Predefined N-body system setups for testing and demonstration
 * -------------------------------------------------------------------------
 *  Each setup function populates the given Body array with positions,
 *  velocities, and masses for a specific configuration.
 *
 *  These can be passed to init_bodies() via a BodySetupFn callback to
 *  generate reproducible test systems for simulation and benchmarking.
 * ------------------------------------------------------------------------- */

/**
 * @brief Initialise a simple model of the Solar System.
 *
 * Populates Sun + 8 planets with approximate masses and orbital velocities.
 * Requires at least 9 bodies.
 *
 * @param b Array of Body structs.
 * @param n Number of bodies (must be ≥ 9).
 */
void setup_solar_system(Body *b, size_t n);

/**
 * @brief Initialise a randomised gravitational cluster ("mini galaxy").
 *
 * One heavy core at the centre, with lighter orbiting bodies around it.
 * Requires at least 2 bodies.
 *
 * @param b Array of Body structs.
 * @param n Number of bodies (must be ≥ 2).
 */
void setup_random_cluster(Body *b, size_t n);

/**
 * @brief Initialise a two-body binary star system.
 *
 * Creates two stars orbiting a common centre of mass in circular motion.
 * Requires at least 2 bodies.
 *
 * @param b Array of Body structs.
 * @param n Number of bodies (must be ≥ 2).
 */
void setup_binary_star(Body *b, size_t n);

/**
 * @brief Initialise a three-body trinary star system.
 *
 * Three equal-mass stars arranged at 120° intervals in a rotating triangle.
 * Requires at least 3 bodies.
 *
 * @param b Array of Body structs.
 * @param n Number of bodies (must be ≥ 3).
 */
void setup_trinary_star(Body *b, size_t n);

/**
 * @brief Initialise bodies arranged in a uniform cubic lattice grid.
 *
 * Requires a perfect cube number of bodies (e.g., 8, 27, 64).
 * All bodies start stationary and equally spaced.
 *
 * @param b Array of Body structs.
 * @param n Number of bodies (must be a perfect cube).
 */
void setup_lattice_grid(Body *b, size_t n);

#endif /* NBODY_TESTS_H */
