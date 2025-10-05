#ifndef NBODY_H
#define NBODY_H

#include <stddef.h>

#define G     1e-3  // Gravitational constant      
#define EPS2  1e-2  // Softening factor to avoid singularities    
#define DT    1e-3  // physics dt per substep
#define M_PI 3.14159265358979323846 // Pi constant

typedef struct {
    double mass;
    double x, y, z;
    double vx, vy, vz;
    double fx, fy, fz;
} Body;

/* Type for a user-defined setup callback */
typedef void (*BodySetupFn)(Body *bodies, size_t n);

/* Initialize bodies with values */
void init_bodies(Body *bodies, size_t n, BodySetupFn custom_init);

/* Compute forces between bodies */
double compute_forces(Body *bodies, size_t n);

/* Update positions/velocities based on forces */
double update_bodies(Body *bodies, size_t n, double dt);

#endif /* NBODY_H */
