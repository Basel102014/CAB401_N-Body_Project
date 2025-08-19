#ifndef NBODY_H
#define NBODY_H

#include <stddef.h>

#define G     1.0   // Gravitational constant      
#define EPS2  1e-2  // Softening factor to avoid singularities    

typedef struct {
    double mass;
    double x, y, z;
    double vx, vy, vz;
    double fx, fy, fz;
} Body;

/* Initialize bodies with values */
void init_bodies(Body *bodies, size_t n);

/* Compute forces between bodies */
void compute_forces(Body *bodies, size_t n);

/* Update positions/velocities based on forces */
void update_bodies(Body *bodies, size_t n, double dt);

#endif /* NBODY_H */
