#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nbody.h"

void init_bodies(Body *bodies, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        bodies[i].mass = 1.0; // Placeholder mass
        bodies[i].x = (double)i; // Placeholder position
        bodies[i].y = 0.0;
        bodies[i].z = 0.0;
        bodies[i].vx = 0.0;
        bodies[i].vy = 0.0;
        bodies[i].vz = 0.0;
        bodies[i].fx = 0.0;
        bodies[i].fy = 0.0;
        bodies[i].fz = 0.0;
    }
    printf("Initialized %zu bodies.\n", n);
    fflush(stdout);
}

void compute_forces(Body *bodies, size_t n) {
    // Pairwise interactions (upper triangle), apply equal & opposite forces
    for (size_t i = 0; i < n; ++i) {
        const double mix = bodies[i].mass;
        const double xi  = bodies[i].x, yi = bodies[i].y, zi = bodies[i].z;

        for (size_t j = i + 1; j < n; ++j) {
            const double dx = bodies[j].x - xi;
            const double dy = bodies[j].y - yi;
            const double dz = bodies[j].z - zi;

            // Softened distance^2 (avoids singularities / huge forces)
            const double r2   = dx*dx + dy*dy + dz*dz + EPS2;
            const double invr = 1.0 / sqrt(r2);
            const double invr3 = invr * invr * invr;

            // Scalar for vector force: G * m_i * m_j / r^3
            const double s = G * mix * bodies[j].mass * invr3;

            const double fx = s * dx;
            const double fy = s * dy;
            const double fz = s * dz;

            bodies[i].fx +=  fx; bodies[i].fy +=  fy; bodies[i].fz +=  fz;
            bodies[j].fx -=  fx; bodies[j].fy -=  fy; bodies[j].fz -=  fz;
        }
    }
    printf("Computed forces for %zu bodies.\n", n);
    fflush(stdout);
}


void update_bodies(Body *bodies, size_t n, double dt) {
    for (size_t i = 0; i < n; ++i) {
        const double m = bodies[i].mass;
        if (m <= 0.0) continue;            // guard against invalid mass

        // a = F / m
        const double ax = bodies[i].fx / m;
        const double ay = bodies[i].fy / m;
        const double az = bodies[i].fz / m;

        // semi-implicit Euler: v(t+dt) = v(t) + a*dt
        bodies[i].vx += ax * dt;
        bodies[i].vy += ay * dt;
        bodies[i].vz += az * dt;

        // then x(t+dt) = x(t) + v(t+dt)*dt
        bodies[i].x  += bodies[i].vx * dt;
        bodies[i].y  += bodies[i].vy * dt;
        bodies[i].z  += bodies[i].vz * dt;
    }

    for (size_t i = 0; i < n; ++i) { bodies[i].fx = bodies[i].fy = bodies[i].fz = 0.0; }
    printf("Updated positions and velocities for %zu bodies.\n", n);
    fflush(stdout);
}


int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s number_of_particles timesteps\n", argv[0]);
        return EXIT_FAILURE;
    }

    size_t n = (size_t)atoi(argv[1]);
    size_t steps = (size_t)atoi(argv[2]);
    Body *bodies = calloc(n, sizeof(Body));
    if (!bodies) {
        perror("calloc");
        return EXIT_FAILURE;
    }

    init_bodies(bodies, n);
    for (size_t step = 0; step < steps; ++step) {
        compute_forces(bodies, n);
        update_bodies(bodies, n, 1.0); /* dt = 1.0 as placeholder */
    }

    printf("Sequential simulation placeholder complete.\n");
    free(bodies);
    return EXIT_SUCCESS;
}
