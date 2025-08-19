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
    // reset force accumulators
    for (size_t i = 0; i < n; ++i) {
        bodies[i].fx = 0.0;
        bodies[i].fy = 0.0;
        bodies[i].fz = 0.0;
    }

    // pairwise forces with Newton's 3rd law
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            double dx = bodies[j].x - bodies[i].x;
            double dy = bodies[j].y - bodies[i].y;
            double dz = bodies[j].z - bodies[i].z;

            // softened distance squared
            double r2 = dx*dx + dy*dy + dz*dz + EPS2;

            // inverse r and inverse r^3
            double inv_r  = 1.0 / sqrt(r2);
            double inv_r3 = inv_r * inv_r * inv_r;

            // scalar magnitude for the vector force
            double s = G * bodies[i].mass * bodies[j].mass * inv_r3;

            double fx = s * dx;
            double fy = s * dy;
            double fz = s * dz;

            // apply +F to i and -F to j
            bodies[i].fx += fx;
            bodies[i].fy += fy;
            bodies[i].fz += fz;

            bodies[j].fx -= fx;
            bodies[j].fy -= fy;
            bodies[j].fz -= fz;
        }
    }
}


void update_bodies(Body *bodies, size_t n, double dt) {
    (void)bodies;
    (void)n;
    (void)dt;
    /* TODO: update velocities and positions */
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
