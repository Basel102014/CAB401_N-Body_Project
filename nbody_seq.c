#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nbody.h"

/* Placeholder implementations */
void init_bodies(Body *bodies, size_t n) {
    (void)bodies;
    (void)n;
    /* TODO: initialize masses and positions */
}

void compute_forces(Body *bodies, size_t n) {
    (void)bodies;
    (void)n;
    /* TODO: compute gravitational forces */
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
