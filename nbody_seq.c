#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nbody.h"
#include "viewer.h" // <- the new module

void init_bodies(Body *bodies, size_t n)
{
    for (size_t i = 0; i < n; ++i)
    {
        bodies[i].mass = 1.0;
        bodies[i].x = (double)i - (double)(n - 1) * 0.5;
        bodies[i].y = 0.0;
        bodies[i].z = 0.0;
        bodies[i].vx = bodies[i].vy = bodies[i].vz = 0.0;
        bodies[i].fx = bodies[i].fy = bodies[i].fz = 0.0;
    }
}

void compute_forces(Body *bodies, size_t n)
{
    for (size_t i = 0; i < n; ++i)
    {
        const double mi = bodies[i].mass;
        const double xi = bodies[i].x, yi = bodies[i].y, zi = bodies[i].z;
        for (size_t j = i + 1; j < n; ++j)
        {
            const double dx = bodies[j].x - xi;
            const double dy = bodies[j].y - yi;
            const double dz = bodies[j].z - zi;
            const double r2 = dx * dx + dy * dy + dz * dz + EPS2;
            const double invr = 1.0 / sqrt(r2);
            const double invr3 = invr * invr * invr;
            const double s = G * mi * bodies[j].mass * invr3;
            const double fx = s * dx, fy = s * dy, fz = s * dz;
            bodies[i].fx += fx;
            bodies[i].fy += fy;
            bodies[i].fz += fz;
            bodies[j].fx -= fx;
            bodies[j].fy -= fy;
            bodies[j].fz -= fz;
        }
    }
}

void update_bodies(Body *bodies, size_t n, double dt)
{
    for (size_t i = 0; i < n; ++i)
    {
        const double m = bodies[i].mass;
        if (m <= 0.0)
            continue;
        const double ax = bodies[i].fx / m;
        const double ay = bodies[i].fy / m;
        const double az = bodies[i].fz / m;
        bodies[i].vx += ax * dt;
        bodies[i].vy += ay * dt;
        bodies[i].vz += az * dt;
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
        bodies[i].z += bodies[i].vz * dt;
        bodies[i].fx = bodies[i].fy = bodies[i].fz = 0.0;
    }
}

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        fprintf(stderr, "Usage: %s number_of_particles timesteps\n", argv[0]);
        return EXIT_FAILURE;
    }

    size_t n = strtoull(argv[1], NULL, 10);
    size_t steps = strtoull(argv[2], NULL, 10); // 0 = run until close

    Body *bodies = calloc(n, sizeof(Body));
    if (!bodies)
    {
        perror("calloc");
        return EXIT_FAILURE;
    }

    init_bodies(bodies, n);

    // hand off to the viewer (which will call your physics each frame)
    viewer_run(bodies, n, steps);

    free(bodies);
    return EXIT_SUCCESS;
}
