#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "nbody.h"
#include "viewer.h"

typedef struct { double fx, fy, fz; } Force;

/* Initialise all bodies */
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

/* Compute gravitational forces for a given body index range */
static void compute_forces_range(const Body *bodies, size_t n, size_t start, size_t end, Force *acc)
{
    for (size_t i = start; i < end; ++i)
    {
        const double mi = bodies[i].mass;
        const double xi = bodies[i].x, yi = bodies[i].y, zi = bodies[i].z;

        for (size_t j = i + 1; j < n; ++j)
        {
            const double dx = bodies[j].x - xi;
            const double dy = bodies[j].y - yi;
            const double dz = bodies[j].z - zi;
            const double r2 = dx*dx + dy*dy + dz*dz + EPS2;
            const double invr  = 1.0 / sqrt(r2);
            const double invr3 = invr * invr * invr;
            const double s  = G * mi * bodies[j].mass * invr3;

            const double fx = s * dx;
            const double fy = s * dy;
            const double fz = s * dz;

            acc[i].fx += fx;  acc[i].fy += fy;  acc[i].fz += fz;
            acc[j].fx -= fx;  acc[j].fy -= fy;  acc[j].fz -= fz;
        }
    }
}

/* Sequential wrapper that uses one scratch buffer */
void compute_forces(Body *bodies, size_t n)
{
    Force *scratch = malloc(n * sizeof(Force));
    if (!scratch) { perror("malloc"); exit(EXIT_FAILURE); }
    memset(scratch, 0, n * sizeof(Force));

    compute_forces_range(bodies, n, 0, n, scratch);

    for (size_t i = 0; i < n; ++i)
    {
        bodies[i].fx = scratch[i].fx;
        bodies[i].fy = scratch[i].fy;
        bodies[i].fz = scratch[i].fz;
    }
    free(scratch);
}

/* Integrate positions and velocities using leapfrog */
void update_bodies(Body *bodies, size_t n, double dt)
{
    for (size_t i = 0; i < n; ++i)
    {
        const double m = bodies[i].mass;
        if (m <= 0.0) continue;
        bodies[i].vx += (bodies[i].fx / m) * (0.5 * dt);
        bodies[i].vy += (bodies[i].fy / m) * (0.5 * dt);
        bodies[i].vz += (bodies[i].fz / m) * (0.5 * dt);
    }

    for (size_t i = 0; i < n; ++i)
    {
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
        bodies[i].z += bodies[i].vz * dt;
    }

    compute_forces(bodies, n);

    for (size_t i = 0; i < n; ++i)
    {
        const double m = bodies[i].mass;
        if (m <= 0.0) continue;
        bodies[i].vx += (bodies[i].fx / m) * (0.5 * dt);
    }
}

/* Run one full timestep */
static inline void sim_step(Body *b, size_t n, double dt)
{
    compute_forces(b, n);
    update_bodies(b, n, dt);
}

/* Get seconds since CLOCK_MONOTONIC */
static inline double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        fprintf(stderr, "Usage: %s number_of_particles timesteps [stride]\n", argv[0]);
        return EXIT_FAILURE;
    }

    const size_t n = strtoull(argv[1], NULL, 10);
    const size_t steps = strtoull(argv[2], NULL, 10);
    const size_t stride = (argc >= 4) ? strtoull(argv[3], NULL, 10) : 10;
    const double dt = 1e-3;

    Body *bodies = calloc(n, sizeof(Body));
    if (!bodies) { perror("calloc"); return EXIT_FAILURE; }
    init_bodies(bodies, n);

    const size_t frames = (steps + stride - 1) / stride;
    Snapshots snaps = { .xyz = malloc(frames * n * 3 * sizeof(float)), .frames = frames, .n = n };
    if (!snaps.xyz) { perror("malloc"); free(bodies); return EXIT_FAILURE; }

    double *masses = malloc(n * sizeof(double));
    if (!masses) { perror("malloc"); free(snaps.xyz); free(bodies); return EXIT_FAILURE; }
    for (size_t i = 0; i < n; ++i) masses[i] = bodies[i].mass;

    /* --- Timing start --- */
    double start_time = now_sec();

    size_t f = 0;
    for (size_t s = 0; s < steps; ++s)
    {
        sim_step(bodies, n, dt);

        if ((s % stride) == 0)
        {
            float *dst = snaps.xyz + (size_t)f * n * 3u;
            for (size_t i = 0; i < n; ++i)
            {
                dst[i*3 + 0] = (float)bodies[i].x;
                dst[i*3 + 1] = (float)bodies[i].y;
                dst[i*3 + 2] = (float)bodies[i].z;
            }
            ++f;
        }
    }

    /* --- Timing end --- */
    double elapsed = now_sec() - start_time;
    printf("\n--- Simulation Complete ---\n");
    printf("Bodies: %zu\nSteps: %zu\nTime: %.4f seconds\n", n, steps, elapsed);
    printf("Average time per step: %.6f s\n", elapsed / (double)steps);

    viewer_play(&snaps, masses);

    free(masses);
    free(snaps.xyz);
    free(bodies);
    return EXIT_SUCCESS;
}
