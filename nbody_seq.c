#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include "nbody.h"
#include "viewer.h"

typedef struct { double fx, fy, fz; } Force;

/* Return monotonic seconds (high precision) */
static inline double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Initialise all bodies */
#include <time.h>  // for time()

void init_bodies(Body *bodies, size_t n)
{
    srand((unsigned)time(NULL)); // seed random generator

    const double range = 10.0;     //  cube
    const double half = range / 2; // half width for centering

    for (size_t i = 0; i < n; ++i)
    {
        bodies[i].mass = 1.0;

        // Random positions in [-5.0, +5.0]
        bodies[i].x = ((double)rand() / RAND_MAX) * range - half;
        bodies[i].y = ((double)rand() / RAND_MAX) * range - half;
        bodies[i].z = ((double)rand() / RAND_MAX) * range - half;

        // Random small initial velocities
        bodies[i].vx = 0.0; 
        bodies[i].vy = 0.0;  
        bodies[i].vz = 0.0;  

        // Clear forces
        bodies[i].fx = bodies[i].fy = bodies[i].fz = 0.0;
    }
}

/* Compute gravitational forces for [start,end) */
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

/* Sequential wrapper (timed) */
double compute_forces(Body *bodies, size_t n)
{
    double t0 = now_sec();

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

    return now_sec() - t0;  // elapsed seconds for this force computation
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

/* Run one full timestep (timed) */
double sim_step(Body *b, size_t n, double dt, double *compute_time)
{
    double t0 = now_sec();
    *compute_time = compute_forces(b, n);
    update_bodies(b, n, dt);
    return now_sec() - t0;
}

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        fprintf(stderr, "Usage: %s number_of_particles timesteps [stride] [--csv]\n", argv[0]);
        return EXIT_FAILURE;
    }

    const size_t n = strtoull(argv[1], NULL, 10);
    const size_t steps = strtoull(argv[2], NULL, 10);
    const size_t stride = (argc >= 4 && argv[3][0] != '-') ? strtoull(argv[3], NULL, 10) : 10;
    bool save_csv = (argc > 3 && strcmp(argv[argc - 1], "--csv") == 0);
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

    double total_sim_time = 0.0, total_force_time = 0.0;
    double min_sim = 1e9, max_sim = 0.0;
    double min_force = 1e9, max_force = 0.0;

    FILE *fp = NULL;
    if (save_csv)
    {
        fp = fopen("timing_results.csv", "w");
        if (!fp) { perror("fopen csv"); save_csv = false; }
        else fprintf(fp, "step,sim_time,force_time\n");
    }

    size_t f = 0;
    for (size_t s = 0; s < steps; ++s)
    {
        double compute_t = 0.0;
        double sim_t = sim_step(bodies, n, dt, &compute_t);

        total_sim_time += sim_t;
        total_force_time += compute_t;
        if (sim_t < min_sim) min_sim = sim_t;
        if (sim_t > max_sim) max_sim = sim_t;
        if (compute_t < min_force) min_force = compute_t;
        if (compute_t > max_force) max_force = compute_t;

        if (save_csv && fp)
            fprintf(fp, "%zu,%.9f,%.9f\n", s, sim_t, compute_t);

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

    if (fp) fclose(fp);

    printf("\n--- Timing Results ---\n");
    printf("Bodies: %zu\nSteps: %zu\n", n, steps);
    printf("Average sim step time : %.8f s\n", total_sim_time / steps);
    printf("Average force time    : %.8f s\n", total_force_time / steps);
    printf("Min / Max sim step    : %.8f / %.8f s\n", min_sim, max_sim);
    printf("Min / Max force time  : %.8f / %.8f s\n", min_force, max_force);
    printf("Total sim time        : %.8f s\n", total_sim_time);

    if (save_csv)
        printf("Per-step timing saved to timing_results.csv\n");

    viewer_play(&snaps, masses);

    free(masses);
    free(snaps.xyz);
    free(bodies);
    return EXIT_SUCCESS;
}
