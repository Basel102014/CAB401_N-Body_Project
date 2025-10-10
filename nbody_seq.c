#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include "nbody_seq.h"

typedef struct
{
    double fx, fy, fz;
} Force;

/* High precision timer */
static inline double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Initialise bodies with random positions and small velocities */
void init_bodies(Body *bodies, size_t n, BodySetupFn custom_init)
{
    if (custom_init)
    {
        /* If user provided their own setup, use that instead */
        custom_init(bodies, n);
        return;
    }

    /* Default random setup */
    srand((unsigned)time(NULL));

    const double range = 10.0;
    const double half = range / 2.0;
    unsigned int seed = 42U;

    for (size_t i = 0; i < n; ++i)
    {
        bodies[i].mass = 10.0;
        bodies[i].x = ((double)rand_r(&seed) / RAND_MAX) * range - half;
        bodies[i].y = ((double)rand_r(&seed) / RAND_MAX) * range - half;
        bodies[i].z = ((double)rand_r(&seed) / RAND_MAX) * range - half;
        bodies[i].vx = ((double)rand_r(&seed) / RAND_MAX) * 0.01 - 0.005;
        bodies[i].vy = ((double)rand_r(&seed) / RAND_MAX) * 0.01 - 0.005;
        bodies[i].vz = ((double)rand_r(&seed) / RAND_MAX) * 0.01 - 0.005;
        bodies[i].fx = bodies[i].fy = bodies[i].fz = 0.0;
    }
}

/* Compute gravitational forces */
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
            const double r2 = dx * dx + dy * dy + dz * dz + EPS2;
            const double invr = 1.0 / sqrt(r2);
            const double invr3 = invr * invr * invr;
            const double s = G * mi * bodies[j].mass * invr3;

            const double fx = s * dx;
            const double fy = s * dy;
            const double fz = s * dz;

            acc[i].fx += fx;
            acc[i].fy += fy;
            acc[i].fz += fz;
            acc[j].fx -= fx;
            acc[j].fy -= fy;
            acc[j].fz -= fz;
        }
    }
}

double compute_forces(Body *bodies, size_t n)
{
    double t0 = now_sec();
    Force *scratch = calloc(n, sizeof(Force));
    if (!scratch)
    {
        perror("calloc");
        exit(EXIT_FAILURE);
    }

    compute_forces_range(bodies, n, 0, n, scratch);

    for (size_t i = 0; i < n; ++i)
    {
        bodies[i].fx = scratch[i].fx;
        bodies[i].fy = scratch[i].fy;
        bodies[i].fz = scratch[i].fz;
    }

    free(scratch);
    return now_sec() - t0;
}

/* Integrate velocities and positions */
double update_bodies(Body *b, size_t n, double dt)
{
    double t0 = now_sec();
    for (size_t i = 0; i < n; ++i)
    {
        const double inv_m = 1.0 / b[i].mass;
        b[i].vx += b[i].fx * inv_m * dt; // v(t+dt) = v(t) + a(t)*dt
        b[i].vy += b[i].fy * inv_m * dt;
        b[i].vz += b[i].fz * inv_m * dt;

        b[i].x += b[i].vx * dt; // x(t+dt) = x(t) + v(t+dt)*dt
        b[i].y += b[i].vy * dt;
        b[i].z += b[i].vz * dt;
    }
    return now_sec() - t0;
}

/* Run one simulation step */
double sim_step(Body *b, size_t n, double dt,
                double *force_time, double *update_time)
{
    double t0 = now_sec();

    // 1) forces from x(t)
    *force_time = compute_forces(b, n);

    // 2) Euler update using a(t)
    *update_time = update_bodies(b, n, dt);

    return now_sec() - t0; // total step time
}

// /* Main entry point */
// int main(int argc, char **argv)
// {
//     Options opt;
//     if (parse_cli(argc, argv, &opt) != 0)
//         return EXIT_FAILURE;

//     const size_t n = opt.n;
//     const size_t steps = opt.steps;
//     const size_t stride = opt.stride;
//     const bool save_csv = opt.csv;
//     const bool no_view = opt.noview;
//     const double dt = 1e-3;
//     BodySetupFn custom_setup = opt.setup;

//     Body *bodies = calloc(n, sizeof(Body));
//     if (!bodies)
//     {
//         perror("calloc");
//         return EXIT_FAILURE;
//     }

//     // Initialize bodies (custom setup or random)
//     init_bodies(bodies, n, custom_setup);

//     const size_t frames = (steps + stride - 1) / stride;
//     Snapshots snaps = {
//         .xyz = malloc(frames * n * 3 * sizeof(float)),
//         .frames = frames,
//         .n = n};
//     if (!snaps.xyz)
//     {
//         perror("malloc");
//         free(bodies);
//         return EXIT_FAILURE;
//     }

//     double *masses = malloc(n * sizeof(double));
//     if (!masses)
//     {
//         perror("malloc");
//         free(snaps.xyz);
//         free(bodies);
//         return EXIT_FAILURE;
//     }
//     for (size_t i = 0; i < n; ++i)
//         masses[i] = bodies[i].mass;

//     double total_sim = 0.0, total_force = 0.0, total_update = 0.0;
//     double min_force = 1e9, max_force = 0.0;
//     double min_update = 1e9, max_update = 0.0;
//     double min_total = 1e9, max_total = 0.0;

//     FILE *fp = NULL;
//     if (save_csv)
//     {
//         fp = fopen("timing_results.csv", "w");
//         if (fp)
//             fprintf(fp, "step,bodies,force_time,update_time,total_time\n");
//         else
//             perror("fopen csv");
//     }

//     size_t f = 0;
//     for (size_t s = 0; s < steps; ++s)
//     {
//         double force_t = 0.0, update_t = 0.0;
//         double total_t = sim_step(bodies, n, dt, &force_t, &update_t);

//         total_force += force_t;
//         total_update += update_t;
//         total_sim += total_t;

//         if (force_t < min_force)
//             min_force = force_t;
//         if (force_t > max_force)
//             max_force = force_t;
//         if (update_t < min_update)
//             min_update = update_t;
//         if (update_t > max_update)
//             max_update = update_t;
//         if (total_t < min_total)
//             min_total = total_t;
//         if (total_t > max_total)
//             max_total = total_t;

//         if (save_csv && fp)
//             fprintf(fp, "%zu,%zu,%.9f,%.9f,%.9f\n",
//                     s, n, force_t, update_t, total_t);

//         if ((s % stride) == 0)
//         {
//             float *dst = snaps.xyz + f * n * 3u;
//             for (size_t i = 0; i < n; ++i)
//             {
//                 dst[i * 3 + 0] = (float)bodies[i].x;
//                 dst[i * 3 + 1] = (float)bodies[i].y;
//                 dst[i * 3 + 2] = (float)bodies[i].z;
//             }
//             ++f;
//         }
//     }

//     if (fp)
//         fclose(fp);

//     // --- Final timing summary ---
//     printf("\n--- Timing Results ---\n");
//     printf("Bodies: %zu\nSteps: %zu\n", n, steps);
//     printf("Avg Force Time : %.8f s\n", total_force / steps);
//     printf("Avg Update Time: %.8f s\n", total_update / steps);
//     printf("Avg Total Step : %.8f s\n", total_sim / steps);
//     printf("Min/Max Force  : %.8f / %.8f s\n", min_force, max_force);
//     printf("Min/Max Update : %.8f / %.8f s\n", min_update, max_update);
//     printf("Min/Max Total  : %.8f / %.8f s\n", min_total, max_total);
//     printf("Total Time     : %.8f s\n", total_sim);

//     if (save_csv)
//         printf("Per-step timing saved to timing_results.csv\n");

//     if (!no_view)
//         viewer_play(&snaps, masses);

//     free(masses);
//     free(snaps.xyz);
//     free(bodies);
//     return EXIT_SUCCESS;
// }
