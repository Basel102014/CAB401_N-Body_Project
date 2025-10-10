#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include "nbody_seq.h"
#include "nbody_parallel.h"

/* Existing now_sec() from your serial code */
static inline double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv)
{
    Options opt;
    if (parse_cli(argc, argv, &opt) != 0)
        return EXIT_FAILURE;

    const size_t n = opt.n;
    const size_t steps = opt.steps;
    const size_t stride = opt.stride;
    const bool save_csv = opt.csv;
    const bool no_view = opt.noview;
    const size_t threads = opt.threads;
    const double dt = 1e-3;
    BodySetupFn custom_setup = opt.setup;

    Body *bodies = calloc(n, sizeof(Body));
    if (!bodies)
    {
        perror("calloc bodies");
        return EXIT_FAILURE;
    }
    init_bodies(bodies, n, custom_setup);

    const size_t frames = (steps + stride - 1) / stride;
    Snapshots snaps = {
        .xyz = malloc(frames * n * 3 * sizeof(float)),
        .frames = frames,
        .n = n};
    if (!snaps.xyz)
    {
        perror("malloc snaps");
        free(bodies);
        return EXIT_FAILURE;
    }

    double *masses = malloc(n * sizeof(double));
    if (!masses)
    {
        perror("malloc masses");
        free(snaps.xyz);
        free(bodies);
        return EXIT_FAILURE;
    }
    for (size_t i = 0; i < n; ++i)
        masses[i] = bodies[i].mass;

    FILE *fp = NULL;
    if (save_csv)
    {
        fp = fopen("timing_results.csv", "w");
        if (fp)
            fprintf(fp, "step,bodies,threads,force_time,update_time,total_time\n");
        else
            perror("fopen csv");
    }

    /* --------- Sequential or Parallel Path --------- */
    if (threads > 1)
    {
        printf("Running in parallel mode with %zu threads\n", threads);
        double avg_force = 0.0, avg_update = 0.0, total_time = 0.0;
        run_nbody_parallel(bodies, n, steps, dt, threads,
                           &avg_force, &avg_update,
                           stride, &snaps, &total_time);

        if (fp)
            fprintf(fp, "0,%zu,%zu,%.9f,%.9f,%.9f,%.9f\n",
                    n, threads, avg_force, avg_update,
                    avg_force + avg_update, total_time);

        printf("\n--- Timing Results (Parallel) ---\n");
        printf("Bodies: %zu | Steps: %zu | Threads: %zu\n", n, steps, threads);
        printf("Avg Force Time : %.8f s\n", avg_force);
        printf("Avg Update Time: %.8f s\n", avg_update);
        printf("Avg Total Step : %.8f s\n", avg_force + avg_update);
        printf("Total Simulation Time: %.8f s\n", total_time);
    }
    else
    {
        printf("Running in sequential mode\n");
        double total_force = 0.0, total_update = 0.0, total_sim = 0.0;

        size_t f = 0;
        for (size_t s = 0; s < steps; ++s)
        {
            double force_t = 0.0, update_t = 0.0;
            double total_t = sim_step(bodies, n, dt, &force_t, &update_t);

            total_force += force_t;
            total_update += update_t;
            total_sim += total_t;

            if (fp)
                fprintf(fp, "%zu,%zu,%zu,%.9f,%.9f,%.9f\n", s, n, threads, force_t, update_t, total_t);

            if ((s % stride) == 0)
            {
                float *dst = snaps.xyz + f * n * 3u;
                for (size_t i = 0; i < n; ++i)
                {
                    dst[i * 3 + 0] = (float)bodies[i].x;
                    dst[i * 3 + 1] = (float)bodies[i].y;
                    dst[i * 3 + 2] = (float)bodies[i].z;
                }
                ++f;
            }
        }

        printf("\n--- Timing Results (Sequential) ---\n");
        printf("Bodies: %zu | Steps: %zu\n", n, steps);
        printf("Avg Force Time : %.8f s\n", total_force / steps);
        printf("Avg Update Time: %.8f s\n", total_update / steps);
        printf("Avg Total Step : %.8f s\n", total_sim / steps);
        printf("Total Simulation Time: %.8f s\n", total_sim);
    }

    if (fp)
        fclose(fp);

    if (!no_view)
        viewer_play(&snaps, masses);

    free(masses);
    free(snaps.xyz);
    free(bodies);
    return EXIT_SUCCESS;
}
