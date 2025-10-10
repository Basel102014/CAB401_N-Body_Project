#define _POSIX_C_SOURCE 200809L
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include "nbody_parallel.h"

/* ------------------------ Timing ------------------------ */
static inline double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* --------------------- Worker Thread -------------------- */
static void *worker(void *arg)
{
    ThreadData *td = (ThreadData *)arg;
    Body *b = td->b;
    Force *acc = td->acc;
    size_t n = td->n, start = td->start, end = td->end;
    WorkSync *sync = td->sync;

    for (;;)
    {
        /* Wait for signal from main */
        pthread_mutex_lock(&sync->mutex);
        while (!sync->step_ready && !sync->stop)
            pthread_cond_wait(&sync->cond_start, &sync->mutex);
        if (sync->stop)
        {
            pthread_mutex_unlock(&sync->mutex);
            break;
        }
        pthread_mutex_unlock(&sync->mutex);

        /* Compute forces for [start, end) */
        for (size_t i = start; i < end; ++i)
        {
            const double mi = b[i].mass;
            const double xi = b[i].x, yi = b[i].y, zi = b[i].z;
            double fx = 0.0, fy = 0.0, fz = 0.0;

            for (size_t j = 0; j < n; ++j)
            {
                if (j == i)
                    continue;
                const double dx = b[j].x - xi;
                const double dy = b[j].y - yi;
                const double dz = b[j].z - zi;
                const double r2 = dx * dx + dy * dy + dz * dz + EPS2;
                const double invr = 1.0 / sqrt(r2);
                const double invr3 = invr * invr * invr;
                const double sF = G * mi * b[j].mass * invr3;
                fx += sF * dx;
                fy += sF * dy;
                fz += sF * dz;
            }

            acc[i].fx = fx;
            acc[i].fy = fy;
            acc[i].fz = fz;
        }

        /* Notify main thread */
        pthread_mutex_lock(&sync->mutex);
        if (++sync->workers_done == td->num_threads)
        {
            sync->step_ready = 0;
            pthread_cond_signal(&sync->cond_done);
        }
        pthread_mutex_unlock(&sync->mutex);
    }

    return NULL;
}

/* ----------------- Public Entry Point -------------------- */
int run_nbody_parallel(Body *bodies,
                       size_t n,
                       size_t steps,
                       double dt,
                       size_t num_threads,
                       double *avg_force_time,
                       double *avg_update_time,
                       size_t stride,
                       Snapshots *snaps,
                       double *total_time_out)
{
    if (!bodies || n == 0 || steps == 0 || num_threads == 0)
        return -1;

    pthread_t *threads = calloc(num_threads, sizeof(*threads));
    ThreadData *td = calloc(num_threads, sizeof(*td));
    Force *acc = calloc(n, sizeof(*acc));
    WorkSync sync;

    if (!threads || !td || !acc)
    {
        perror("calloc");
        return -1;
    }

    pthread_mutex_init(&sync.mutex, NULL);
    pthread_cond_init(&sync.cond_start, NULL);
    pthread_cond_init(&sync.cond_done, NULL);
    sync.step_ready = 0;
    sync.workers_done = 0;
    sync.stop = 0;

    /* Partition bodies */
    const size_t chunk = (n + num_threads - 1) / num_threads;

    for (size_t t = 0; t < num_threads; ++t)
    {
        size_t start = t * chunk;
        size_t end = (start + chunk > n) ? n : start + chunk;
        td[t] = (ThreadData){
            .b = bodies,
            .acc = acc,
            .n = n,
            .start = start,
            .end = end,
            .num_threads = num_threads,
            .sync = &sync};
        pthread_create(&threads[t], NULL, worker, &td[t]);
    }

    /* Simulation loop */
    double total_force_time = 0.0, total_update_time = 0.0;
    const double t0 = now_sec();

    for (size_t s = 0; s < steps; ++s)
    {
        const double tf0 = now_sec();

        pthread_mutex_lock(&sync.mutex);
        sync.step_ready = 1;
        sync.workers_done = 0;
        pthread_cond_broadcast(&sync.cond_start);
        pthread_mutex_unlock(&sync.mutex);

        pthread_mutex_lock(&sync.mutex);
        while (sync.workers_done < num_threads)
            pthread_cond_wait(&sync.cond_done, &sync.mutex);
        pthread_mutex_unlock(&sync.mutex);

        total_force_time += (now_sec() - tf0);

        /* Copy acc → bodies */
        for (size_t i = 0; i < n; ++i)
        {
            bodies[i].fx = acc[i].fx;
            bodies[i].fy = acc[i].fy;
            bodies[i].fz = acc[i].fz;
        }

        /* Euler update */
        const double tu0 = now_sec();
        for (size_t i = 0; i < n; ++i)
        {
            const double inv_m = 1.0 / bodies[i].mass;
            bodies[i].vx += bodies[i].fx * inv_m * dt;
            bodies[i].vy += bodies[i].fy * inv_m * dt;
            bodies[i].vz += bodies[i].fz * inv_m * dt;
            bodies[i].x += bodies[i].vx * dt;
            bodies[i].y += bodies[i].vy * dt;
            bodies[i].z += bodies[i].vz * dt;
        }
        total_update_time += (now_sec() - tu0);

        /* Save snapshot */
        if (snaps && (s % stride) == 0)
        {
            size_t frame_index = s / stride;
            float *dst = snaps->xyz + frame_index * n * 3u;
            for (size_t i = 0; i < n; ++i)
            {
                dst[i * 3 + 0] = (float)bodies[i].x;
                dst[i * 3 + 1] = (float)bodies[i].y;
                dst[i * 3 + 2] = (float)bodies[i].z;
            }
        }
    }

    const double total_time = now_sec() - t0;

    /* Clean shutdown */
    pthread_mutex_lock(&sync.mutex);
    sync.stop = 1;
    sync.step_ready = 1;
    pthread_cond_broadcast(&sync.cond_start);
    pthread_mutex_unlock(&sync.mutex);

    for (size_t t = 0; t < num_threads; ++t)
        pthread_join(threads[t], NULL);

    /* Output results */
    if (avg_force_time)
        *avg_force_time = total_force_time / (double)steps;
    if (avg_update_time)
        *avg_update_time = total_update_time / (double)steps;
    if (total_time_out)
        *total_time_out = total_time;

    pthread_mutex_destroy(&sync.mutex);
    pthread_cond_destroy(&sync.cond_start);
    pthread_cond_destroy(&sync.cond_done);
    free(acc);
    free(td);
    free(threads);
    return 0;
}
