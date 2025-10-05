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

/* Monotonic clock (high precision) */
static inline double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/*
 * setup_solar_system():
 * Creates a simplified but correctly proportioned Solar System.
 * Units:
 *   - Distance ≈ astronomical units (AU)
 *   - Velocity ≈ km/s, scaled for stability
 *   - Mass scaled to prevent extreme forces
 *
 * Stable with G = 1, dt ≈ 1e-3 and your current integration scheme.
 */

void setup_solar_system(Body *b, size_t n)
{
    if (n < 9) {
        fprintf(stderr, "Solar system setup requires at least 9 bodies (Sun + 8 planets).\n");
        exit(EXIT_FAILURE);
    }

    // ---- constants ----
    const double mass_sun = 1.989e30;
    const double MASS_SCALE = 1e-25;
    const double DIST_SCALE = 2.0; // 1 AU = 2 units

    // ---- planet data ----
    struct Planet {
        const char *name;
        double mass, a, e; // kg, AU, eccentricity
    } P[] = {
        {"Mercury", 3.301e23, 0.39, 0.2056},
        {"Venus",   4.867e24, 0.72, 0.0068},
        {"Earth",   5.972e24, 1.00, 0.0167},
        {"Mars",    6.417e23, 1.52, 0.0934},
        {"Jupiter", 1.898e27, 5.20, 0.0489},
        {"Saturn",  5.683e26, 9.58, 0.0565},
        {"Uranus",  8.681e25, 19.20, 0.0463},
        {"Neptune", 1.024e26, 30.05, 0.0097}
    };

    const size_t planet_count = sizeof(P) / sizeof(P[0]);

    // ---- Sun ----
    const double Ms = mass_sun * MASS_SCALE;
    b[0].mass = Ms;
    b[0].x = b[0].y = b[0].z = 0.0;
    b[0].vx = b[0].vy = b[0].vz = 0.0;
    b[0].fx = b[0].fy = b[0].fz = 0.0;

    // ---- planets ----
    for (size_t i = 0; i < planet_count; ++i) {
        size_t idx = i + 1;

        double a = P[i].a * DIST_SCALE; // semi-major axis (scaled)
        double e = P[i].e;
        double theta = (M_PI / 4.0) * i; // offset around orbit

        // compute current orbital radius (r = a(1 - e^2)/(1 + e cos θ))
        double r = a * (1 - e * e) / (1 + e * cos(theta));

        // orbital speed (vis-viva)
        double v = sqrt(G * Ms * (2.0 / r - 1.0 / a));

        // position (XZ plane)
        double x = r * cos(theta);
        double z = r * sin(theta);

        // velocity (tangent, perpendicular to radius)
        double vx = -v * sin(theta);
        double vz =  v * cos(theta);

        // mass scaled
        double m = P[i].mass * MASS_SCALE;

        b[idx].mass = m;
        b[idx].x = x;
        b[idx].y = 0.0;
        b[idx].z = z;
        b[idx].vx = vx;
        b[idx].vy = 0.0;
        b[idx].vz = vz;
        b[idx].fx = b[idx].fy = b[idx].fz = 0.0;
    }

    // ---- momentum correction (keep center of mass stable) ----
    double px = 0.0, py = 0.0, pz = 0.0;
    for (size_t i = 0; i <= planet_count; ++i) {
        px += b[i].mass * b[i].vx;
        py += b[i].mass * b[i].vy;
        pz += b[i].mass * b[i].vz;
    }
    b[0].vx -= px / Ms;
    b[0].vy -= py / Ms;
    b[0].vz -= pz / Ms;
}

/* Initialise bodies with random positions and small velocities */
void init_bodies(Body *bodies, size_t n, BodySetupFn custom_init)
{
    if (custom_init) {
        /* If user provided their own setup, use that instead */
        custom_init(bodies, n);
        return;
    }

    /* Default random setup */
    srand((unsigned)time(NULL));

    const double range = 10.0;
    const double half  = range / 2.0;

    for (size_t i = 0; i < n; ++i) {
        bodies[i].mass = 10.0;
        bodies[i].x = ((double)rand() / RAND_MAX) * range - half;
        bodies[i].y = ((double)rand() / RAND_MAX) * range - half;
        bodies[i].z = ((double)rand() / RAND_MAX) * range - half;
        bodies[i].vx = ((double)rand() / RAND_MAX) * 0.01 - 0.005;
        bodies[i].vy = ((double)rand() / RAND_MAX) * 0.01 - 0.005;
        bodies[i].vz = ((double)rand() / RAND_MAX) * 0.01 - 0.005;
        bodies[i].fx = bodies[i].fy = bodies[i].fz = 0.0;
    }
}

/* Compute gravitational forces for [start,end) */
static void compute_forces_range(const Body *bodies, size_t n, size_t start, size_t end, Force *acc) {
    for (size_t i = start; i < end; ++i) {
        const double mi = bodies[i].mass;
        const double xi = bodies[i].x, yi = bodies[i].y, zi = bodies[i].z;

        for (size_t j = i + 1; j < n; ++j) {
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
double compute_forces(Body *bodies, size_t n) {
    double t0 = now_sec();

    Force *scratch = calloc(n, sizeof(Force));
    if (!scratch) { perror("calloc"); exit(EXIT_FAILURE); }

    compute_forces_range(bodies, n, 0, n, scratch);

    for (size_t i = 0; i < n; ++i) {
        bodies[i].fx = scratch[i].fx;
        bodies[i].fy = scratch[i].fy;
        bodies[i].fz = scratch[i].fz;
    }

    free(scratch);
    return now_sec() - t0;
}

/* Integrate positions and velocities using leapfrog (timed) */
double update_bodies(Body *bodies, size_t n, double dt) {
    double t0 = now_sec();

    for (size_t i = 0; i < n; ++i) {
        const double m = bodies[i].mass;
        if (m <= 0.0) continue;
        bodies[i].vx += (bodies[i].fx / m) * (0.5 * dt);
        bodies[i].vy += (bodies[i].fy / m) * (0.5 * dt);
        bodies[i].vz += (bodies[i].fz / m) * (0.5 * dt);
    }

    for (size_t i = 0; i < n; ++i) {
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
        bodies[i].z += bodies[i].vz * dt;
    }

    compute_forces(bodies, n);

    for (size_t i = 0; i < n; ++i) {
        const double m = bodies[i].mass;
        if (m <= 0.0) continue;
        bodies[i].vx += (bodies[i].fx / m) * (0.5 * dt);
        bodies[i].vy += (bodies[i].fy / m) * (0.5 * dt);
        bodies[i].vz += (bodies[i].fz / m) * (0.5 * dt);
    }

    return now_sec() - t0;
}

/* Run one full timestep (timed for each phase) */
double sim_step(Body *b, size_t n, double dt,
                double *force_time, double *update_time)
{
    double t0 = now_sec();

    *force_time = compute_forces(b, n);
    double tf1 = now_sec();

    *update_time = update_bodies(b, n, dt);
    double total = now_sec() - t0;

    // sanity: recompute total manually if needed
    (void)tf1;
    return total;
}

/* Main simulation entry */
int main(int argc, char **argv)
{
    if (argc < 3) {
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
    init_bodies(bodies, n, setup_solar_system);

    const size_t frames = (steps + stride - 1) / stride;
    Snapshots snaps = { .xyz = malloc(frames * n * 3 * sizeof(float)), .frames = frames, .n = n };
    if (!snaps.xyz) { perror("malloc"); free(bodies); return EXIT_FAILURE; }

    double *masses = malloc(n * sizeof(double));
    if (!masses) { perror("malloc"); free(snaps.xyz); free(bodies); return EXIT_FAILURE; }
    for (size_t i = 0; i < n; ++i) masses[i] = bodies[i].mass;

    double total_sim = 0.0, total_force = 0.0, total_update = 0.0;
    double min_force = 1e9, max_force = 0.0;
    double min_update = 1e9, max_update = 0.0;
    double min_total = 1e9, max_total = 0.0;

    FILE *fp = NULL;
    if (save_csv) {
        fp = fopen("timing_results.csv", "w");
        if (!fp) { perror("fopen csv"); save_csv = false; }
        else fprintf(fp, "step,force_time,update_time,total_time\n");
    }

    size_t f = 0;
    for (size_t s = 0; s < steps; ++s)
    {
        double force_t = 0.0, update_t = 0.0;
        double total_t = sim_step(bodies, n, dt, &force_t, &update_t);

        total_force += force_t;
        total_update += update_t;
        total_sim += total_t;

        if (force_t < min_force) min_force = force_t;
        if (force_t > max_force) max_force = force_t;
        if (update_t < min_update) min_update = update_t;
        if (update_t > max_update) max_update = update_t;
        if (total_t < min_total) min_total = total_t;
        if (total_t > max_total) max_total = total_t;

        if (save_csv && fp)
            fprintf(fp, "%zu,%.9f,%.9f,%.9f\n", s, force_t, update_t, total_t);

        if ((s % stride) == 0) {
            float *dst = snaps.xyz + f * n * 3u;
            for (size_t i = 0; i < n; ++i) {
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
    printf("Avg Force Time : %.8f s\n", total_force / steps);
    printf("Avg Update Time: %.8f s\n", total_update / steps);
    printf("Avg Total Step : %.8f s\n", total_sim / steps);
    printf("Min/Max Force  : %.8f / %.8f s\n", min_force, max_force);
    printf("Min/Max Update : %.8f / %.8f s\n", min_update, max_update);
    printf("Min/Max Total  : %.8f / %.8f s\n", min_total, max_total);
    printf("Total Time     : %.8f s\n", total_sim);

    if (save_csv)
        printf("Per-step timing saved to timing_results.csv\n");

    viewer_play(&snaps, masses);

    free(masses);
    free(snaps.xyz);
    free(bodies);
    return EXIT_SUCCESS;
}
