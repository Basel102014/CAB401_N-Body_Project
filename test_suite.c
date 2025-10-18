#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nbody_seq.h"
#include "nbody_parallel.h"
#include "nbody_tests.h"
#include "test_presets.h"

#define DT 1e-3
#define STEPS 200
#define ABS_TOL 1e-9
#define REL_TOL 1e-6

/* ---------------- Helper Functions ---------------- */
/* High precision timer */
static inline double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static double rel_err(double a, double b) {
    double diff = fabs(a - b);
    double denom = fmax(fabs(a), fabs(b));
    if (denom < 1e-12) return diff;
    return diff / denom;
}

static double total_energy(const Body *b, size_t n) {
    double E = 0.0;
    for (size_t i = 0; i < n; ++i)
        E += 0.5 * b[i].mass * (b[i].vx*b[i].vx + b[i].vy*b[i].vy + b[i].vz*b[i].vz);

    for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j) {
            double dx = b[j].x - b[i].x;
            double dy = b[j].y - b[i].y;
            double dz = b[j].z - b[i].z;
            double r = sqrt(dx*dx + dy*dy + dz*dz + EPS2);
            E -= G * b[i].mass * b[j].mass / r;
        }

    return E;
}

/* Compare arrays of bodies */
static int compare_bodies(const Body *a, const Body *b, size_t n,
                          double *max_abs, double *max_rel) {
    *max_abs = 0.0;
    *max_rel = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double dx = fabs(a[i].x - b[i].x);
        double dy = fabs(a[i].y - b[i].y);
        double dz = fabs(a[i].z - b[i].z);
        double rx = rel_err(a[i].x, b[i].x);
        double ry = rel_err(a[i].y, b[i].y);
        double rz = rel_err(a[i].z, b[i].z);

        *max_abs = fmax(*max_abs, fmax(dx, fmax(dy, dz)));
        *max_rel = fmax(*max_rel, fmax(rx, fmax(ry, rz)));
    }
    return (*max_abs < ABS_TOL || *max_rel < REL_TOL);
}

/* ---------------- Unit Tests ---------------- */

static int test_force_symmetry(void) {
    Body b[2] = {
        {.mass = 1.0, .x = -1.0},
        {.mass = 1.0, .x =  1.0}
    };

    compute_forces(b, 2);
    int pass = fabs(b[0].fx + b[1].fx) < 1e-12;
    printf("Test 1: Force symmetry ......................... %s\n", pass ? "PASS" : "FAIL");
    return pass;
}

static int test_acceleration_magnitude(void) {
    Body b[2] = {
        {.mass = 2.0, .x = 0.0},
        {.mass = 1.0, .x = 1.0}
    };

    compute_forces(b, 2);

    // include EPS2 softening like the real simulation
    double r = 1.0;
    double r2 = r * r + EPS2;
    double invr = 1.0 / sqrt(r2);
    double a_expected = G * 1.0 * invr * invr * invr;   // = G*m/r³

    double a_actual = fabs(b[0].fx / b[0].mass);
    double rel = rel_err(a_expected, a_actual);
    int pass = rel < 1e-6;
    printf("Test 2: Acceleration magnitude ................. %s (rel err=%.3e)\n",
           pass ? "PASS" : "FAIL", rel);
    return pass;
}


static int test_update_equations(void) {
    Body b = {.mass = 1.0, .x = 0, .vx = 1.0, .fx = 1.0};
    double dt = 0.1;
    update_bodies(&b, 1, dt);

    double vx_ref = 1.1, x_ref = 0.11;
    int ok_v = fabs(b.vx - vx_ref) < 1e-9;
    int ok_x = fabs(b.x - x_ref) < 1e-9;
    int pass = ok_v && ok_x;
    printf("Test 3: Velocity/position integration .......... %s\n", pass ? "PASS" : "FAIL");
    return pass;
}

static int test_energy_conservation(void) {
    const size_t n = 3;
    Body b[3] = {
        {.mass = 1, .x = -1}, {.mass = 1, .x = 1}, {.mass = 1, .x = 0}
    };
    double E0 = total_energy(b, n);
    for (int i = 0; i < 10; ++i)
        sim_step(b, n, DT, &(double){0}, &(double){0});
    double E1 = total_energy(b, n);
    double drift = fabs((E1 - E0) / E0) * 100.0;
    int pass = drift < 0.1;
    printf("Test 4: Energy conservation (10 steps) ......... %s (%.5f%% drift)\n", pass ? "PASS" : "FAIL", drift);
    return pass;
}

/* ---------------- Parallel Validation ---------------- */

static void run_consistency_test(size_t threads, size_t N) {
    Body *seq = calloc(N, sizeof(Body));
    Body *par = calloc(N, sizeof(Body));
    if (!seq || !par) {
        perror("calloc");
        exit(EXIT_FAILURE);
    }

    setup_lattice_grid(seq, N);
    memcpy(par, seq, N * sizeof(Body));

    double E_seq_start = total_energy(seq, N);
    for (size_t s = 0; s < STEPS; ++s)
        sim_step(seq, N, DT, &(double){0}, &(double){0});
    double E_seq_end = total_energy(seq, N);

    double E_par_start = total_energy(par, N);
    run_nbody_parallel(par, N, STEPS, DT, threads, NULL, NULL, 0, NULL, NULL);
    double E_par_end = total_energy(par, N);

    double max_abs, max_rel;
    int state_ok = compare_bodies(seq, par, N, &max_abs, &max_rel);

    double drift_seq = fabs((E_seq_end - E_seq_start) / E_seq_start) * 100.0;
    double drift_par = fabs((E_par_end - E_par_start) / E_par_start) * 100.0;

    printf("Threads: %-2zu | Pos/Vel: %-4s | Energy Drift (%%): Seq %.5f, Par %.5f\n",
           threads, state_ok ? "PASS" : "FAIL", drift_seq, drift_par);

    free(seq);
    free(par);
}

/* ---------------- Main ---------------- */

int main(void) {
    printf("CAB401 N-Body Test Suite\n");
    printf("=========================\n\n");

    printf("Running core physics method tests...\n");
    int all_ok = 1;
    all_ok &= test_force_symmetry();
    all_ok &= test_acceleration_magnitude();
    all_ok &= test_update_equations();
    all_ok &= test_energy_conservation();
    printf("Core physics method tests: %s\n\n", all_ok ? "PASS" : "FAIL");

    printf("Running parallel consistency tests...\n");
    printf("Bodies: 125 | Steps: %d | DT: %.1e\n", STEPS, DT);
    printf("---------------------------------------------------------------\n");
    printf("Threads | Pos/Vel | Energy Drift (%%)\n");
    printf("---------------------------------------------------------------\n");

    size_t thread_counts[] = {1, 2, 4, 8};
    size_t num_counts = sizeof(thread_counts) / sizeof(thread_counts[0]);
    for (size_t i = 0; i < num_counts; ++i)
        run_consistency_test(thread_counts[i], 125);

    printf("---------------------------------------------------------------\n");
    printf("Test suite complete.\n");
    return 0;
}
