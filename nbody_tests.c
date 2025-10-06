#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nbody.h"
#include "nbody_tests.h"

#define MASS_SCALE 1e-25
#define DIST_SCALE 2.0

/* -------------------- Solar System -------------------- */
void setup_solar_system(Body *b, size_t n)
{
    if (n < 9)
    {
        fprintf(stderr, "Solar system setup requires at least 9 bodies (Sun + 8 planets).\n");
        exit(EXIT_FAILURE);
    }

    const double mass_sun = 1.989e30;
    const double Ms = mass_sun * MASS_SCALE;

    struct Planet
    {
        const char *name;
        double mass, a, e;
    } P[] = {
        {"Mercury", 3.301e23, 0.39, 0.2056},
        {"Venus", 4.867e24, 0.72, 0.0068},
        {"Earth", 5.972e24, 1.00, 0.0167},
        {"Mars", 6.417e23, 1.52, 0.0934},
        {"Jupiter", 1.898e27, 5.20, 0.0489},
        {"Saturn", 5.683e26, 9.58, 0.0565},
        {"Uranus", 8.681e25, 19.20, 0.0463},
        {"Neptune", 1.024e26, 30.05, 0.0097}};

    const size_t planet_count = sizeof(P) / sizeof(P[0]);

    b[0].mass = Ms;
    b[0].x = b[0].y = b[0].z = 0.0;
    b[0].vx = b[0].vy = b[0].vz = 0.0;

    for (size_t i = 0; i < planet_count; ++i)
    {
        size_t idx = i + 1;
        double a = P[i].a * DIST_SCALE;
        double e = P[i].e;
        double theta = (M_PI / 4.0) * i;

        double r = a * (1 - e * e) / (1 + e * cos(theta));
        double v = sqrt(G * Ms * (2.0 / r - 1.0 / a));

        b[idx].mass = P[i].mass * MASS_SCALE;
        b[idx].x = r * cos(theta);
        b[idx].y = 0.0;
        b[idx].z = r * sin(theta);
        b[idx].vx = -v * sin(theta);
        b[idx].vy = 0.0;
        b[idx].vz = v * cos(theta);
    }

    double px = 0, py = 0, pz = 0;
    for (size_t i = 0; i <= planet_count; ++i)
    {
        px += b[i].mass * b[i].vx;
        py += b[i].mass * b[i].vy;
        pz += b[i].mass * b[i].vz;
    }
    b[0].vx -= px / Ms;
    b[0].vy -= py / Ms;
    b[0].vz -= pz / Ms;
}

/* -------------------- Binary Star System -------------------- */
void setup_binary_star(Body *b, size_t n)
{
    if (n < 2)
    {
        fprintf(stderr, "Binary star system requires at least 2 bodies.\n");
        exit(EXIT_FAILURE);
    }

    const double m = 1.0e30 * MASS_SCALE;
    const double d = 4.0;
    const double v = sqrt(G * m / (d / 2.0));

    b[0].mass = m;
    b[0].x = -d / 2.0; b[0].y = 0.0; b[0].z = 0.0;
    b[0].vx = 0.0; b[0].vy = 0.0; b[0].vz = v;

    b[1].mass = m;
    b[1].x = d / 2.0; b[1].y = 0.0; b[1].z = 0.0;
    b[1].vx = 0.0; b[1].vy = 0.0; b[1].vz = -v;
}

/* -------------------- Trinary Star System -------------------- */
void setup_trinary_star(Body *b, size_t n)
{
    if (n < 3)
    {
        fprintf(stderr, "Trinary system requires at least 3 bodies.\n");
        exit(EXIT_FAILURE);
    }

    const double m = 1.0e30 * MASS_SCALE;
    const double r = 5.0;
    const double v = sqrt(G * m / (r * sqrt(3.0)));

    for (int i = 0; i < 3; ++i)
    {
        double ang = i * (2.0 * M_PI / 3.0);
        b[i].mass = m;
        b[i].x = r * cos(ang);
        b[i].y = 0.0;
        b[i].z = r * sin(ang);
        b[i].vx = -v * sin(ang);
        b[i].vy = 0.0;
        b[i].vz = v * cos(ang);
    }
}

/* -------------------- Random Cluster (Mini Galaxy) -------------------- */
void setup_random_cluster(Body *b, size_t n)
{
    if (n < 2)
    {
        fprintf(stderr, "Cluster test requires at least 2 bodies.\n");
        exit(EXIT_FAILURE);
    }

    srand(42); // deterministic

    const double core_mass = 1e31 * MASS_SCALE;
    const double radius = 20.0;

    b[0].mass = core_mass;
    b[0].x = b[0].y = b[0].z = 0.0;
    b[0].vx = b[0].vy = b[0].vz = 0.0;

    for (size_t i = 1; i < n; ++i)
    {
        double r = ((double)rand() / RAND_MAX) * radius;
        double theta = ((double)rand() / RAND_MAX) * 2 * M_PI;
        double phi = ((double)rand() / RAND_MAX) * M_PI;

        double x = r * sin(phi) * cos(theta);
        double y = r * sin(phi) * sin(theta);
        double z = r * cos(phi);

        double v = sqrt(G * core_mass / (r + 1.0));

        b[i].mass = 1e29 * MASS_SCALE;
        b[i].x = x; b[i].y = y; b[i].z = z;
        b[i].vx = -v * sin(theta);
        b[i].vy = v * cos(theta);
        b[i].vz = 0.0;
    }
}

/* -------------------- Lattice Grid -------------------- */
void setup_lattice_grid(Body *b, size_t n)
{
    const int grid = cbrt(n);
    if (grid * grid * grid != (int)n)
    {
        fprintf(stderr, "Lattice test requires perfect cube body count (e.g., 8, 27, 64).\n");
        exit(EXIT_FAILURE);
    }

    const double spacing = 4.0;
    const double mass = 1e29 * MASS_SCALE;

    int idx = 0;
    for (int x = 0; x < grid; ++x)
    {
        for (int y = 0; y < grid; ++y)
        {
            for (int z = 0; z < grid; ++z)
            {
                b[idx].mass = mass;
                b[idx].x = (x - grid / 2) * spacing;
                b[idx].y = (y - grid / 2) * spacing;
                b[idx].z = (z - grid / 2) * spacing;
                b[idx].vx = b[idx].vy = b[idx].vz = 0.0;
                idx++;
            }
        }
    }
}
