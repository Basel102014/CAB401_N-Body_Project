#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nbody.h"
#include "nbody_tests.h"

void setup_solar_system(Body *b, size_t n)
{
    if (n < 9)
    {
        fprintf(stderr, "Solar system setup requires at least 9 bodies (Sun + 8 planets).\n");
        exit(EXIT_FAILURE);
    }

    const double mass_sun = 1.989e30;
    const double MASS_SCALE = 1e-25;
    const double DIST_SCALE = 2.0;

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
    const double Ms = mass_sun * MASS_SCALE;

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
