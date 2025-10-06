#ifndef NBODY_TESTS_H
#define NBODY_TESTS_H

#include "nbody.h"

void setup_solar_system(Body *b, size_t n);
void setup_random_cluster(Body *b, size_t n);
void setup_binary_star(Body *b, size_t n);
void setup_trinary_star(Body *b, size_t n);
void setup_lattice_grid(Body *b, size_t n); 

#endif
