#include <string.h>
#include <stdio.h>
#include "nbody_tests.h"
#include "test_presets.h"

// ---- Define all available test presets ----
const TestPreset PRESETS[] = {
    { "solar",   setup_solar_system, 9,     1000000, 10 },
    { "binary",  setup_binary_star,  2,      500000, 10 },
    { "trinary", setup_trinary_star, 3,     1000000, 10 },
    { "cluster", setup_random_cluster, 50,  1000000, 10 },
    { "lattice", setup_lattice_grid, 27,    1000000, 10 },
};

const size_t PRESET_COUNT = sizeof(PRESETS)/sizeof(PRESETS[0]);

// ---- Lookup function ----
const TestPreset *find_preset(const char *name) {
    for (size_t i = 0; i < sizeof(PRESETS) / sizeof(PRESETS[0]); ++i) {
        if (strcmp(PRESETS[i].name, name) == 0)
            return &PRESETS[i];
    }
    return NULL;
}
