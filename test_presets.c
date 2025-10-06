#include <string.h>
#include <stdio.h>
#include "nbody_tests.h"
#include "test_presets.h"

// ---- Define all available test presets ----
static const TestPreset PRESETS[] = {
    { "solar", setup_solar_system, 9, 1000000, 10 },
    // Add future tests here:
    // { "binary", setup_binary_star, 2, 100000, 1 },
};

// ---- Lookup function ----
const TestPreset *find_preset(const char *name) {
    for (size_t i = 0; i < sizeof(PRESETS) / sizeof(PRESETS[0]); ++i) {
        if (strcmp(PRESETS[i].name, name) == 0)
            return &PRESETS[i];
    }
    return NULL;
}
