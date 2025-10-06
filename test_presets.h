#ifndef TEST_PRESETS_H
#define TEST_PRESETS_H

#include <stddef.h>
#include "nbody.h"
#include "cli_helpers.h" // for TestPreset

// Returns a pointer to the preset if found, otherwise NULL
const TestPreset *find_preset(const char *name);

#endif
