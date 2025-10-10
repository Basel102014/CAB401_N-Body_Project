#ifndef TEST_PRESETS_H
#define TEST_PRESETS_H

#include <stddef.h>
#include "nbody.h"
#include "cli_helpers.h"  // For TestPreset definition

/* -------------------------------------------------------------------------
 *  Test preset lookup
 * -------------------------------------------------------------------------
 *  Provides an interface for resolving named test configurations
 *  (e.g., "solar", "binary", "cluster") into preset parameter sets.
 * ------------------------------------------------------------------------- */

/**
 * @brief Find a predefined test preset by name.
 *
 * Looks up a preset (e.g., "solar", "cluster", "binary") and
 * returns a pointer to its TestPreset structure.
 *
 * @param name The name of the preset to look for.
 * @return Pointer to TestPreset if found, otherwise NULL.
 */
const TestPreset *find_preset(const char *name);

#endif /* TEST_PRESETS_H */
