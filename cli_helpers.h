#ifndef CLI_HELPERS_H
#define CLI_HELPERS_H

#include <stddef.h>
#include <stdbool.h>
#include "nbody.h"  // for BodySetupFn

/* -------------------------------------------------------------------------
 *  Command-line interface helpers
 * -------------------------------------------------------------------------
 *  Provides structures and functions for parsing CLI arguments
 *  and configuring N-body simulation options.
 * ------------------------------------------------------------------------- */

/**
 * @brief Parsed command-line options.
 *
 * This structure holds all user-provided runtime settings
 * parsed from command-line arguments.
 */
typedef struct {
    size_t n;           // Number of bodies
    size_t steps;       // Total number of simulation steps
    size_t stride;      // Frame stride for snapshot saving
    bool csv;           // Enable CSV timing output
    bool noview;        // Disable visualisation mode
    BodySetupFn setup;  // Optional body setup function (preset)
} Options;

/**
 * @brief Predefined test preset configuration.
 *
 * Used to register preset simulation setups (e.g., solar system, cluster)
 * with default parameters.
 */
typedef struct {
    const char *name;   // Preset name (e.g., "solar", "cluster")
    BodySetupFn setup;  // Associated setup function
    size_t def_bodies;  // Default number of bodies
    size_t def_steps;   // Default number of steps
    size_t def_stride;  // Default snapshot stride
} TestPreset;

/* -------------------------------------------------------------------------
 *  CLI utility functions
 * ------------------------------------------------------------------------- */

/**
 * @brief Print program usage information.
 *
 * @param prog The name of the executable (argv[0]).
 */
void print_usage(const char *prog);

/**
 * @brief Parse command-line arguments into an Options struct.
 *
 * Supports options for:
 *  --n <value>        Number of bodies
 *  --steps <value>    Total steps
 *  --stride <value>   Snapshot stride
 *  --csv              Enable per-step timing output to CSV
 *  --noview           Disable visualisation (headless mode)
 *  --preset <name>    Select predefined test setup
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @param opt  Pointer to Options struct to populate.
 * @return 0 on success, non-zero on failure.
 */
int parse_cli(int argc, char **argv, Options *opt);

#endif /* CLI_HELPERS_H */
