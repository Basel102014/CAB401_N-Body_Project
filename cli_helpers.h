#ifndef CLI_HELPERS_H
#define CLI_HELPERS_H

#include <stddef.h>
#include <stdbool.h>
#include "nbody.h"  // for BodySetupFn

// Parsed command-line options
typedef struct {
    size_t n;
    size_t steps;
    size_t stride;
    bool csv;
    bool noview;
    BodySetupFn setup;
} Options;

// Test preset definition
typedef struct {
    const char *name;
    BodySetupFn setup;
    size_t def_bodies;
    size_t def_steps;
    size_t def_stride;
} TestPreset;

// CLI utility functions
void print_usage(const char *prog);
int parse_cli(int argc, char **argv, Options *opt);

#endif
