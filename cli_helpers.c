#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "nbody.h"
#include "nbody_tests.h"
#include "cli_helpers.h"
#include "test_presets.h"

void print_usage(const char *prog) {
    fprintf(stderr,
        "Usage:\n"
        "  %s --solar [--csv] [--stride=N] [--bodies=N] [--steps=N]\n"
        "  %s --bodies=N --steps=N [--stride=N] [--csv]\n"
        "  (legacy) %s <bodies> <steps> [stride] [name] [--csv]\n",
        prog, prog, prog);
}

int parse_cli(int argc, char **argv, Options *opt) {
    opt->n = 0;
    opt->steps = 0;
    opt->stride = 10;
    opt->csv = false;
    opt->setup = NULL;

    const TestPreset *chosen = NULL;

    for (int i = 1; i < argc; ++i) {
        const char *arg = argv[i];

        if (strncmp(arg, "--bodies=", 9) == 0) {
            opt->n = strtoull(arg + 9, NULL, 10);
            continue;
        }
        if (strncmp(arg, "--steps=", 8) == 0) {
            opt->steps = strtoull(arg + 8, NULL, 10);
            continue;
        }
        if (strncmp(arg, "--stride=", 9) == 0) {
            opt->stride = strtoull(arg + 9, NULL, 10);
            continue;
        }
        if (strcmp(arg, "--csv") == 0) {
            opt->csv = true;
            continue;
        }

        if (strncmp(arg, "--", 2) == 0) {
            const char *name = arg + 2;
            const TestPreset *p = find_preset(name);
            if (!p) {
                fprintf(stderr, "Unknown option/preset: %s\n", arg);
                return -1;
            }
            chosen = p;
            opt->setup = p->setup;
            continue;
        }

        // legacy positional parsing
        char *end = NULL;
        if (opt->n == 0) {
            opt->n = strtoull(arg, &end, 10);
            if (end && *end == '\0')
                continue;
        }
        if (opt->steps == 0) {
            opt->steps = strtoull(arg, &end, 10);
            if (end && *end == '\0')
                continue;
        }
        if (opt->stride == 10) {
            size_t s = strtoull(arg, &end, 10);
            if (end && *end == '\0') {
                opt->stride = s;
                continue;
            }
        }

        const TestPreset *p = find_preset(arg);
        if (p) {
            chosen = p;
            opt->setup = p->setup;
            continue;
        }

        fprintf(stderr, "Unknown argument: %s\n", arg);
        return -1;
    }

    if (chosen) {
        if (opt->n == 0) opt->n = chosen->def_bodies;
        if (opt->steps == 0) opt->steps = chosen->def_steps;
        if (chosen->def_stride > 0 && opt->stride == 10)
            opt->stride = chosen->def_stride;
    }

    if (opt->n == 0 || opt->steps == 0) {
        print_usage(argv[0]);
        return -1;
    }
    if (opt->stride == 0) {
        fprintf(stderr, "Error: --stride must be >= 1\n");
        return -1;
    }

    return 0;
}
