#include <stddef.h>
#include "nbody.h"

typedef struct {
    // packed [frame][i][xyz] -> length = frames * n * 3
    float  *xyz;
    size_t  frames;
    size_t  n;
} Snapshots;

// Interactive (does physics inside the loop)
void viewer_run(Body *bodies, size_t n, size_t steps);

// Playback-only (no physics): show precomputed positions
// masses: array of length n (used for size mapping)
void viewer_play(const Snapshots *snaps, const double *masses);