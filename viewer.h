#ifndef VIEWER_H
#define VIEWER_H

#include <stddef.h>
#include "nbody.h"

/* -------------------------------------------------------------------------
 *  N-body Viewer Interface
 * -------------------------------------------------------------------------
 *  Provides playback-only 3D visualisation using SDL2.
 *  The viewer renders precomputed positions from simulation snapshots.
 * ------------------------------------------------------------------------- */

/**
 * @brief Stores precomputed position data for playback.
 *
 * Positions are packed as [frame][body][xyz], with total length:
 *   frames * n * 3 floats
 */
typedef struct {
    float  *xyz;      // Interleaved position data
    size_t  frames;   // Number of recorded frames
    size_t  n;        // Number of bodies per frame
} Snapshots;

/**
 * @brief Playback-only renderer for simulation snapshots.
 *
 * Displays precomputed positions of bodies in 3D using SDL2.
 * No physics or integration occurs during playback.
 *
 * @param snaps   Pointer to Snapshots structure (must be valid).
 * @param masses  Array of body masses (used for size and brightness mapping).
 */
void viewer_play(const Snapshots *snaps, const double *masses);

#endif /* VIEWER_H */
