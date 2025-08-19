# N-Body Simulation (3D)

A simple, fast 3D N-Body gravitational simulation in C.
Includes a sequential solver. The sequential build has an SDL2 viewer for real-time visualization with a free-fly camera.

## Features
- 3D Newtonian gravity with softening (EPS2)
- Sequential solver (reference implementation)
- Parallel solver using PThreads (pairwise O(N^2) split)
- SDL2 viewer (depth-sorted sprites, perspective camera, gentle depth dim/size)
- Tweakable integrator params: DT, SUBSTEPS, G, EPS2

## Repository Layout
nbody.h          — Body struct, G/EPS2 defines, prototypes
nbody_seq.c      — physics + main (sequential)
viewer.h / viewer.c — SDL2 renderer and camera controls
nbody_parallel.c — parallel version (not implemented)
Makefile         — builds seq (with SDL2) and keeps tree clean (no .o files)

## Prerequisites
- GCC or Clang
- Make
- SDL2 dev headers (viewer only)
  - Ubuntu/Debian: sudo apt-get install libsdl2-dev
  - macOS (Homebrew): brew install sdl2

## Building

Using the provided Makefile (recommended):
make

One-liner without Makefile (Linux/macOS):
gcc -O2 -std=c11 -Wall -Wextra -Wpedantic nbody_seq.c viewer.c -lm $(sdl2-config --cflags --libs) -o nbody_seq

## Running

Sequential:
./nbody_seq <number_of_particles> <timesteps>

Notes:
- timesteps = 0 → run until window close
- Otherwise, the viewer advances physics until exactly <timesteps> substeps have been simulated (bundled as SUBSTEPS per frame), then exits.

Examples:
./nbody_seq 200 0
./nbody_seq 1500 100000

Parallel (not implemented yet):
./nbody_parallel <number_of_particles> <timesteps> <threads>

## Controls (Viewer)
- W / S: forward / backward
- A / D: strafe left / right
- Q / E: move down / up
- Arrow keys: look around (yaw/pitch)
- Space: pause/resume simulation
- R: reset camera to a look-at-origin pose (about −30 degrees)
- Esc: quit

## Physics Model
Bodies are point masses with gravitational force:
F_ij = G * m_i * m_j / (r^2 + EPS2)^(3/2) * r_vec

Integrator: semi-implicit (symplectic) Euler
- v(t+dt) = v(t) + a(t) * dt
- x(t+dt) = x(t) + v(t+dt) * dt

Forces are computed pairwise with Newton’s 3rd law to avoid double-counting.

## Tuning and Defaults
- G (in nbody.h): gravitational constant (default 1.0)
- EPS2 (in nbody.h): softening (try 1e-6 to 1e-2)
- DT (in viewer): per-substep timestep (default 1e-3)
- SUBSTEPS (in viewer): physics substeps per frame (default 8)

Viewer sprite sizing:
- RMIN_PX, RMAX_PX: base pixel radius range (mass to pixels)
- Distance scaling uses a reference Z_REF so near objects can appear larger and far ones smaller, without blowing up in size.

Stability tips:
- If close encounters explode, reduce DT or increase SUBSTEPS, and/or increase EPS2.
- For prettier demos, add small random transverse velocities in init_bodies.

## License
MIT

## Author
Bailey Rossiter
Queensland University of Technology
