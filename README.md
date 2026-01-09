# N‑Body Simulation (3D)

A high‑performance 3D N‑Body gravitational simulation written in C.

This project was developed as part of a university parallel and systems programming unit and represents a **complete implementation**, including:
- a reference **sequential solver**
- a **parallel solver using POSIX threads**
- a real‑time **SDL2 visualiser**
- benchmarking and correctness testing (documented in the accompanying report)

The project focuses on low‑level data layout, numerical integration, and performance scaling on modern CPUs.

---

## Overview

The simulation models a system of point masses interacting via Newtonian gravity in three dimensions. Forces are computed pairwise with softening to avoid singularities, and motion is integrated forward in time using a symplectic method.

The codebase was intentionally written at a low level (plain C, no external physics libraries) to make data flow, performance bottlenecks, and parallelisation strategies explicit.

---

## Features

- 3D Newtonian gravity with softening (EPS2)
- Fully implemented **sequential solver** (reference implementation)
- Fully implemented **parallel solver using PThreads** (pairwise O(N²) workload split)
- Real‑time SDL2 viewer with free‑fly camera
- Depth‑sorted sprite rendering with perspective scaling
- Tweakable integrator parameters: `DT`, `SUBSTEPS`, `G`, `EPS2`
- Deterministic test presets and benchmarking support

---

## Repository Layout

- `nbody.h` — Body struct, physical constants, shared prototypes
- `nbody_seq.c` — Sequential physics solver and simulation loop
- `nbody_parallel.c` — Parallel solver using POSIX threads
- `viewer.h` / `viewer.c` — SDL2 renderer and camera controls
- `main.c` — Program entry point and CLI handling
- `test_suite.c`, `nbody_tests.c`, `test_presets.c` — Correctness tests and presets
- `benchmark.sh` — Performance benchmarking script
- `Makefile` — Build targets for sequential and parallel binaries

---

## Prerequisites

- GCC or Clang
- Make
- SDL2 development headers (viewer only)

### Installing SDL2

**Ubuntu / Debian**
```bash
sudo apt-get install libsdl2-dev
```

**macOS (Homebrew)**
```bash
brew install sdl2
```

---

## Building

Using the provided Makefile (recommended):
```bash
make
```

This produces:
- `nbody_seq` — sequential solver with SDL2 visualiser
- `nbody_parallel` — parallel solver

Manual build (Linux/macOS, sequential only):
```bash
gcc -O2 -std=c11 -Wall -Wextra -Wpedantic \
    nbody_seq.c viewer.c main.c -lm \
    $(sdl2-config --cflags --libs) -o nbody_seq
```

---

## Running

### Sequential

```bash
./nbody_seq <number_of_particles> <timesteps>
```

Notes:
- `timesteps = 0` runs the simulation until the window is closed
- Otherwise, the simulation advances until exactly `<timesteps>` substeps have been simulated, then exits

Examples:
```bash
./nbody_seq 200 0
./nbody_seq 1500 100000
```

### Parallel

```bash
./nbody_parallel <number_of_particles> <timesteps> <threads>
```

The parallel solver partitions the O(N²) force computation across worker threads and synchronises per‑step updates.

---

## Viewer Controls

- **W / S** — forward / backward
- **A / D** — strafe left / right
- **Q / E** — move down / up
- **Arrow keys** — look around (yaw / pitch)
- **Space** — pause / resume simulation
- **R** — reset camera to look‑at‑origin pose
- **Esc** — quit

---

## Physics Model

Bodies are modelled as point masses interacting under Newtonian gravity:

```
F_ij = G · m_i · m_j / (r² + EPS2)^(3/2) · r⃗
```

Softening (`EPS2`) prevents numerical instability during close encounters.

### Integration

A semi‑implicit (symplectic) Euler integrator is used:

```
v(t + Δt) = v(t) + a(t) · Δt
x(t + Δt) = x(t) + v(t + Δt) · Δt
```

Forces are computed pairwise using Newton’s third law to avoid double‑counting.

---

## Tuning and Defaults

- `G` (in `nbody.h`) — gravitational constant (default: `1.0`)
- `EPS2` (in `nbody.h`) — softening factor (typical range: `1e‑6` to `1e‑2`)
- `DT` (in `viewer`) — per‑substep timestep (default: `1e‑3`)
- `SUBSTEPS` (in `viewer`) — physics substeps per frame (default: `8`)

### Stability Tips

- If close encounters cause instability, reduce `DT`, increase `SUBSTEPS`, and/or increase `EPS2`
- For visually pleasing demos, initialise bodies with small transverse velocities

---

## Performance Notes

Benchmarking shows strong scaling with increasing thread counts, including near‑superlinear speedup in some configurations due to improved cache locality on modern CPUs.

Detailed performance analysis and correctness validation are documented in the accompanying project report.

---

## Licence

MIT

---

## Author

**Bailey Rossiter**  
Queensland University of Technology

