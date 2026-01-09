# CAB401 Parallelisation Project – N-Body Simulation (3D)

A high-performance **3D N-Body gravitational simulation** written in **C**, featuring both **sequential** and **parallel (POSIX Threads)** implementations. This project was developed for CAB401 and represents a **complete, end-to-end parallelisation study**, including correctness validation, benchmarking, and performance analysis.

The simulation models Newtonian gravitational interactions between many bodies and demonstrates how low-level data layout and explicit parallel decomposition impact scalability on modern CPUs.

All simulations are **deterministic**: runs with identical parameters produce identical results due to fixed seeding during initialisation.

---

## Overview

The program is built as a **single binary** that selects execution mode at runtime via a thread count:

- **Sequential mode**: `--threads=1` (reference implementation)
- **Parallel mode**: `--threads>1` (POSIX Threads-based solver)

Optional components included in the project:

- **Validation suite** (`test_suite.c`) — verifies numerical equivalence between sequential and parallel solvers
- **SDL2 Viewer** (`viewer.c`) — optional real-time 3D visualisation with free-fly camera
- **Test presets** (`test_presets.c/h`) — predefined starting configurations (solar system, binary/trinary systems, clusters, lattice grids)
- **Benchmarking scripts** — automated performance measurement and CSV output

The sequential solver was implemented first and serves as the **ground-truth reference**, making correctness verification of the parallel solver explicit and robust.

---

## Features

- 3D Newtonian gravity with softening (`EPS2`)
- Deterministic initialisation (fixed random seed)
- Fully implemented **sequential solver**
- Fully implemented **parallel solver using POSIX Threads**
- Pairwise O(N²) force computation with explicit workload partitioning
- Real-time SDL2 visualiser with perspective camera
- Depth-sorted sprite rendering with distance-based scaling and dimming
- Configurable integrator parameters (`DT`, `SUBSTEPS`, `G`, `EPS2`)
- Automated validation and benchmarking support

---

## Repository Layout

- `main.c` — Program entry point and CLI handling
- `nbody.h` — Body structure, physical constants, shared prototypes
- `nbody_seq.c` — Sequential physics solver (reference implementation)
- `nbody_parallel.c` — Parallel solver using POSIX Threads
- `viewer.h` / `viewer.c` — SDL2 renderer and camera controls
- `test_suite.c` — Validation harness for solver equivalence
- `nbody_tests.c` — Core correctness tests
- `test_presets.c/h` — Predefined simulation setups
- `benchmark.sh` — Automated benchmarking script
- `Makefile` — Build targets for sequential, parallel, and test binaries
- `cab401-report.pdf` — Detailed performance analysis and discussion

---

## Hardware Requirements

### Minimum
- Dual-core x86_64 CPU
- 4 GB RAM
- ~100 MB available storage
- Linux-based OS (tested on Ubuntu 22.04)

### Recommended
- 8-core x86_64 CPU
- ≥8 GB RAM
- Modern cache-heavy CPUs benefit parallel performance
- SDL2 development libraries for visualisation
- Valgrind / Callgrind for profiling (optional)

---

## Software Requirements

- **Compiler**: GCC 13.2.0+ or Clang (C11)
- **Libraries**:
  - POSIX Threads (`libpthread`)
  - Math library (`libm`)
  - SDL2 (optional, viewer only)

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

## Compilation

Using the provided Makefile (recommended):
```bash
make
```

This produces a single executable (e.g. `nbody_sim`).

### Build and Run Tests
```bash
make test
```

### Clean
```bash
make clean
```

> **Note:** Binary names may vary depending on the Makefile configuration. The examples below assume `./nbody_sim`.

---

## Running the Simulation

Execution mode is selected at runtime using `--threads`.

### Common Flags

- `--bodies=N`   number of bodies
- `--steps=N`    total simulation steps
- `--stride=N`   snapshot stride (default: 10)
- `--threads=N`  1 = sequential, >1 = parallel (default: 1)
- `--csv`        write per-run timing to `timing_results.csv`
- `--noview`     disable SDL2 viewer

### Presets

- `--solar`
- `--binary`
- `--trinary`
- `--cluster`
- `--lattice`

### Examples

Sequential baseline:
```bash
./nbody_sim --bodies=2000 --steps=5000 --threads=1 --noview
```

Parallel run (8 threads):
```bash
./nbody_sim --bodies=2000 --steps=5000 --threads=8 --noview
```

Solar system preset with CSV timing:
```bash
./nbody_sim --solar --csv
```

---

## Predefined Simulation Presets

Preset defaults are defined in `test_presets.c`. Command-line arguments override preset defaults.

| Preset      | Bodies | Steps     | Stride | Description |
|------------|--------|-----------|--------|-------------|
| `--solar`   | 9      | 1,000,000 | 10     | Sun + 8 planets |
| `--binary`  | 2      | 500,000   | 10     | Two-body orbit |
| `--trinary` | 3      | 1,000,000 | 10     | Three-body rotating system |
| `--cluster` | 50     | 1,000,000 | 10     | Randomised cluster |
| `--lattice` | 3375   | 10,000    | 10     | Perfect-cube lattice |

---

## Viewer Controls (SDL2)

- **W / S** — forward / backward
- **A / D** — strafe left / right
- **Q / E** — move down / up
- **Arrow keys** — look around (yaw / pitch)
- **Space** — pause / resume simulation
- **R** — reset camera to look-at-origin pose
- **Esc** — quit

---

## Physics Model

Bodies are treated as point masses interacting under Newtonian gravity:

```
F_ij = G · m_i · m_j / (r² + EPS2)^(3/2) · r⃗
```

Softening (`EPS2`) prevents numerical instability during close encounters.

### Integration Scheme

A semi-implicit (symplectic) Euler integrator is used:

```
v(t + Δt) = v(t) + a(t) · Δt
x(t + Δt) = x(t) + v(t + Δt) · Δt
```

Forces are computed pairwise using Newton’s third law to avoid double-counting.

---

## Performance and Scaling

Benchmarking demonstrates strong parallel scaling with increasing thread counts, including **near-superlinear speedup** in some configurations. This behaviour is attributed to improved cache locality as the workload is partitioned across cores on modern CPUs.

Detailed performance results, validation methodology, and discussion of scaling limits are provided in `cab401-report.pdf`.

---

## Licence

MIT

---

## Author

**Bailey Rossiter**  
Student ID: 11326158  
Queensland University of Technology

