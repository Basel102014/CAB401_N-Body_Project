# N-Body Simulation

## Overview
This project implements a **3D N-Body gravitational simulation** in C.  
It starts with a **sequential version** and includes a **parallel version** using **PThreads** to improve performance on multi-core CPUs.  
The simulation models gravitational interactions among particles in 3D space.

## Features
- Sequential and parallel implementations  
- Parallelism using **PThreads** and shared memory  
- Multiple test cases (e.g., galaxy cluster, uniform distributions)  
- Benchmarking and performance measurement tools  
- Optional visualization support

## Getting Started

### Prerequisites
- **GCC** or compatible C compiler  
- **PThreads** library (usually included)  
- **Make** (optional)

### Building
Compile the sequential version:  
```bash
gcc -o nbody_seq nbody_seq.c -lm
```

Compile the parallel version:  
```bash
gcc -o nbody_parallel nbody_parallel.c -lpthread -lm
```

### Running
Basic usage:  
```bash
./nbody_seq [number_of_particles] [timesteps]
./nbody_parallel [number_of_particles] [timesteps] [threads]
```

Example:  
```bash
./nbody_parallel 1000 100 8
```

## Testing
Includes test cases such as:
- Single massive body with orbiting smaller bodies  
- Equal mass random distributions  
- Clustered bodies  

## Performance
Benchmarking tools and instructions are included to compare sequential and parallel versions.  
Expect **significant speedup** on multi-core CPUs.

## License
This project is licensed under the **MIT License**.

## Author
Bailey Rossiter
Queensland University of Technology
