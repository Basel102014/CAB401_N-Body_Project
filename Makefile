# ============================================================
#  Fast build for numeric workloads (multi-threaded, math-heavy)
# ============================================================

CC       = gcc
CSTD     = -std=c17
WARN     = -Wall -Wextra -Wpedantic -Wno-unused-parameter -Wno-unused-variable

# --- Optimisation flags ---
OPT      = -O3 -march=native -ffast-math -funroll-loops -fno-math-errno \
           -fno-trapping-math -fomit-frame-pointer -ftree-vectorize \
           -falign-functions=32 -falign-loops=32 -fstrict-aliasing

# --- Thread & SIMD safety ---
CFLAGS   = $(CSTD) $(WARN) $(OPT) -pthread

# --- SDL2 (viewer) ---
SDL_CFLAGS := $(shell sdl2-config --cflags 2>/dev/null || pkg-config --cflags sdl2 2>/dev/null)
SDL_LIBS   := $(shell sdl2-config --libs   2>/dev/null || pkg-config --libs   sdl2 2>/dev/null)

LDLIBS = -lm -lpthread

TARGET = nbody_sim
SRC    = main.c nbody_seq.c nbody_parallel.c nbody_tests.c cli_helpers.c test_presets.c viewer.c
HEADERS = nbody.h nbody_seq.h nbody_parallel.h nbody_tests.h cli_helpers.h test_presets.h viewer.h

all: $(TARGET)

$(TARGET): $(SRC) $(HEADERS)
	@echo "Building N-Body Simulation (Optimised, Parallel)..."
	$(CC) $(CFLAGS) $(SDL_CFLAGS) -o $@ $(SRC) $(SDL_LIBS) $(LDLIBS)

clean:
	@echo "Cleaning build outputs..."
	rm -f $(TARGET)

.PHONY: all clean
