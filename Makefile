# ============================================================
#  N-Body Simulation (Sequential + Parallel Unified Build)
# ============================================================

# --- Compiler and flags ---
CC       = gcc
CSTD     = -std=c11
WARN     = -Wall -Wextra -Wpedantic
OPT      = -O2
CFLAGS   = $(CSTD) $(WARN) $(OPT)

# --- SDL2 detection (fallback to pkg-config if sdl2-config missing) ---
SDL_CFLAGS := $(shell sdl2-config --cflags 2>/dev/null || pkg-config --cflags sdl2 2>/dev/null)
SDL_LIBS   := $(shell sdl2-config --libs   2>/dev/null || pkg-config --libs   sdl2 2>/dev/null)

# --- Libraries ---
LDLIBS = -lm -lpthread

# --- Target ---
TARGET = nbody_sim

# --- Source files ---
SRC = main.c nbody_seq.c nbody_parallel.c nbody_tests.c cli_helpers.c test_presets.c viewer.c
HEADERS = nbody.h nbody_seq.h nbody_parallel.h nbody_tests.h cli_helpers.h test_presets.h viewer.h

# --- Default target ---
all: $(TARGET)

# ============================================================
#  Build Rules
# ============================================================

$(TARGET): $(SRC) $(HEADERS)
	@echo "Building unified N-Body simulation (sequential + parallel)..."
	$(CC) $(CFLAGS) $(SDL_CFLAGS) -o $@ $(SRC) $(SDL_LIBS) $(LDLIBS)

# ============================================================
#  Utility Targets
# ============================================================

clean:
	@echo "Cleaning build outputs..."
	rm -f $(TARGET)

.PHONY: all clean
