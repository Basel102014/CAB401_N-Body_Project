# ============================================================
#  N-Body Simulation (Sequential Version)
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
LDLIBS = -lm

# --- Targets ---
TARGET_SEQ = nbody_seq
TARGET_PAR = nbody_par

# --- Source files ---
COMMON_SRC   = nbody_tests.c cli_helpers.c test_presets.c viewer.c
SEQ_SRC      = nbody_seq.c
PAR_SRC      = nbody_parallel.c
HEADERS      = nbody.h nbody_tests.h cli_helpers.h test_presets.h viewer.h nbody_parallel.h nbody_seq.h

# --- Default target ---
all: $(TARGET_SEQ)

# ============================================================
#  Build Targets
# ============================================================

# Sequential build
$(TARGET_SEQ): $(SEQ_SRC) $(COMMON_SRC) $(HEADERS)
	@echo "Building sequential version..."
	$(CC) $(CFLAGS) $(SDL_CFLAGS) -o $@ $(SEQ_SRC) $(COMMON_SRC) $(SDL_LIBS) $(LDLIBS)

# Parallel build (optional future extension)
$(TARGET_PAR): $(PAR_SRC) $(COMMON_SRC) $(HEADERS)
	@echo "Building parallel version..."
	$(CC) $(CFLAGS) $(SDL_CFLAGS) -o $@ $(PAR_SRC) $(COMMON_SRC) -lpthread $(SDL_LIBS) $(LDLIBS)

# ============================================================
#  Utility Targets
# ============================================================

# Clean build artifacts
clean:
	@echo "Cleaning build outputs..."
	rm -f $(TARGET_SEQ) $(TARGET_PAR)

# Phony targets
.PHONY: all clean
