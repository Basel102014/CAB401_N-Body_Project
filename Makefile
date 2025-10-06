# ===== Makefile (no leftover .o files) =====

CC       = gcc
WARN     = -Wall -Wextra -Wpedantic
CSTD     = -std=c11
OPT      = -O2
CFLAGS   = $(WARN) $(CSTD) $(OPT)

# SDL2 flags (uses sdl2-config, falls back to pkg-config)
SDL_CFLAGS := $(shell sdl2-config --cflags 2>/dev/null || pkg-config --cflags sdl2 2>/dev/null)
SDL_LIBS   := $(shell sdl2-config --libs   2>/dev/null || pkg-config --libs   sdl2 2>/dev/null)

LDLIBS = -lm
TARGET = nbody_seq

# ---- Source files ----
SRC = nbody_seq.c nbody_tests.c cli_helpers.c test_presets.c viewer.c
HEADERS = nbody.h nbody_tests.h cli_helpers.h test_presets.h viewer.h

# ---- Default target ----
all: $(TARGET)

# ---- Build main executable ----
$(TARGET): $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) $(SDL_CFLAGS) -o $@ $(SRC) $(SDL_LIBS) $(LDLIBS)

# ---- Clean up build outputs ----
clean:
	rm -f $(TARGET)

.PHONY: all clean
