# Makefile (no leftover .o files)

CC       = gcc
WARN     = -Wall -Wextra -Wpedantic
CSTD     = -std=c11
OPT      = -O2
CFLAGS   = $(WARN) $(CSTD) $(OPT)

# SDL2 flags (uses sdl2-config, falls back to pkg-config)
SDL_CFLAGS := $(shell sdl2-config --cflags 2>/dev/null || pkg-config --cflags sdl2 2>/dev/null)
SDL_LIBS   := $(shell sdl2-config --libs   2>/dev/null || pkg-config --libs   sdl2 2>/dev/null)

LDLIBS   = -lm

TARGETS = nbody_seq

all: $(TARGETS)

# Build in a single command — no intermediate .o files are produced
nbody_seq: nbody_seq.c viewer.c nbody.h viewer.h
	$(CC) $(CFLAGS) $(SDL_CFLAGS) -o $@ nbody_seq.c viewer.c $(SDL_LIBS) $(LDLIBS)

clean:
	rm -f $(TARGETS)

.PHONY: all clean
