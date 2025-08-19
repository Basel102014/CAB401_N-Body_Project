# Compiler & flags
CC       = gcc
WARN     = -Wall -Wextra -Wpedantic
CSTD     = -std=c11
OPT      = -O2
CFLAGS   = $(WARN) $(CSTD) $(OPT)

# Detect SDL2 (try sdl2-config, then pkg-config)
SDL_CFLAGS := $(shell sdl2-config --cflags 2>/dev/null || pkg-config --cflags sdl2 2>/dev/null)
SDL_LIBS   := $(shell sdl2-config --libs   2>/dev/null || pkg-config --libs   sdl2 2>/dev/null)

# Math & threads
LDLIBS   = -lm
THREADS  = -pthread

# Targets
all: nbody_seq nbody_parallel

# Sequential build (SDL2 viewer)
# If SDL2 isn't found, stop with a helpful error
ifeq ($(strip $(SDL_LIBS)),)
$(error SDL2 development libs not found. Install `libsdl2-dev` (Debian/Ubuntu) or `brew install sdl2` (macOS).)
endif

nbody_seq: nbody_seq.c nbody.h
	$(CC) $(CFLAGS) $(SDL_CFLAGS) -o $@ $< $(SDL_LIBS) $(LDLIBS)

# Parallel (pthread) build
nbody_parallel: nbody_parallel.c nbody.h
	$(CC) $(CFLAGS) -o $@ $< $(THREADS) $(LDLIBS)

clean:
	rm -f nbody_seq nbody_parallel *.o

.PHONY: all clean
