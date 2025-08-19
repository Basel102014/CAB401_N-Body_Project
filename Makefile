# Makefile
CC       = gcc
WARN     = -Wall -Wextra -Wpedantic
CSTD     = -std=c11
OPT      = -O2
CFLAGS   = $(WARN) $(CSTD) $(OPT)

SDL_CFLAGS := $(shell sdl2-config --cflags 2>/dev/null || pkg-config --cflags sdl2 2>/dev/null)
SDL_LIBS   := $(shell sdl2-config --libs   2>/dev/null || pkg-config --libs   sdl2 2>/dev/null)

LDLIBS   = -lm
THREADS  = -pthread

all: nbody_seq

nbody_seq: nbody_seq.o viewer.o
	$(CC) $(CFLAGS) -o $@ $^ $(SDL_LIBS) $(LDLIBS)

nbody_seq.o: nbody_seq.c nbody.h viewer.h
	$(CC) $(CFLAGS) -c $<

viewer.o: viewer.c viewer.h nbody.h
	$(CC) $(CFLAGS) $(SDL_CFLAGS) -c $<

clean:
	rm -f nbody_seq *.o

.PHONY: all clean
