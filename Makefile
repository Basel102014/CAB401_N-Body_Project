CC=gcc
CFLAGS=-Wall -Wextra -std=c11 -O2
LDFLAGS=-lm
PTHREAD=-lpthread

all: nbody_seq nbody_parallel

nbody_seq: nbody_seq.c nbody.h
	$(CC) $(CFLAGS) -o $@ nbody_seq.c $(LDFLAGS)

nbody_parallel: nbody_parallel.c nbody.h
	$(CC) $(CFLAGS) -o $@ nbody_parallel.c $(PTHREAD) $(LDFLAGS)

clean:
	rm -f nbody_seq nbody_parallel

.PHONY: all clean
