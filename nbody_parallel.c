#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "nbody.h"

/* Placeholder implementations shared with sequential version */
void init_bodies(Body *bodies, size_t n) { (void)bodies; (void)n; /* TODO: initialize masses and positions */ }
void compute_forces(Body *bodies, size_t n) { (void)bodies; (void)n; /* TODO: compute forces */ }
void update_bodies(Body *bodies, size_t n, double dt) { (void)bodies; (void)n; (void)dt; /* TODO: update velocities and positions */ }


typedef struct {
    Body *bodies;
    size_t n;
    size_t start;
    size_t end;
} WorkerArgs;

/* Placeholder worker function */
static void *worker(void *arg) {
    WorkerArgs *args = (WorkerArgs *)arg;
    (void)args;
    /* TODO: compute forces for subset and update bodies */
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s number_of_particles timesteps threads\n", argv[0]);
        return EXIT_FAILURE;
    }

    size_t n = (size_t)atoi(argv[1]);
    size_t steps = (size_t)atoi(argv[2]);
    size_t thread_count = (size_t)atoi(argv[3]);

    Body *bodies = calloc(n, sizeof(Body));
    pthread_t *threads = calloc(thread_count, sizeof(pthread_t));
    WorkerArgs *args = calloc(thread_count, sizeof(WorkerArgs));
    if (!bodies || !threads || !args) {
        perror("calloc");
        free(bodies);
        free(threads);
        free(args);
        return EXIT_FAILURE;
    }

    init_bodies(bodies, n);
    for (size_t step = 0; step < steps; ++step) {
        /* TODO: divide work among threads */
        for (size_t t = 0; t < thread_count; ++t) {
            args[t].bodies = bodies;
            args[t].n = n;
            args[t].start = 0; /* TODO */
            args[t].end = n;   /* TODO */
            pthread_create(&threads[t], NULL, worker, &args[t]);
        }
        for (size_t t = 0; t < thread_count; ++t) {
            pthread_join(threads[t], NULL);
        }
    }

    printf("Parallel simulation placeholder complete.\n");
    free(args);
    free(threads);
    free(bodies);
    return EXIT_SUCCESS;
}
