#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <SDL2/SDL.h>
#include "nbody.h"

// ------ Tweakables ------
#define WIN_W  960
#define WIN_H  600
#define DOT_SZ 4          // pixel size of each body
#define SUBSTEPS 8        // physics substeps per frame for stability
// ------------------------

void init_bodies(Body *bodies, size_t n)
{
    for (size_t i = 0; i < n; ++i)
    {
        bodies[i].mass = 1.0;    // Placeholder mass
        bodies[i].x = (double)i; // Placeholder position
        bodies[i].y = 0.0;
        bodies[i].z = 0.0;
        bodies[i].vx = 0.0;
        bodies[i].vy = 0.0;
        bodies[i].vz = 0.0;
        bodies[i].fx = 0.0;
        bodies[i].fy = 0.0;
        bodies[i].fz = 0.0;
    }
    printf("Initialized %zu bodies.\n", n);
    fflush(stdout);
}

void compute_forces(Body *bodies, size_t n)
{
    // Pairwise interactions, apply equal & opposite forces
    for (size_t i = 0; i < n; ++i)
    {
        const double mix = bodies[i].mass;
        const double xi = bodies[i].x, yi = bodies[i].y, zi = bodies[i].z;

        for (size_t j = i + 1; j < n; ++j)
        {
            const double dx = bodies[j].x - xi;
            const double dy = bodies[j].y - yi;
            const double dz = bodies[j].z - zi;

            // Softened distance^2 (avoids singularities / huge forces)
            const double r2 = dx * dx + dy * dy + dz * dz + EPS2;
            const double invr = 1.0 / sqrt(r2);
            const double invr3 = invr * invr * invr;

            // Scalar for vector force: G * m_i * m_j / r^3
            const double s = G * mix * bodies[j].mass * invr3;

            const double fx = s * dx;
            const double fy = s * dy;
            const double fz = s * dz;

            bodies[i].fx += fx;
            bodies[i].fy += fy;
            bodies[i].fz += fz;
            bodies[j].fx -= fx;
            bodies[j].fy -= fy;
            bodies[j].fz -= fz;
        }
    }
    // printf("Computed forces for %zu bodies.\n", n);
    // fflush(stdout);
}

void update_bodies(Body *bodies, size_t n, double dt)
{
    for (size_t i = 0; i < n; ++i)
    {
        const double m = bodies[i].mass;
        if (m <= 0.0)
            continue; // guard against invalid mass

        // a = F / m
        const double ax = bodies[i].fx / m;
        const double ay = bodies[i].fy / m;
        const double az = bodies[i].fz / m;

        // semi-implicit Euler: v(t+dt) = v(t) + a*dt
        bodies[i].vx += ax * dt;
        bodies[i].vy += ay * dt;
        bodies[i].vz += az * dt;

        // then x(t+dt) = x(t) + v(t+dt)*dt
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
        bodies[i].z += bodies[i].vz * dt;
    }

    for (size_t i = 0; i < n; ++i)
    {
        bodies[i].fx = bodies[i].fy = bodies[i].fz = 0.0;
    }
    // printf("Updated positions and velocities for %zu bodies.\n", n);
    // fflush(stdout);
}

static inline void world_to_screen(double x, double y, double scale,
                                   double panx, double pany, int *sx, int *sy) {
    *sx = (int)lround(panx + x * scale);
    *sy = (int)lround(pany - y * scale); // world +y up -> screen -y
}


int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "Usage: %s number_of_particles timesteps\n", argv[0]);
        return EXIT_FAILURE;
    }

    size_t n = (size_t)atoi(argv[1]);
    size_t steps = (size_t)atoi(argv[2]);  // only used if you want a fixed run
    (void)steps;                            // viewer runs until window close

    Body *bodies = calloc(n, sizeof(Body));
    if (!bodies) { perror("calloc"); return EXIT_FAILURE; }

    // --- init physics ---
    init_bodies(bodies, n);

    // --- init SDL ---
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "SDL_Init error: %s\n", SDL_GetError());
        return EXIT_FAILURE;
    }

    SDL_Window *win = SDL_CreateWindow("N-Body (SDL2)",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        WIN_W, WIN_H, SDL_WINDOW_SHOWN);
    if (!win) { fprintf(stderr, "SDL_CreateWindow: %s\n", SDL_GetError()); return EXIT_FAILURE; }

    SDL_Renderer *ren = SDL_CreateRenderer(win, -1,
        SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!ren) { fprintf(stderr, "SDL_CreateRenderer: %s\n", SDL_GetError()); return EXIT_FAILURE; }

    // --- view state ---
    double scale = 60.0;                 // pixels per world unit
    double panx  = WIN_W * 0.5;          // center world origin on screen
    double pany  = WIN_H * 0.5;
    bool paused  = false;

    // physics timestep (world seconds)
    // keep small; we substep for stability at 60 FPS
    const double dt = 1e-3;

    // --- main loop ---
    bool running = true;
    while (running) {
        // events
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) {
            if (ev.type == SDL_QUIT) running = false;
            if (ev.type == SDL_KEYDOWN) {
                switch (ev.key.keysym.sym) {
                    case SDLK_ESCAPE: running = false; break;
                    case SDLK_SPACE:  paused = !paused; break;
                    case SDLK_r:      // reset view
                        scale = 60.0; panx = WIN_W*0.5; pany = WIN_H*0.5; break;
                    case SDLK_EQUALS: // '+' key
                    case SDLK_PLUS:   scale *= 1.1; break;
                    case SDLK_MINUS:  scale /= 1.1; break;
                    case SDLK_UP:     pany += 20; break;
                    case SDLK_DOWN:   pany -= 20; break;
                    case SDLK_LEFT:   panx += 20; break;
                    case SDLK_RIGHT:  panx -= 20; break;
                    default: break;
                }
            }
        }

        // physics
        if (!paused) {
            for (int s = 0; s < SUBSTEPS; ++s) {
                compute_forces(bodies, n);
                update_bodies(bodies, n, dt);
            }
        }

        // render
        SDL_SetRenderDrawColor(ren, 10, 12, 16, 255); // background
        SDL_RenderClear(ren);

        // optional: draw crosshair at world origin
        int cx, cy;
        world_to_screen(0.0, 0.0, scale, panx, pany, &cx, &cy);
        SDL_SetRenderDrawColor(ren, 40, 40, 50, 255);
        SDL_RenderDrawLine(ren, 0, cy, WIN_W, cy);
        SDL_RenderDrawLine(ren, cx, 0, cx, WIN_H);

        // draw bodies
        SDL_SetRenderDrawColor(ren, 230, 230, 230, 255);
        for (size_t i = 0; i < n; ++i) {
            int sx, sy;
            world_to_screen(bodies[i].x, bodies[i].y, scale, panx, pany, &sx, &sy);
            SDL_Rect r = { sx - DOT_SZ/2, sy - DOT_SZ/2, DOT_SZ, DOT_SZ };
            SDL_RenderFillRect(ren, &r);
        }

        SDL_RenderPresent(ren);
        // VSYNC caps to ~60 FPS; no manual delay needed
    }

    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
    free(bodies);
    return EXIT_SUCCESS;
}
