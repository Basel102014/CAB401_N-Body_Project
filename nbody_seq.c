#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <SDL2/SDL.h>
#include "nbody.h"

// ------ Tweakables ------
#define WIN_W 960
#define WIN_H 600
#define SUBSTEPS 8 // physics substeps per frame for stability

// mass -> pixel size
#define RMIN_PX 2  // smallest dot (pixels)
#define RMAX_PX 12 // largest dot (pixels)
// ------------------------

// draw a filled circle (no SDL_gfx needed)
static void draw_filled_circle(SDL_Renderer *ren, int cx, int cy, int r)
{
    if (r <= 0)
        return;
    for (int dy = -r; dy <= r; ++dy)
    {
        int dx = (int)floor(sqrt((double)r * r - (double)dy * dy));
        SDL_RenderDrawLine(ren, cx - dx, cy + dy, cx + dx, cy + dy);
    }
}

// map mass -> pixel radius using cube-root (constant-density look)
static int mass_to_radius_px(double m, double mmin, double mmax)
{
    double a = cbrt(fmax(0.0, m));
    double amin = cbrt(fmax(0.0, mmin));
    double amax = cbrt(fmax(0.0, mmax));
    double t = (amax > amin) ? (a - amin) / (amax - amin) : 0.0;
    double r = RMIN_PX + t * (RMAX_PX - RMIN_PX);
    if (r < RMIN_PX)
        r = RMIN_PX;
    if (r > RMAX_PX)
        r = RMAX_PX;
    return (int)lround(r);
}

void init_bodies(Body *bodies, size_t n)
{
    for (size_t i = 0; i < n; ++i)
    {
        bodies[i].mass = 1.0; // tweak these later to see size differences
        bodies[i].x = (double)i;
        bodies[i].y = 0.0;
        bodies[i].z = 0.0;
        bodies[i].vx = bodies[i].vy = bodies[i].vz = 0.0;
        bodies[i].fx = bodies[i].fy = bodies[i].fz = 0.0;
    }
    printf("Initialized %zu bodies.\n", n);
    fflush(stdout);
}

void compute_forces(Body *bodies, size_t n)
{
    for (size_t i = 0; i < n; ++i)
    {
        const double mi = bodies[i].mass;
        const double xi = bodies[i].x, yi = bodies[i].y, zi = bodies[i].z;
        for (size_t j = i + 1; j < n; ++j)
        {
            const double dx = bodies[j].x - xi;
            const double dy = bodies[j].y - yi;
            const double dz = bodies[j].z - zi;
            const double r2 = dx * dx + dy * dy + dz * dz + EPS2;
            const double invr = 1.0 / sqrt(r2);
            const double invr3 = invr * invr * invr;
            const double s = G * mi * bodies[j].mass * invr3;
            const double fx = s * dx, fy = s * dy, fz = s * dz;
            bodies[i].fx += fx;
            bodies[i].fy += fy;
            bodies[i].fz += fz;
            bodies[j].fx -= fx;
            bodies[j].fy -= fy;
            bodies[j].fz -= fz;
        }
    }
}

void update_bodies(Body *bodies, size_t n, double dt)
{
    for (size_t i = 0; i < n; ++i)
    {
        const double m = bodies[i].mass;
        if (m <= 0.0)
            continue;
        const double ax = bodies[i].fx / m;
        const double ay = bodies[i].fy / m;
        const double az = bodies[i].fz / m;
        bodies[i].vx += ax * dt;
        bodies[i].vy += ay * dt;
        bodies[i].vz += az * dt;
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
        bodies[i].z += bodies[i].vz * dt;
        bodies[i].fx = bodies[i].fy = bodies[i].fz = 0.0; // clear for next step
    }
}

static inline void world_to_screen(double x, double y, double scale,
                                   double panx, double pany, int *sx, int *sy)
{
    *sx = (int)lround(panx + x * scale);
    *sy = (int)lround(pany - y * scale); // world +y up -> screen -y
}

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        fprintf(stderr, "Usage: %s number_of_particles timesteps\n", argv[0]);
        return EXIT_FAILURE;
    }

    size_t n = (size_t)atoi(argv[1]);
    size_t steps  = strtoull(argv[2], NULL, 10);

    Body *bodies = calloc(n, sizeof(Body));
    if (!bodies)
    {
        perror("calloc");
        return EXIT_FAILURE;
    }

    init_bodies(bodies, n);

    // compute mass range once (recompute later if your masses change)
    double mmin = bodies[0].mass, mmax = bodies[0].mass;
    for (size_t i = 1; i < n; ++i)
    {
        if (bodies[i].mass < mmin)
            mmin = bodies[i].mass;
        if (bodies[i].mass > mmax)
            mmax = bodies[i].mass;
    }

    if (SDL_Init(SDL_INIT_VIDEO) != 0)
    {
        fprintf(stderr, "SDL_Init error: %s\n", SDL_GetError());
        return EXIT_FAILURE;
    }
    SDL_Window *win = SDL_CreateWindow("N-Body (SDL2)",
                                       SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIN_W, WIN_H, SDL_WINDOW_SHOWN);
    if (!win)
    {
        fprintf(stderr, "SDL_CreateWindow: %s\n", SDL_GetError());
        return EXIT_FAILURE;
    }
    SDL_Renderer *ren = SDL_CreateRenderer(win, -1,
                                           SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!ren)
    {
        fprintf(stderr, "SDL_CreateRenderer: %s\n", SDL_GetError());
        return EXIT_FAILURE;
    }

    double scale = 60.0, panx = WIN_W * 0.5, pany = WIN_H * 0.5;
    bool paused = false;
    const double dt = 1e-3;

    size_t steps_done = 0;
    bool running = true;
    while (running)
    {
        SDL_Event ev;
        while (SDL_PollEvent(&ev))
        {
            if (ev.type == SDL_QUIT)
                running = false;
            if (ev.type == SDL_KEYDOWN)
            {
                switch (ev.key.keysym.sym)
                {
                case SDLK_ESCAPE:
                    running = false;
                    break;
                case SDLK_SPACE:
                    paused = !paused;
                    break;
                case SDLK_r:
                    scale = 60.0;
                    panx = WIN_W * 0.5;
                    pany = WIN_H * 0.5;
                    break;
                case SDLK_EQUALS:
                case SDLK_PLUS:
                    scale *= 1.1;
                    break;
                case SDLK_MINUS:
                    scale /= 1.1;
                    break;
                case SDLK_UP:
                    pany += 20;
                    break;
                case SDLK_DOWN:
                    pany -= 20;
                    break;
                case SDLK_LEFT:
                    panx += 20;
                    break;
                case SDLK_RIGHT:
                    panx -= 20;
                    break;
                }
            }
        }

        if (!paused) {
            // do at most SUBSTEPS each frame, but don’t exceed the requested total
            int substeps_this_frame = SUBSTEPS;
            if (steps > 0 && steps_done + SUBSTEPS > steps) {
                substeps_this_frame = (int)(steps - steps_done);
            }

            for (int s = 0; s < substeps_this_frame; ++s) {
                compute_forces(bodies, n);
                update_bodies(bodies, n, dt);
            }
            steps_done += substeps_this_frame;

            // auto-quit when we’ve hit the target
            if (steps > 0 && steps_done >= steps) {
                running = false;
            }
        }

        SDL_SetRenderDrawColor(ren, 10, 12, 16, 255);
        SDL_RenderClear(ren);

        // crosshair at world origin
        int cx, cy;
        world_to_screen(0.0, 0.0, scale, panx, pany, &cx, &cy);
        SDL_SetRenderDrawColor(ren, 40, 40, 50, 255);
        SDL_RenderDrawLine(ren, 0, cy, WIN_W, cy);
        SDL_RenderDrawLine(ren, cx, 0, cx, WIN_H);

        // draw bodies with mass-scaled pixel radius
        SDL_SetRenderDrawColor(ren, 230, 230, 230, 255);
        for (size_t i = 0; i < n; ++i)
        {
            int sx, sy;
            world_to_screen(bodies[i].x, bodies[i].y, scale, panx, pany, &sx, &sy);
            int rad = mass_to_radius_px(bodies[i].mass, mmin, mmax);
            draw_filled_circle(ren, sx, sy, rad);
        }

        SDL_RenderPresent(ren);
    }

    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
    free(bodies);
    return EXIT_SUCCESS;
}
