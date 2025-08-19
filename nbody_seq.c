#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <SDL2/SDL.h>
#include "nbody.h"

// ------ Tweakables ------
#define WIN_W 960
#define WIN_H 600
#define SUBSTEPS 8 // physics substeps per frame
#define DT 1e-3    // physics dt per substep

// mass -> pixel size
#define RMIN_PX 2
#define RMAX_PX 12

// distance-to-size curve: at zcam == Z_REF, scale = 1.0
#define Z_REF 20.0
// ------------------------

// Simple camera
typedef struct
{
    double x, y, z; // position
    double yaw;     // rotate around Y (left/right)
    double pitch;   // rotate around X (up/down)
    double fov_deg; // vertical field of view (degrees)
} Camera;

// Sprite for sorting by depth
typedef struct
{
    double zcam; // depth in camera space (forward)
    int sx, sy;  // screen position
    size_t i;    // index of body
} Sprite;

// --- Utilities ---
static inline double clampd(double v, double a, double b)
{
    return v < a ? a : (v > b ? b : v);
}

// draw a filled circle (no extra libs)
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
    return (int)lround(clampd(r, RMIN_PX, RMAX_PX));
}

// project a world point to screen using the camera (perspective)
// Returns false if behind camera or outside near/far
static bool project_point(const Camera *cam, double fx, double fy,
                          double nearz, double farz,
                          double wx, double wy, double wz,
                          int *out_sx, int *out_sy, double *zcam_out)
{
    // Translate to camera
    double px = wx - cam->x;
    double py = wy - cam->y;
    double pz = wz - cam->z;

    // Yaw (around Y)
    double c_yaw = cos(cam->yaw), s_yaw = sin(cam->yaw);
    double x1 = c_yaw * px + s_yaw * pz;
    double z1 = -s_yaw * px + c_yaw * pz;

    // Pitch (around X)
    double c_pitch = cos(cam->pitch), s_pitch = sin(cam->pitch);
    double y2 = c_pitch * py - s_pitch * z1;
    double z2 = s_pitch * py + c_pitch * z1;

    if (z2 <= nearz || z2 > farz)
        return false; // cull behind/too far

    // Perspective divide
    int cx = WIN_W / 2, cy = WIN_H / 2;
    *out_sx = (int)lround(cx + (x1 * fx) / z2);
    *out_sy = (int)lround(cy - (y2 * fy) / z2);
    if (zcam_out)
        *zcam_out = z2;
    return true;
}

// depth sort: draw far -> near
static int cmp_sprite_desc(const void *a, const void *b)
{
    double za = ((const Sprite *)a)->zcam;
    double zb = ((const Sprite *)b)->zcam;
    return (zb > za) - (zb < za); // descending by z (far first)
}

// ---- Your physics (unchanged) ----
void init_bodies(Body *bodies, size_t n)
{
    for (size_t i = 0; i < n; ++i)
    {
        bodies[i].mass = 1.0;                            // tweak to see different sizes
        bodies[i].x = (double)i - (double)(n - 1) * 0.5; // center the line
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

// ---- Main with 3D camera/projection ----
int main(int argc, char **argv)
{
    if (argc < 3)
    {
        fprintf(stderr, "Usage: %s number_of_particles timesteps\n", argv[0]);
        return EXIT_FAILURE;
    }

    size_t n = strtoull(argv[1], NULL, 10);
    size_t steps = strtoull(argv[2], NULL, 10); // 0 = run until close

    Body *bodies = calloc(n, sizeof(Body));
    if (!bodies)
    {
        perror("calloc");
        return EXIT_FAILURE;
    }
    init_bodies(bodies, n);

    // mass range for Option A sizing
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
    SDL_Window *win = SDL_CreateWindow("N-Body 3D (SDL2)",
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

    double yaw = -45.0 * (M_PI/180.0);
    double pitch = -30.0 * (M_PI/180.0);
    double R = 20.0;  // distance from origin

    double fx =  cos(pitch) * cos(yaw);
    double fy =  sin(pitch);
    double fz = -cos(pitch) * sin(yaw);

    Camera cam = {
        .x = -R*fx,
        .y = -R*fy,
        .z = -R*fz,
        .yaw = yaw,
        .pitch = pitch,
        .fov_deg = 60.0
    };


    // Projection constants
    const double nearz = 0.1;
    const double farz = 1000.0;

    size_t steps_done = 0;
    bool paused = false, running = true;

    // temp buffers for sprites
    Sprite *sprites = malloc(n * sizeof(Sprite));
    if (!sprites)
    {
        perror("malloc");
        return EXIT_FAILURE;
    }

    while (running)
    {
        // Handle input
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
                    cam.x = 0.0;
                    cam.y = 2.0;
                    cam.z = 15.0;
                    cam.yaw = 0.0;
                    cam.pitch = -0.15;
                    break;
                case SDLK_UP:
                    cam.pitch -= 0.04;
                    break;
                case SDLK_DOWN:
                    cam.pitch += 0.04;
                    break;
                case SDLK_LEFT:
                    cam.yaw -= 0.05;
                    break;
                case SDLK_RIGHT:
                    cam.yaw += 0.05;
                    break;
                default:
                    break;
                }
            }
        }

        // Smooth movement with held keys (WASD + Q/E)
        const Uint8 *keys = SDL_GetKeyboardState(NULL);
        double move_speed = 6.0 * DT * SUBSTEPS;      // world units per visual frame
        double fx = cos(cam.yaw), fz = -sin(cam.yaw); // forward (xz plane)
        double lx = sin(cam.yaw), lz = cos(cam.yaw);  // left
        if (keys[SDL_SCANCODE_W])
        {
            cam.x += fx * move_speed;
            cam.z += fz * move_speed;
        }
        if (keys[SDL_SCANCODE_S])
        {
            cam.x -= fx * move_speed;
            cam.z -= fz * move_speed;
        }
        if (keys[SDL_SCANCODE_A])
        {
            cam.x += lx * move_speed;
            cam.z += lz * move_speed;
        }
        if (keys[SDL_SCANCODE_D])
        {
            cam.x -= lx * move_speed;
            cam.z -= lz * move_speed;
        }
        if (keys[SDL_SCANCODE_Q])
        {
            cam.y -= move_speed;
        }
        if (keys[SDL_SCANCODE_E])
        {
            cam.y += move_speed;
        }

        // Physics
        if (!paused)
        {
            int substeps_this_frame = SUBSTEPS;
            if (steps > 0 && steps_done + SUBSTEPS > steps)
            {
                substeps_this_frame = (int)(steps - steps_done);
            }
            for (int s = 0; s < substeps_this_frame; ++s)
            {
                compute_forces(bodies, n);
                update_bodies(bodies, n, DT);
            }
            steps_done += substeps_this_frame;
            if (steps > 0 && steps_done >= steps)
                running = false;
        }

        // Render
        SDL_SetRenderDrawColor(ren, 10, 12, 16, 255);
        SDL_RenderClear(ren);

        double fov_rad = cam.fov_deg * M_PI / 180.0;
        double fy_proj = 0.5 * WIN_H / tan(0.5 * fov_rad);
        double fx_proj = fy_proj * ((double)WIN_W / (double)WIN_H);

        size_t vis = 0;
        for (size_t i = 0; i < n; ++i)
        {
            int sx, sy;
            double zcam;
            if (project_point(&cam, fx_proj, fy_proj, nearz, farz,
                              bodies[i].x, bodies[i].y, bodies[i].z,
                              &sx, &sy, &zcam))
            {
                if (sx < -100 || sx > WIN_W + 100 || sy < -100 || sy > WIN_H + 100)
                    continue;
                sprites[vis].sx = sx;
                sprites[vis].sy = sy;
                sprites[vis].zcam = zcam;
                sprites[vis].i = i;
                ++vis;
            }
        }

        qsort(sprites, vis, sizeof(Sprite), cmp_sprite_desc);

        for (size_t k = 0; k < vis; ++k)
        {
            size_t i = sprites[k].i;
            int sx = sprites[k].sx, sy = sprites[k].sy;

            // gentle dimming by depth (near=255, far≈191)
            double t = clampd((sprites[k].zcam - 5.0) / 200.0, 0.0, 1.0);
            Uint8 shade = (Uint8)lround(255.0 * (1.0 - 0.25 * t));
            SDL_SetRenderDrawColor(ren, shade, shade, shade, 255);

            // depth size: never bigger than base, shrinks with distance
            int base_rad = mass_to_radius_px(bodies[i].mass, mmin, mmax);
            double depth_scale = Z_REF / sprites[k].zcam;
            depth_scale = clampd(depth_scale, 0.25, 10.00);
            int rad = (int)lround(base_rad * depth_scale);
            if (rad < 1)
                rad = 1;

            draw_filled_circle(ren, sx, sy, rad);
        }

        // crosshair at screen center
        SDL_SetRenderDrawColor(ren, 40, 40, 50, 150);
        SDL_RenderDrawLine(ren, 0, WIN_H / 2, WIN_W, WIN_H / 2);
        SDL_RenderDrawLine(ren, WIN_W / 2, 0, WIN_W / 2, WIN_H);

        SDL_RenderPresent(ren);
    }

    free(sprites);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
    free(bodies);
    return EXIT_SUCCESS;
}