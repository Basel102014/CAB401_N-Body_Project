#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include "nbody.h" // Body, compute_forces, update_bodies
#include "viewer.h"

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

// draw a filled circle
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

void viewer_run(Body *bodies, size_t n, size_t steps)
{
    // mass range (for mass→pixels)
    double mmin = (n ? bodies[0].mass : 1.0), mmax = mmin;
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
        return;
    }
    SDL_Window *win = SDL_CreateWindow("N-Body 3D (SDL2)",
                                       SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIN_W, WIN_H, SDL_WINDOW_SHOWN);
    if (!win)
    {
        fprintf(stderr, "SDL_CreateWindow: %s\n", SDL_GetError());
        SDL_Quit();
        return;
    }
    SDL_Renderer *ren = SDL_CreateRenderer(win, -1,
                                           SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!ren)
    {
        fprintf(stderr, "SDL_CreateRenderer: %s\n", SDL_GetError());
        SDL_DestroyWindow(win);
        SDL_Quit();
        return;
    }

    // Camera settings
    double yaw = -45.0 * (M_PI / 180.0);
    double pitch = -30.0 * (M_PI / 180.0);
    double R = 20.0;
    double fx = cos(pitch) * cos(yaw);
    double fy = sin(pitch);
    double fz = -cos(pitch) * sin(yaw);

    Camera cam = {.x = -R * fx, .y = -R * fy, .z = -R * fz, .yaw = yaw, .pitch = pitch, .fov_deg = 60.0};

    const double nearz = 0.1, farz = 1000.0;

    size_t steps_done = 0;
    bool paused = false, running = true;

    Sprite *sprites = malloc(n * sizeof(Sprite));
    if (!sprites)
    {
        perror("malloc sprites");
        SDL_DestroyRenderer(ren);
        SDL_DestroyWindow(win);
        SDL_Quit();
        return;
    }

    while (running)
    {
        // events
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
                    cam.x = -R * fx;
                    cam.y = -R * fy;
                    cam.z = -R * fz;
                    cam.yaw = yaw;
                    cam.pitch = pitch;
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

        // WASD + Q/E fly
        const Uint8 *keys = SDL_GetKeyboardState(NULL);
        double move_speed = 6.0 * DT * SUBSTEPS;
        double fvx = cos(cam.yaw), fvz = -sin(cam.yaw);
        double lvx = sin(cam.yaw), lvz = cos(cam.yaw);
        if (keys[SDL_SCANCODE_W])
        {
            cam.x += fvx * move_speed;
            cam.z += fvz * move_speed;
        }
        if (keys[SDL_SCANCODE_S])
        {
            cam.x -= fvx * move_speed;
            cam.z -= fvz * move_speed;
        }
        if (keys[SDL_SCANCODE_A])
        {
            cam.x += lvx * move_speed;
            cam.z += lvz * move_speed;
        }
        if (keys[SDL_SCANCODE_D])
        {
            cam.x -= lvx * move_speed;
            cam.z -= lvz * move_speed;
        }
        if (keys[SDL_SCANCODE_Q])
            cam.y -= move_speed;
        if (keys[SDL_SCANCODE_E])
            cam.y += move_speed;

        // physics
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

        // render
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

            // depth size (currently: bigger near, smaller far)
            int base_rad = mass_to_radius_px(bodies[i].mass, mmin, mmax);
            double depth_scale = Z_REF / sprites[k].zcam;
            depth_scale = clampd(depth_scale, 0.25, 10.0);
            int rad = (int)lround(base_rad * depth_scale);
            if (rad < 1)
                rad = 1;

            draw_filled_circle(ren, sx, sy, rad);
        }

        // crosshair
        SDL_SetRenderDrawColor(ren, 40, 40, 50, 150);
        SDL_RenderDrawLine(ren, 0, WIN_H / 2, WIN_W, WIN_H / 2);
        SDL_RenderDrawLine(ren, WIN_W / 2, 0, WIN_W / 2, WIN_H);

        SDL_RenderPresent(ren);
    }

    free(sprites);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
}

void viewer_play(const Snapshots *snaps, const double *masses) {
    if (!snaps || !snaps->xyz || snaps->n == 0 || snaps->frames == 0) return;

    const size_t n = snaps->n;

    // mass range for Option A sizing
    double mmin = masses ? masses[0] : 1.0, mmax = mmin;
    if (masses) {
        for (size_t i = 1; i < n; ++i) {
            if (masses[i] < mmin) mmin = masses[i];
            if (masses[i] > mmax) mmax = masses[i];
        }
    }

    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "SDL_Init error: %s\n", SDL_GetError());
        return;
    }
    SDL_Window *win = SDL_CreateWindow("N-Body 3D (Playback)",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIN_W, WIN_H, SDL_WINDOW_SHOWN);
    if (!win) { fprintf(stderr, "SDL_CreateWindow: %s\n", SDL_GetError()); SDL_Quit(); return; }
    SDL_Renderer *ren = SDL_CreateRenderer(win, -1,
        SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!ren) { fprintf(stderr, "SDL_CreateRenderer: %s\n", SDL_GetError());
        SDL_DestroyWindow(win); SDL_Quit(); return; }

    // same default camera (looking at origin, ~-30°)
    double yaw = -45.0 * (M_PI/180.0), pitch = -30.0 * (M_PI/180.0), R = 20.0;
    double fx =  cos(pitch) * cos(yaw);
    double fy =  sin(pitch);
    double fz = -cos(pitch) * sin(yaw);
    Camera cam = { .x=-R*fx, .y=-R*fy, .z=-R*fz, .yaw=yaw, .pitch=pitch, .fov_deg=60.0 };

    const double nearz = 0.1, farz = 1000.0;

    Sprite *sprites = malloc(n * sizeof(Sprite));
    if (!sprites) { perror("malloc sprites"); SDL_DestroyRenderer(ren); SDL_DestroyWindow(win); SDL_Quit(); return; }

    size_t frame = 0;
    bool paused = false, running = true;

    while (running) {
        // events + simple fly controls (same as viewer_run)
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) {
            if (ev.type == SDL_QUIT) running = false;
            if (ev.type == SDL_KEYDOWN) {
                switch (ev.key.keysym.sym) {
                    case SDLK_ESCAPE: running = false; break;
                    case SDLK_SPACE:  paused = !paused; break;
                    case SDLK_r:      cam.x=-R*fx; cam.y=-R*fy; cam.z=-R*fz; cam.yaw=yaw; cam.pitch=pitch; break;
                    case SDLK_UP:     cam.pitch -= 0.04; break;
                    case SDLK_DOWN:   cam.pitch += 0.04; break;
                    case SDLK_LEFT:   cam.yaw   -= 0.05; break;
                    case SDLK_RIGHT:  cam.yaw   += 0.05; break;
                    default: break;
                }
            }
        }
        const Uint8 *keys = SDL_GetKeyboardState(NULL);
        double move_speed = 0.03; // camera-only in playback
        double fvx = cos(cam.yaw), fvz = -sin(cam.yaw);
        double lvx = sin(cam.yaw), lvz =  cos(cam.yaw);
        if (keys[SDL_SCANCODE_W]) { cam.x += fvx*move_speed; cam.z += fvz*move_speed; }
        if (keys[SDL_SCANCODE_S]) { cam.x -= fvx*move_speed; cam.z -= fvz*move_speed; }
        if (keys[SDL_SCANCODE_A]) { cam.x += lvx*move_speed; cam.z += lvz*move_speed; }
        if (keys[SDL_SCANCODE_D]) { cam.x -= lvx*move_speed; cam.z -= lvz*move_speed; }
        if (keys[SDL_SCANCODE_Q]) cam.y -= move_speed;
        if (keys[SDL_SCANCODE_E]) cam.y += move_speed;

        // advance frame ~once per vsync (change step to control playback speed)
        if (!paused) {
            frame++;
            if (frame >= snaps->frames) frame = snaps->frames - 1; // or loop: frame=0;
        }

        SDL_SetRenderDrawColor(ren, 10, 12, 16, 255);
        SDL_RenderClear(ren);

        double fov_rad = cam.fov_deg * M_PI / 180.0;
        double fy_proj = 0.5 * WIN_H / tan(0.5 * fov_rad);
        double fx_proj = fy_proj * ((double)WIN_W / (double)WIN_H);

        const float *src = snaps->xyz + frame * n * 3;
        size_t vis = 0;
        for (size_t i = 0; i < n; ++i) {
            double wx = (double)src[i*3 + 0];
            double wy = (double)src[i*3 + 1];
            double wz = (double)src[i*3 + 2];
            int sx, sy; double zcam;
            if (project_point(&cam, fx_proj, fy_proj, nearz, farz, wx, wy, wz, &sx, &sy, &zcam)) {
                if (sx < -100 || sx > WIN_W+100 || sy < -100 || sy > WIN_H+100) continue;
                sprites[vis].sx = sx; sprites[vis].sy = sy; sprites[vis].zcam = zcam; sprites[vis].i = i; ++vis;
            }
        }

        qsort(sprites, vis, sizeof(Sprite), cmp_sprite_desc);

        for (size_t k = 0; k < vis; ++k) {
            size_t i = sprites[k].i;
            int sx = sprites[k].sx, sy = sprites[k].sy;
            double t = clampd((sprites[k].zcam - 5.0) / 200.0, 0.0, 1.0);
            Uint8 shade = (Uint8)lround(255.0 * (1.0 - 0.25 * t));
            SDL_SetRenderDrawColor(ren, shade, shade, shade, 255);

            int base_rad = mass_to_radius_px(masses ? masses[i] : 1.0, mmin, mmax);
            double depth_scale = Z_REF / sprites[k].zcam;
            depth_scale = clampd(depth_scale, 0.25, 10.0);
            int rad = (int)lround(base_rad * depth_scale);
            if (rad < 1) rad = 1;

            draw_filled_circle(ren, sx, sy, rad);
        }

        // crosshair
        SDL_SetRenderDrawColor(ren, 40, 40, 50, 150);
        SDL_RenderDrawLine(ren, 0, WIN_H/2, WIN_W, WIN_H/2);
        SDL_RenderDrawLine(ren, WIN_W/2, 0, WIN_W/2, WIN_H);

        SDL_RenderPresent(ren);
    }

    free(sprites);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
}