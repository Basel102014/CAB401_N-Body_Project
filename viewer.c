#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include "nbody.h"
#include "viewer.h"

// --- Window and simulation settings ---
#define WIN_W 960
#define WIN_H 600
#define SUBSTEPS 8
#define DT 1e-3

// --- Rendering constants ---
#define RMIN_PX 2
#define RMAX_PX 12
#define Z_REF 20.0

// --- Camera struct (position + orientation) ---
typedef struct
{
    double x, y, z; // Position in 3D space
    double yaw;     // Horizontal rotation
    double pitch;   // Vertical rotation
    double fov_deg; // Field of view in degrees
} Camera;

// --- Sprite struct for depth sorting ---
typedef struct
{
    double zcam; // Depth from camera
    int sx, sy;  // Screen position
    size_t i;    // Body index
} Sprite;

// Clamp a value between two bounds
static inline double clampd(double v, double a, double b)
{
    return v < a ? a : (v > b ? b : v);
}

// Draw a filled circle (used for bodies)
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

// Convert mass to on-screen radius (keeps visual scaling consistent)
static int mass_to_radius_px(double m, double mmin, double mmax)
{
    double a = cbrt(fmax(0.0, m));
    double amin = cbrt(fmax(0.0, mmin));
    double amax = cbrt(fmax(0.0, mmax));
    double t = (amax > amin) ? (a - amin) / (amax - amin) : 0.0;
    double r = RMIN_PX + t * (RMAX_PX - RMIN_PX);
    return (int)lround(clampd(r, RMIN_PX, RMAX_PX));
}

// Project 3D world coordinates to 2D screen space
static bool project_point(
    const Camera *cam, double fx, double fy,
    double nearz, double farz,
    double wx, double wy, double wz,
    int *out_sx, int *out_sy, double *zcam_out)
{
    double px = wx - cam->x;
    double py = wy - cam->y;
    double pz = wz - cam->z;

    double c_yaw = cos(cam->yaw), s_yaw = sin(cam->yaw);
    double x1 = c_yaw * px + s_yaw * pz;
    double z1 = -s_yaw * px + c_yaw * pz;

    double c_pitch = cos(cam->pitch), s_pitch = sin(cam->pitch);
    double y2 = c_pitch * py - s_pitch * z1;
    double z2 = s_pitch * py + c_pitch * z1;

    if (z2 <= nearz || z2 > farz)
        return false;

    int cx = WIN_W / 2, cy = WIN_H / 2;
    *out_sx = (int)lround(cx + (x1 * fx) / z2);
    *out_sy = (int)lround(cy - (y2 * fy) / z2);
    if (zcam_out)
        *zcam_out = z2;
    return true;
}

// Compare sprites by depth (for rendering back-to-front)
static int cmp_sprite_desc(const void *a, const void *b)
{
    double za = ((const Sprite *)a)->zcam;
    double zb = ((const Sprite *)b)->zcam;
    return (zb > za) - (zb < za);
}

// --- Main viewer function ---
void viewer_play(const Snapshots *snaps, const double *masses)
{
    if (!snaps || !snaps->xyz || snaps->n == 0 || snaps->frames == 0)
        return;

    const size_t n = snaps->n;

    // Determine min/max mass for visual scaling
    double mmin = masses ? masses[0] : 1.0, mmax = mmin;
    if (masses)
    {
        for (size_t i = 1; i < n; ++i)
        {
            if (masses[i] < mmin)
                mmin = masses[i];
            if (masses[i] > mmax)
                mmax = masses[i];
        }
    }

    if (SDL_Init(SDL_INIT_VIDEO) != 0)
    {
        fprintf(stderr, "SDL_Init error: %s\n", SDL_GetError());
        return;
    }

    SDL_Window *win = SDL_CreateWindow(
        "N-Body 3D (Playback)",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        WIN_W, WIN_H, SDL_WINDOW_SHOWN);

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

    SDL_SetRenderDrawBlendMode(ren, SDL_BLENDMODE_BLEND);

    // Default camera setup
    double yaw = -45.0 * (M_PI / 180.0), pitch = -30.0 * (M_PI / 180.0), R = 20.0;
    double fx0 = cos(pitch) * cos(yaw);
    double fy0 = sin(pitch);
    double fz0 = -cos(pitch) * sin(yaw);
    Camera cam = {.x = -R * fx0, .y = -R * fy0, .z = -R * fz0, .yaw = yaw, .pitch = pitch, .fov_deg = 60.0};

    const double nearz = 0.001, farz = 1000.0;
    Sprite *sprites = malloc(n * sizeof(Sprite));
    if (!sprites)
    {
        perror("malloc sprites");
        SDL_DestroyRenderer(ren);
        SDL_DestroyWindow(win);
        SDL_Quit();
        return;
    }

    size_t frame = 0;
    bool paused = false, running = true;

    const int TRAIL_FRAMES_MAX = 200;
    const int TRAIL_STEP = 1;

    // --- Main loop ---
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
                    cam.x = cam.y = cam.z = 0;
                    cam.yaw = cam.pitch = 0;
                    break;
                case SDLK_RIGHT:
                    cam.yaw += 0.05;
                    break;
                case SDLK_LEFT:
                    cam.yaw -= 0.05;
                    break;
                case SDLK_UP:
                    cam.pitch -= 0.04;
                    break;
                case SDLK_DOWN:
                    cam.pitch += 0.04;
                    break;
                }
                const double half_pi = M_PI * 0.5 - 1e-4;
                if (cam.pitch > half_pi)
                    cam.pitch = half_pi;
                if (cam.pitch < -half_pi)
                    cam.pitch = -half_pi;
                if (cam.yaw > M_PI)
                    cam.yaw -= 2.0 * M_PI;
                if (cam.yaw < -M_PI)
                    cam.yaw += 2.0 * M_PI;
            }
        }

        // Camera movement (WASD + Q/E)
        const Uint8 *keys = SDL_GetKeyboardState(NULL);
        double move_speed = 0.3;

        double cosPitch = cos(cam.pitch), sinPitch = sin(cam.pitch);
        double cosYaw = cos(cam.yaw), sinYaw = sin(cam.yaw);

        double fx = sinYaw * cosPitch;
        double fy = sinPitch;
        double fz = -cosYaw * cosPitch;
        double rx = cosYaw;
        double rz = sinYaw;

        if (keys[SDL_SCANCODE_S])
        {
            cam.x += fx * move_speed;
            cam.y -= fy * move_speed;
            cam.z += fz * move_speed;
        }
        if (keys[SDL_SCANCODE_W])
        {
            cam.x -= fx * move_speed;
            cam.y += fy * move_speed;
            cam.z -= fz * move_speed;
        }
        if (keys[SDL_SCANCODE_A])
        {
            cam.x -= rx * move_speed;
            cam.z -= rz * move_speed;
        }
        if (keys[SDL_SCANCODE_D])
        {
            cam.x += rx * move_speed;
            cam.z += rz * move_speed;
        }
        if (keys[SDL_SCANCODE_E])
            cam.y += move_speed;
        if (keys[SDL_SCANCODE_Q])
            cam.y -= move_speed;

        if (!paused)
        {
            frame++;
            if (frame >= snaps->frames)
                frame = snaps->frames - 1;
        }

                SDL_SetRenderDrawColor(ren, 10, 12, 16, 255);
        SDL_RenderClear(ren);

        double fov_rad = cam.fov_deg * M_PI / 180.0;
        double fy_proj = 0.5 * WIN_H / tan(0.5 * fov_rad);
        double fx_proj = fy_proj * ((double)WIN_W / (double)WIN_H);

        // ======================================================
        // 🔹 1. Perspective Grid on XZ Plane (y = 0)
        // ======================================================
        const double grid_extent = 100.0; // world units
        const double grid_step   = 2.0;

        for (double x = -grid_extent; x <= grid_extent; x += grid_step) {
            int sx1, sy1, sx2, sy2;
            double zcam;
            if (project_point(&cam, fx_proj, fy_proj, nearz, farz, x, 0.0, -grid_extent, &sx1, &sy1, &zcam) &&
                project_point(&cam, fx_proj, fy_proj, nearz, farz, x, 0.0,  grid_extent, &sx2, &sy2, &zcam)) {

                double fade = 1.0 - fabs(x) / grid_extent;
                Uint8 alpha = (Uint8)(40 + 60 * fade);
                if (fabs(x) < 1e-6)
                    SDL_SetRenderDrawColor(ren, 160, 100, 100, 180); // highlight Z-axis
                else
                    SDL_SetRenderDrawColor(ren, 80, 80, 100, alpha);
                SDL_RenderDrawLine(ren, sx1, sy1, sx2, sy2);
            }
        }

        for (double z = -grid_extent; z <= grid_extent; z += grid_step) {
            int sx1, sy1, sx2, sy2;
            double zcam;
            if (project_point(&cam, fx_proj, fy_proj, nearz, farz, -grid_extent, 0.0, z, &sx1, &sy1, &zcam) &&
                project_point(&cam, fx_proj, fy_proj, nearz, farz,  grid_extent, 0.0, z, &sx2, &sy2, &zcam)) {

                double fade = 1.0 - fabs(z) / grid_extent;
                Uint8 alpha = (Uint8)(40 + 60 * fade);
                if (fabs(z) < 1e-6)
                    SDL_SetRenderDrawColor(ren, 100, 160, 100, 180); // highlight X-axis
                else
                    SDL_SetRenderDrawColor(ren, 80, 80, 100, alpha);
                SDL_RenderDrawLine(ren, sx1, sy1, sx2, sy2);
            }
        }

        // ======================================================
        // 🔹 2. Trails (fade with time)
        // ======================================================
        int trail_frames = TRAIL_FRAMES_MAX;
        if (trail_frames > (int)frame) trail_frames = (int)frame;

        if (trail_frames > 1) {
            for (size_t i = 0; i < n; ++i) {
                int prev_sx = 0, prev_sy = 0;
                bool has_prev = false;

                for (int t = trail_frames; t > 0; t -= TRAIL_STEP) {
                    size_t fidx = frame - t;
                    const float *src = snaps->xyz + fidx * n * 3;
                    double wx = src[i * 3 + 0];
                    double wy = src[i * 3 + 1];
                    double wz = src[i * 3 + 2];

                    int sx, sy;
                    double zcam;
                    if (project_point(&cam, fx_proj, fy_proj, nearz, farz, wx, wy, wz, &sx, &sy, &zcam)) {
                        if (sx >= -100 && sx <= WIN_W + 100 && sy >= -100 && sy <= WIN_H + 100) {
                            Uint8 alpha = (Uint8)(30 + (trail_frames - t) * (150.0 / fmax(1, trail_frames)));
                            SDL_SetRenderDrawColor(ren, 200, 200, 220, alpha);
                            if (has_prev)
                                SDL_RenderDrawLine(ren, prev_sx, prev_sy, sx, sy);
                            prev_sx = sx;
                            prev_sy = sy;
                            has_prev = true;
                        } else has_prev = false;
                    } else has_prev = false;
                }
            }
        }

        // ======================================================
        // 🔹 3. Bodies (depth shading)
        // ======================================================
        const float *src_now = snaps->xyz + frame * n * 3;
        size_t vis = 0;
        for (size_t i = 0; i < n; ++i) {
            double wx = src_now[i * 3 + 0];
            double wy = src_now[i * 3 + 1];
            double wz = src_now[i * 3 + 2];
            int sx, sy;
            double zcam;
            if (project_point(&cam, fx_proj, fy_proj, nearz, farz, wx, wy, wz, &sx, &sy, &zcam)) {
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

        for (size_t k = 0; k < vis; ++k) {
            size_t i = sprites[k].i;
            int sx = sprites[k].sx, sy = sprites[k].sy;

            // Depth-based brightness (brighter near camera)
            double depth_norm = clampd((sprites[k].zcam - 5.0) / 300.0, 0.0, 1.0);
            double brightness = pow(1.0 - depth_norm, 2.0);
            Uint8 shade = (Uint8)lround(50.0 + 205.0 * brightness);

            SDL_SetRenderDrawColor(ren, shade, shade, shade, 255);

            int base_rad = mass_to_radius_px(masses ? masses[i] : 1.0, mmin, mmax);
            double depth_scale = Z_REF / sprites[k].zcam;
            depth_scale = clampd(depth_scale, 0.25, 10.0);
            int rad = (int)lround(base_rad * depth_scale);
            if (rad < 1) rad = 1;

            draw_filled_circle(ren, sx, sy, rad);
        }

        // Crosshair for orientation
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