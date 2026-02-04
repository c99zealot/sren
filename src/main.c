// @TODO start using a Scene object to organise everything instead of just rendering it separately
// @TODO normal mapping
// @TODO serious OBJ/MTL parser
// @TODO make matrices structs so we can use =
// @TODO view frustum, proper perspective transform
// @TODO logging
// @TODO dev console
// @TODO Renderer object
// @TODO move arena to sren.c
// @TODO multiple viewports (for splitscreen etc.)
// @TODO coloured shadows
// @TODO PNG textures
// @TODO screenshot bind
// @TODO perf
// @TODO cleaner fontset handling (probably just have init_renderer load it)
// @TODO specific functions for rendering models with transparent textures,
//       skip alpha blending when rendering models with opaque textures for performance

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <float.h>

#include <SDL2/SDL.h>

#include "sren.h"
#include "arena.h"

#define TIME(x) {                                                   \
        clock_t start = clock();                                    \
        x;                                                          \
        printf("%s : %.4fms\n", #x,                                 \
                (double)(clock() - start) / CLOCKS_PER_SEC * 1000); \
}

const int g_window_width = 800;
const int g_window_height = 800;

//
// render_axes - renders portions of the xyz axes to the framebuffer
//
void render_axes(Camera *cam) {
        Vec3 o = VEC3(0.0, 0.0, 0.0);
        Vec3 x = VEC3(0.1, 0.0, 0.0);
        Vec3 y = VEC3(0.0, 0.1, 0.0);
        Vec3 z = VEC3(0.0, 0.0, 0.1);

        o = m4v3_mul(cam->view_mat, o);
        x = m4v3_mul(cam->view_mat, x);
        y = m4v3_mul(cam->view_mat, y);
        z = m4v3_mul(cam->view_mat, z);

        if (o.z > 0 || x.z > 0 || y.z > 0 || z.z > 0) {
                return;
        }

        o = m4v3_mul(g_viewport, persp(o));
        x = m4v3_mul(g_viewport, persp(x));
        y = m4v3_mul(g_viewport, persp(y));
        z = m4v3_mul(g_viewport, persp(z));

        if (out_of_view(o) || out_of_view(x) || out_of_view(y) || out_of_view(z)) {
                return;
        }

        line(o, x, VEC4(1, 0, 0, 1));
        line(o, y, VEC4(0, 1, 0, 1));
        line(o, z, VEC4(0, 0, 1, 1));
}

int main(int argc, char **argv) {
        if (argc < 3) {
                fprintf(stderr, "usage: %s <obj> <texture>\n", argv[0]);
                return -1;
        }

        if (SDL_Init(SDL_INIT_VIDEO) != 0) {
                fprintf(stderr, "SDL_Init failure: %s\n", SDL_GetError());
                return -1;
        }

        SDL_Window *window = SDL_CreateWindow("sren", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, g_window_width, g_window_height, SDL_WINDOW_SHOWN);
        if (window == NULL) {
                fprintf(stderr, "SDL_CreateWindow failure: %s\n", SDL_GetError());
                SDL_Quit();
                return -1;
        }

        SDL_Surface *screen = SDL_GetWindowSurface(window);
        if (screen == NULL) {
                fprintf(stderr, "SDL_GetWindowSurface failure: %s\n", SDL_GetError());
                SDL_Quit();
                return -1;
        }
        
        Arena render_arena;
        arena_init(&render_arena, ARENA_RESERVE_DEFAULT, 1 << 24, 1);

        init_renderer(&render_arena, screen->pixels, g_window_width, g_window_height);

        Material plastic = {0.18, 0.65, 0.55, 55};

        Model *main_model = load_model(&render_arena, argv[1], argv[2], NULL, &plastic, 1024, 1024, 0, 0);
        Model *floor_model = load_model(&render_arena, "assets/floor.obj", "assets/floor.tex", NULL, &plastic, 1024, 1024, 0, 0);
        Model *ceiling_model = load_model(&render_arena, "assets/ceiling.obj", "assets/chainlink.tex", NULL, &plastic, 1024, 1024, 0, 0);
        //Scene *scene = mkscene(&render_arena);
        
        Texture *fontset = load_texture(&render_arena, "assets/fontset.tex", 7392, 128);
 
        Light light = {
                .colour = VEC4(1.0, 1.0, 1.0, 1.0),
                .ambient = 0.2,
                .diffuse = 1,
                .specular = 0.7,
                .dropoff = 0.032,
                .subject = VEC3(0, -1, 0),
                .up = VEC3(0, 1, 0),
        };
        init_light(&render_arena, &light, 1024, 1024);

        Camera cam = {
                .pos = VEC3(0.1, 0.4, 1),
                .subject = VEC3(0, 0, 0),
                .up = VEC3(0, 1, 0)
        };
        init_cam(&cam);

        SDL_WarpMouseInWindow(window, g_window_width/2, g_window_height/2);
        SDL_SetRelativeMouseMode(SDL_TRUE);

        const double mouse_sens = 0.02;

        double curs_dx = 0;
        double curs_dy = 0;

        Vec3 cam_vel = VEC3(0, 0, 0);

        int frames_drawn = 0;
        double fps = 0;
        double fps_high = 0;
        double fps_low = DBL_MAX;
        double frame_time = 0;
        double frame_time_high = 0;
        double frame_time_low = DBL_MAX;

        double i = 0.0;
        for (;;) {
                curs_dx = 0;
                curs_dy = 0;

                SDL_Event event;
                while (SDL_PollEvent(&event)) {
                        if (event.type == SDL_QUIT) {
                                goto _exit;
                        }

                        if (event.type == SDL_MOUSEMOTION) {
                                curs_dx = mouse_sens * (double)event.motion.xrel;
                                curs_dy = mouse_sens * (double)event.motion.yrel;
                        }
                }

                const uint8_t *keys = SDL_GetKeyboardState(NULL);

                cam_vel.x = keys[SDL_SCANCODE_A] * -0.04;
                cam_vel.y = keys[SDL_SCANCODE_LCTRL] * -0.04;
                cam_vel.z = keys[SDL_SCANCODE_W] * -0.04;

                cam_vel.x += keys[SDL_SCANCODE_D] * 0.04;
                cam_vel.y += keys[SDL_SCANCODE_LSHIFT] * 0.04;
                cam_vel.z += keys[SDL_SCANCODE_S] * 0.04;

                if (keys[SDL_SCANCODE_ESCAPE]) {
                        SDL_SetRelativeMouseMode(SDL_FALSE);
                }

                move_light_to(&light, VEC3(sin(i), 1.1, cos(i)));
#if 0
                light.colour = VEC4(
                        0.8f + 0.2*sin(8*i),
                        0.8f + 0.2*sin(8*i + 2.094),
                        0.8f + 0.2*sin(8*i + 4.188),
                        1.0
                );
#endif

                double cx = cos(-curs_dx);
                double cy = cos(-curs_dy);
                double sx = sin(-curs_dx);
                double sy = sin(-curs_dy);

                Mat3 cam_rot = {
                        cx,     0,   sx,
                        sy*sx,  cy, -sy*cx,
                        -cy*sx, sy,  cy*cx
                };
                move_cam(&cam, cam_rot, cam_vel);

                clock_t start = clock();

                clear_screen(0x21);

                // @TODO same thing as the shadow mapping - handle implicitly
                reset_dbuf();

                if (frames_drawn % 2) {
                        // @TODO this is a bit low-level for shadow map handling, tie it to moving the light source and organise everything in a Scene so it can be re-rendered to the shadow map
                        reset_smap(light.shadow_map);
                        //render_model_smap(floor_model, &light);
                        render_model_smap(main_model, &light);
                        //render_model_smap(ceiling_model, &light);
                }

                render_model(floor_model, &cam, &light);
                render_model(main_model, &cam, &light);
                //render_model(ceiling_model, &cam, &light);

                if ((frames_drawn % 64) == 0) {
                        frame_time = (double)(clock() - start) / CLOCKS_PER_SEC * 1000;
                        fps = 1000 / frame_time;
                        fps_high = fps > fps_high ? fps : fps_high;
                        fps_low = fps < fps_low ? fps : fps_low;
                        frame_time_high = frame_time > frame_time_high ? frame_time : frame_time_high;
                        frame_time_low = frame_time < frame_time_low ? frame_time : frame_time_low;
                }

                //fog(CLAMP(0, 1, 0.5*cam.pos.y), VEC4(0.8, 0.8, 0.8, 1));

                Vec3 light_pos_proj = persp(m4v3_mul(cam.view_mat, light.pos));
                Vec3 light_pos_screen = m4v3_mul(g_viewport, light_pos_proj);
                if (!out_of_view(light_pos_screen) && dbuf_read(light_pos_screen.x, light_pos_screen.y) > 1.0/light_pos_screen.z) {
                        point(light_pos_screen.x, light_pos_screen.y, light.colour);
                        render_text(light_pos_proj, light.colour, 0.1, fontset, "light %v", &light.pos);
                }

                render_axes(&cam);

                time_t rawtime;
                time(&rawtime);
                struct tm *timeinfo = localtime(&rawtime);

                render_text(
                        VEC3(-0.95, 0.9, 0), VEC4(0, 1, 0, 0.4), 0.15, fontset,

                        "-SRen-\n\n"

                        "%s\n"
                        "%dx%d @ %f fps [%f/%f]\n"
                        "%fms frame time [%f/%f]\n"
                        "%d frames drawn\n"
                        "main model: \"%s\" (%d faces)\n"
                        "texture: \"%s\" (%dx%d)\n"
                        "mem: %f MiB\n"
                        "camera @ %v\n",

                        asctime(timeinfo),
                        g_window_width, g_window_height,
                        fps, fps_high, fps_low,
                        frame_time, frame_time_high, frame_time_low,
                        frames_drawn,
                        argv[1], main_model->obj->face_count,
                        argv[2], main_model->texture->width, main_model->texture->height,
                        (double)arena_get_usage(&render_arena)/(1024*1024),
                        &cam.pos
                );

                SDL_UpdateWindowSurface(window);

                i += 0.005;

                ++frames_drawn;
        }

_exit:
        arena_deinit(&render_arena);
        SDL_Quit();
        return 0;
}
