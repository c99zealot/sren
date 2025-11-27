// @TODO shadow mapping, build a shadow map for each Light, use Scene object so the shadow map can be built for all models in one go
// @TODO also start using a Scene object to organise everything instead of just rendering it separately

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

#include <SDL2/SDL.h>

#include "sren.h"
#include "arena.h"

#define SWAP(T, x, y) { \
        T _tmp = x;     \
        (x) = (y);      \
        (y) = _tmp;     \
}

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

#define RGB(x, y, z)  (((x) << 16) | ((y) << 8) | (z))
#define RGBf(x, y, z) (((uint8_t)(x*255) << 16) | ((uint8_t)(y*255) << 8) | (uint8_t)(z*255))

#define TIME(x) {                                                   \
        clock_t start = clock();                                    \
        x;                                                          \
        printf("%s : %.4fms\n", #x,                                 \
                (double)(clock() - start) / CLOCKS_PER_SEC * 1000); \
}

SDL_Window  *g_window;
SDL_Surface *g_screen;
uint32_t *g_frame_buffer;
double *g_depth_buffer;

Light g_light = {
        .colour = {1.0, 1.0, 1.0},
        .intensity = 1.0
};

const int g_window_width = 800;
const int g_window_height = 800;

Mat4 g_viewport = {
        (g_window_width - 1)/2,   0.0,                   0.0, 0.0,
        0.0,                    (g_window_height/2 - 1), 0.0, 0.0,
        0.0,                    0.0,                     1.0, 0.0,
        0.0,                    0.0,                     0.0, 1.0
};

static inline void point(int x, int y, int colour) {
        g_frame_buffer[g_window_width * (g_window_height/2 - y) - g_window_width/2 + x] = colour;
}

static inline void dbuf_write(int x, int y, double depth) {
        g_depth_buffer[g_window_width * (g_window_height/2 - y) - g_window_width/2 + x] = depth;
}

static inline double dbuf_read(int x, int y) {
        return g_depth_buffer[g_window_width * (g_window_height/2 - y) - g_window_width/2 + x];
}

void line(int ax, int ay, int bx, int by, int colour) {
        int steep = abs(ay - by) > abs(ax - bx);
        if (steep) {
                SWAP(int, ax, ay);
                SWAP(int, bx, by);
        }

        if (ax > bx) {
                SWAP(int, ax, bx);
                SWAP(int, ay, by);
        }

        float y = ay;
        float m = (by - ay) / (float)(bx - ax);
        for (int x = ax; x <= bx; ++x) {
                if (steep) {
                        point(y, x, colour);
                } else {
                        point(x, y, colour);
                }
                y += m;
        }
}

int out_of_view(Vec3 v) {
        return (v.x > g_window_width/2 - 1) ||
               (v.x < -g_window_width/2) ||
               (v.y > g_window_height/2) ||
               (v.y < -g_window_height/2);
}

void render_face(Model *model, Camera *cam, Vec3 a, Vec3 b, Vec3 c, Face *f, Vec3 light_pos) {
        Vec3 *v = model->obj->vertices;
        Vec3 *vt = model->obj->uvs;
        Vec3 *vn = model->obj->norms;

        Vec3 ta = vt[f->v0[UV]];
        Vec3 tb = vt[f->v1[UV]];
        Vec3 tc = vt[f->v2[UV]];

        Vec3 na = m4v3_mul(cam->inv_tr, vn[f->v0[NORM]]);
        Vec3 nb = m4v3_mul(cam->inv_tr, vn[f->v1[NORM]]);
        Vec3 nc = m4v3_mul(cam->inv_tr, vn[f->v2[NORM]]);

        Vec3 a_to_light = vec3_sub(light_pos, a);
        Vec3 b_to_light = vec3_sub(light_pos, b);
        Vec3 c_to_light = vec3_sub(light_pos, c);

        double light_a = attenuate(vec3_dot(na, unit(a_to_light)), vec3_norm(a_to_light));
        double light_b = attenuate(vec3_dot(nb, unit(b_to_light)), vec3_norm(b_to_light));
        double light_c = attenuate(vec3_dot(nc, unit(c_to_light)), vec3_norm(c_to_light));

        a = m4v3_mul(g_viewport, a);
        b = m4v3_mul(g_viewport, b);
        c = m4v3_mul(g_viewport, c);

        if (out_of_view(a) || out_of_view(b) || out_of_view(c)) {
                return;
        }

        double recip_area = 1.0/signed_tri_area2(a, b, c);

        int min_x = MIN(MIN(a.x, b.x), c.x);
        int max_x = MAX(MAX(a.x, b.x), c.x);
        int min_y = MIN(MIN(a.y, b.y), c.y);
        int max_y = MAX(MAX(a.y, b.y), c.y);

        double Aa = b.y - c.y;
        double Ab = c.y - a.y;
        double Ac = a.y - b.y;
                                   
        double Ba = c.x - b.x;
        double Bb = a.x - c.x;
        double Bc = b.x - a.x;

        double Ca = b.x*c.y - c.x*b.y;
        double Cb = c.x*a.y - a.x*c.y;
        double Cc = a.x*b.y - b.x*a.y;
              
        double Ea = Aa*min_x + Ba*max_y + Ca;
        double Eb = Ab*min_x + Bb*max_y + Cb;
        double Ec = Ac*min_x + Bc*max_y + Cc;

        double Ea_x0 = Ea;
        double Eb_x0 = Eb;
        double Ec_x0 = Ec;

        for (int y = max_y; y >= min_y; --y) {
                for (int x = min_x; x <= max_x; ++x) {
                        if (Ea >= 0 && Eb >= 0 && Ec >= 0) {
                                Vec3 uv = vec3_scale(
                                        vec3_add(
                                                vec3_scale(ta, Ea),
                                                vec3_add(
                                                        vec3_scale(tb, Eb),
                                                        vec3_scale(tc, Ec)
                                                )
                                        ),
                                        recip_area
                                );

                                Vec3 colour = sample_texture(model->texture, uv.x, uv.y);
                                double interp_light = recip_area * (Ea*light_a + Eb*light_b + Ec*light_c);
                                colour = vec3_scale(colour, interp_light);
                                uint32_t lit_colour = RGBf(colour.x, colour.y, colour.z);

                                double depth = recip_area * (Ea*a.z + Eb*b.z + Ec*c.z);
                                if (depth > dbuf_read(x, y)) {
                                        point(x, y, lit_colour);
                                        dbuf_write(x, y, depth);
                                }
                        }

                        Ea += Aa;
                        Eb += Ab;
                        Ec += Ac;
                }

                Ea = Ea_x0 -= Ba;
                Eb = Eb_x0 -= Bb;
                Ec = Ec_x0 -= Bc;
        }
}

static inline Vec3 persp(Vec3 v, double f) {
        return vec3_scale(v, 1.0/(1 - v.z/f));
}

void render_model(Model *model, Camera *cam) {
        Obj *obj = model->obj;
        Vec3 *v = obj->vertices;
        Face *f = obj->faces;

        Vec3 light_in_view = m4v3_mul(cam->view_mat, g_light.pos);

        for (int i = 0; i < obj->face_count; ++i) {
                Vec3 a = v[f[i].v0[VERTEX]];
                Vec3 b = v[f[i].v1[VERTEX]];
                Vec3 c = v[f[i].v2[VERTEX]];

                a = m4v3_mul(cam->view_mat, a);
                b = m4v3_mul(cam->view_mat, b);
                c = m4v3_mul(cam->view_mat, c);

#ifdef PERSP
                double foc = vec3_norm(vec3_sub(cam->subject, cam->pos));
                a = persp(a, foc);
                b = persp(b, foc);
                c = persp(c, foc);
#endif

                Vec3 face_norm = cross(vec3_sub(b, a), vec3_sub(c, a));
                if (face_norm.z <= 0) {
                        continue;
                }

#ifdef WIREFRAME
                line(a.x, a.y, b.x, b.y, 0xff);
                line(b.x, b.y, c.x, c.y, 0xff);
                line(c.x, c.y, a.x, a.y, 0xff);
                continue;
#endif
                render_face(model, cam, a, b, c, &f[i], light_in_view);
        }
}

void reset_dbuf(void) {
        for (int i = 0; i < g_window_width*g_window_height; ++i) {
                g_depth_buffer[i] = -1024.0;
        }
}

int main(int argc, char **argv) {
        if (argc < 3) {
                fprintf(stderr, "usage: %s <obj> <texture>\n", argv[0]);
                return -1;
        }

        if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
                fprintf(stderr, "SDL_Init failure: %s\n", SDL_GetError());
                return -1;
        }

        g_window = SDL_CreateWindow("software renderer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, g_window_width, g_window_height, SDL_WINDOW_SHOWN);
        if (g_window == NULL) {
                fprintf(stderr, "SDL_CreateWindow failure: %s\n", SDL_GetError());
                SDL_Quit();
                return -1;
        }

        Arena render_arena;
        arena_init(&render_arena, ARENA_RESERVE_DEFAULT, 1 << 16, sizeof(double));
        g_depth_buffer = arena_alloc(&render_arena, g_window_width*g_window_height);

        Model *main_model = load_model(&render_arena, argv[1], argv[2], NULL, 1024, 1024, 0, 0);
        Model *floor_model = load_model(&render_arena, "assets/floor.obj", "assets/floor_texture.data", NULL, 1024, 1024, 0, 0);
        //Scene *scene = mkscene(&render_arena);

        g_screen = SDL_GetWindowSurface(g_window);
        if (g_screen == NULL) {
                fprintf(stderr, "SDL_GetWindowSurface failure: %s\n", SDL_GetError());
                SDL_Quit();
                return -1;
        }
        g_frame_buffer = g_screen->pixels;

        double i = 0.0;
        for (;;) {
                SDL_Event event;
                while (SDL_PollEvent(&event)) {
                        switch (event.type) {
                                case SDL_QUIT:
                                        goto _exit;
                        }
                }

                g_light.pos = (Vec3){sin(i), sin(i), cos(i)};

                Camera cam = {
                        .pos     = {cos(0.5*i), 0.5, sin(0.5*i)},
                        .subject = {0, 0, 0},
                        .up      = {0, 1, 0}
                };
                set_view(&cam);

                reset_dbuf();

                clock_t start = clock();

                render_model(floor_model, &cam);
                render_model(main_model, &cam);

                g_light.pos = m4v3_mul(cam.view_mat, g_light.pos);
                g_light.pos = m4v3_mul(g_viewport, g_light.pos);

                if (!out_of_view(g_light.pos) && dbuf_read(g_light.pos.x, g_light.pos.y) < g_light.pos.z) {
                        point(g_light.pos.x, g_light.pos.y, RGBf(1.0, 1.0, 1.0));
                }

                double elapsed_ms = (double)(clock() - start) / CLOCKS_PER_SEC * 1000;
                printf("frame time: %.4fms\nfps: %f\n\n", elapsed_ms, 1000 / elapsed_ms);

                SDL_UpdateWindowSurface(g_window);

                memset(g_frame_buffer, 0, g_window_width*g_window_height*sizeof(uint32_t));
                i += 0.005;
        }

_exit:
        arena_deinit(&render_arena);
        SDL_Quit();
        return 0;
}
