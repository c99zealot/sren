// @TODO shadow mapping, build a shadow map for each Light, use Scene object so the shadow map can be built for all models in one go
// @TODO make all the texture-y things Textures; depth buffer, shadow maps, etc. and write a general sampling functions for 0,0 @ bottom left / 0,0 @ center
// @TODO also start using a Scene object to organise everything instead of just rendering it separately
// @TODO handle asset loading failures with default assets

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
#define MAX3(x, y, z) (MAX(MAX(x, y), z))
#define MIN3(x, y, z) (MIN(MIN(x, y), z))

#define RGB(x, y, z)  (((x) << 16) | ((y) << 8) | (z))
#define RGBf(x, y, z) (((uint8_t)(x*255) << 16) | ((uint8_t)(y*255) << 8) | (uint8_t)(z*255))

#define SCREEN(x, y) (g_window_width * (g_window_height/2 - y) - g_window_width/2 + x)

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
        .colour = VEC3(1.0, 1.0, 1.0),
        .subject = VEC3(0, 0, 0),
        .intensity = 1.0
};

const int g_window_width = 600;
const int g_window_height = 600;
const int g_min_x = -g_window_width  / 2;
const int g_min_y = -g_window_height / 2;
const int g_max_x =  g_window_width  / 2 - 1;
const int g_max_y =  g_window_height / 2 - 1;

Mat4 g_viewport = {
        (g_window_width - 1)/2, 0.0,                     0.0, 0.0,
        0.0,                    (g_window_height/2 - 1), 0.0, 0.0,
        0.0,                    0.0,                     1.0, 0.0,
        0.0,                    0.0,                     0.0, 1.0
};

//
// point - draws a point at (x, y) in the frame buffer
//
static inline void point(int x, int y, int colour) {
        g_frame_buffer[SCREEN(x, y)] = colour;
}

//
// dbuf_write - writes a depth value at (x, y) in the depth buffer
//
static inline void dbuf_write(int x, int y, double depth) {
        g_depth_buffer[SCREEN(x, y)] = depth;
}

//
// dbuf_read - reads the depth value at (x, y) in the depth buffer
//
static inline double dbuf_read(int x, int y) {
        return g_depth_buffer[SCREEN(x, y)];
}

//
// smap_write - writes a depth value at (x, y) in a shadow map
//
static inline void smap_write(Shadow_Map *shadow_map, int x, int y, double depth) {
        size_t i = (shadow_map->width * (shadow_map->height/2 - y) - shadow_map->width/2 + x);
        shadow_map->data[i] = depth;
}

//
// smap_read - reads the depth value at (x, y) in a shadow map
//
static inline double smap_read(Shadow_Map *shadow_map, int x, int y) {
        size_t i = (shadow_map->width * (shadow_map->height/2 - y) - shadow_map->width/2 + x);
        return shadow_map->data[i];
}

//
// line - draws a line in the frame buffer from (ax, ay) to (bx, by)
//
void line(Vec3 a, Vec3 b, int colour) {
        int steep = abs(a.y - b.y) > abs(a.x - b.x);
        if (steep) {
                SWAP(int, a.x, a.y);
                SWAP(int, b.x, b.y);
        }

        if (a.x > b.x) {
                SWAP(int, a.x, b.x);
                SWAP(int, a.y, b.y);
        }

        float _y = a.y;
        float m = (b.y - a.y) / (float)(b.x - a.x);
        for (int _x = a.x; _x <= b.x; ++_x) {
                if (steep) {
                        point(_y, _x, colour);
                } else {
                        point(_x, _y, colour);
                }
                _y += m;
        }
}

//
// persp - basic perspective projection
//
static inline Vec3 persp(Vec3 v) {
        double recip_z = -1.0/v.z;
        return VEC3(v.x * recip_z, v.y * recip_z, v.z);
}

//
// out_of_view - checks whether a point resides outside of the screen
//
int out_of_view(Vec3 v) {
        return (v.x > g_window_width/2 - 1) ||
               (v.x < -g_window_width/2) ||
               (v.y > g_window_height/2 - 1) ||
               (v.y < -g_window_height/2);
}

//
// render_face - renders a single face in a model to the framebuffer, performs texture mapping & Gouraud shading
//
void render_face(Model *model, Camera *cam, Vec3 a, Vec3 b, Vec3 c, Vec3 la, Vec3 lb, Vec3 lc, Face *f, Vec3 light_pos, Light *light) {
        Vec3 *vt = model->obj->uvs;
        Vec3 *vn = model->obj->norms;

        Vec3 ta = vt[f->v0[UV]];
        Vec3 tb = vt[f->v1[UV]];
        Vec3 tc = vt[f->v2[UV]];

        double a_recip_z = 1.0/a.z;
        double b_recip_z = 1.0/b.z;
        double c_recip_z = 1.0/c.z;

        Vec3 ta_proj = vec3_scale(ta, a_recip_z);
        Vec3 tb_proj = vec3_scale(tb, b_recip_z);
        Vec3 tc_proj = vec3_scale(tc, c_recip_z);

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

        Vec3 a_light_proj = vec3_scale(la, a_recip_z);
        Vec3 b_light_proj = vec3_scale(lb, b_recip_z);
        Vec3 c_light_proj = vec3_scale(lc, c_recip_z);

        double recip_area = 1.0/signed_tri_area2(a, b, c);

        int min_x = MIN3(a.x, b.x, c.x);
        int max_x = MAX3(a.x, b.x, c.x);
        int min_y = MIN3(a.y, b.y, c.y);
        int max_y = MAX3(a.y, b.y, c.y);

        min_x = MAX(min_x, -g_window_width/2);
        min_y = MAX(min_y, -g_window_height/2);
        max_x = MIN(max_x, g_window_width/2 - 1);
        max_y = MIN(max_y, g_window_height/2 - 1);

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
                                double interp_recip_z = recip_area * (Ea*a_recip_z + Eb*b_recip_z + Ec*c_recip_z);

                                double depth = -interp_recip_z;
                                if (depth <= dbuf_read(x, y)) {
                                        continue;
                                }

                                Vec3 uv = vec3_scale(
                                        vec3_add(
                                                vec3_scale(ta_proj, Ea),
                                                vec3_add(
                                                        vec3_scale(tb_proj, Eb),
                                                        vec3_scale(tc_proj, Ec)
                                                )
                                        ),
                                        recip_area
                                );
                                uv = vec3_scale(uv, 1.0/interp_recip_z);

                                Vec3 interp_light_proj = vec3_scale(
                                        vec3_add(
                                                vec3_scale(a_light_proj, Ea),
                                                vec3_add(
                                                        vec3_scale(b_light_proj, Eb),
                                                        vec3_scale(c_light_proj, Ec)
                                                )
                                        ),
                                        recip_area
                                );
                                interp_light_proj = vec3_scale(interp_light_proj, 1.0/interp_recip_z);
                                interp_light_proj = m4v3_mul(g_viewport, persp(interp_light_proj));

                                double light_pixel = 0.1;
                                if (smap_read(light->shadow_map, interp_light_proj.x, interp_light_proj.y) <= interp_light_proj.z + 0.05) {
                                        light_pixel = recip_area * (Ea*light_a + Eb*light_b + Ec*light_c);
                                }

                                Vec3 colour = sample_texture(model->texture, uv.x, uv.y);
                                colour = vec3_scale(colour, light_pixel);
                                uint32_t lit_colour = RGBf(colour.x, colour.y, colour.z);

                                point(x, y, lit_colour);
                                dbuf_write(x, y, depth);
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

//
// render_smap_face - renders a single face in a model to a shadow map
//
void render_smap_face(Light *light, Vec3 a, Vec3 b, Vec3 c) {
        a = m4v3_mul(g_viewport, a);
        b = m4v3_mul(g_viewport, b);
        c = m4v3_mul(g_viewport, c);

        double recip_area = 1.0/signed_tri_area2(a, b, c);

        int min_x = MIN3(a.x, b.x, c.x);
        int max_x = MAX3(a.x, b.x, c.x);
        int min_y = MIN3(a.y, b.y, c.y);
        int max_y = MAX3(a.y, b.y, c.y);

        min_x = MAX(min_x, g_min_x);
        min_y = MAX(min_y, g_min_y);
        max_x = MIN(max_x, g_max_x);
        max_y = MIN(max_y, g_max_y);

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
                                double depth = recip_area * (Ea*a.z + Eb*b.z + Ec*c.z);
                                if (depth > smap_read(light->shadow_map, x, y)) {
                                        smap_write(light->shadow_map, x, y, depth);
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

//
// render_model - renders all the faces in a model
//
void render_model(Model *model, Camera *cam) {
        Obj *obj = model->obj;
        Vec3 *v = obj->vertices;
        Face *f = obj->faces;

        Vec3 light_in_view = m4v3_mul(cam->view_mat, g_light.pos);

        for (int i = 0; i < obj->face_count; ++i) {
                Vec3 a = v[f[i].v0[VERTEX]];
                Vec3 b = v[f[i].v1[VERTEX]];
                Vec3 c = v[f[i].v2[VERTEX]];

                Vec3 la = m4v3_mul(g_light.view_mat, a);
                Vec3 lb = m4v3_mul(g_light.view_mat, b);
                Vec3 lc = m4v3_mul(g_light.view_mat, c);

                a = m4v3_mul(cam->view_mat, a);
                b = m4v3_mul(cam->view_mat, b);
                c = m4v3_mul(cam->view_mat, c);

                if (a.z > 0 || b.z > 0 || c.z > 0) {
                        continue;
                }

                a = persp(a);
                b = persp(b);
                c = persp(c);

                Vec3 face_norm = cross(vec3_sub(b, a), vec3_sub(c, a));
                if (face_norm.z <= 0) {
                        continue;
                }

#ifdef WIREFRAME
                line(a, b, 0xff);
                line(b, c, 0xff);
                line(c, a, 0xff);
#else
                render_face(model, cam, a, b, c, la, lb, lc, &f[i], light_in_view, &g_light);
#endif
        }
}

//
// render_model_smap - renders depths of a model's faces to a shadow map
//
void render_model_smap(Model *model, Light *light) {
        Obj *obj = model->obj;
        Vec3 *v = obj->vertices;
        Face *f = obj->faces;

        for (int i = 0; i < obj->face_count; ++i) {
                Vec3 a = persp(m4v3_mul(light->view_mat, v[f[i].v0[VERTEX]]));
                Vec3 b = persp(m4v3_mul(light->view_mat, v[f[i].v1[VERTEX]]));
                Vec3 c = persp(m4v3_mul(light->view_mat, v[f[i].v2[VERTEX]]));

                if (a.z > 0 || b.z > 0 || c.z > 0) {
                        continue;
                }

                Vec3 face_norm = cross(vec3_sub(b, a), vec3_sub(c, a));
                if (face_norm.z <= 0) {
                        continue;
                }

                render_smap_face(light, a, b, c);
        }
}

//
// reset_dbuf - resets the depth buffer to the maximum depth value
//
void reset_dbuf(void) {
        for (int i = 0; i < g_window_width*g_window_height; ++i) {
                g_depth_buffer[i] = -1024.0;
        }
}

//
// reset_smap - resets a Light's shadow map to the maximum depth value
//
void reset_smap(Light *light) {
        Shadow_Map *smap = light->shadow_map;
        for (int i = 0; i < smap->width*smap->height; ++i) {
                smap->data[i] = -1024.0;
        }
}

void draw_smap(Shadow_Map *smap, uint8_t *pixels) {
        for (int i = 0; i < smap->width*smap->height; ++i) {
                double d = (smap->data[i] + 4) / 8;
                uint32_t colour = RGBf(d, d, d);
                memcpy(pixels, &colour, 4);
                pixels += 4;
        }
}

void render_axes(Camera *cam) {
        Vec3 o = VEC3(0, 0, 0), x = VEC3(0.1, 0, 0), y = VEC3(0, 0.1, 0), z = VEC3(0, 0, 0.1);

        o = m4v3_mul(g_viewport, persp(m4v3_mul(cam->view_mat, o)));
        x = m4v3_mul(g_viewport, persp(m4v3_mul(cam->view_mat, x)));
        y = m4v3_mul(g_viewport, persp(m4v3_mul(cam->view_mat, y)));
        z = m4v3_mul(g_viewport, persp(m4v3_mul(cam->view_mat, z)));

        line(o, x, 0xff0000);
        line(o, y, 0x00ff00);
        line(o, z, 0x0000ff);
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

        g_window = SDL_CreateWindow("sren", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, g_window_width, g_window_height, SDL_WINDOW_SHOWN);
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

        g_light.shadow_map = mk_smap(&render_arena, g_window_width, g_window_height);

#if 0
        SDL_Window *smap_window = SDL_CreateWindow("sren_smap", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, g_window_width, g_window_height, SDL_WINDOW_SHOWN);
        if (smap_window == NULL) {
                fprintf(stderr, "SDL_CreateWindow failure: %s\n", SDL_GetError());
                SDL_Quit();
                return -1;
        }

        SDL_Surface *smap_screen = SDL_GetWindowSurface(smap_window);
        if (smap_screen == NULL) {
                fprintf(stderr, "SDL_GetWindowSurface failure: %s\n", SDL_GetError());
                SDL_Quit();
                return -1;
        }
#endif

        double i = 0.0;
        for (;;) {
                SDL_Event event;
                while (SDL_PollEvent(&event)) {
                        switch (event.type) {
                                case SDL_QUIT:
                                        goto _exit;
                        }
                }

                g_light.pos = VEC3(0.4*sin(i), 1.5, 0.4*cos(i));
                g_light.subject = VEC3(0, 0, 0);
                g_light.up = VEC3(0, 1, 0);
                set_light_view(&g_light);
                reset_smap(&g_light);

                Camera cam = {
                        .pos     = VEC3(1.5*cos(0.4*i), 0.8*sin(0.8*i), 1.5*sin(0.4*i)),
                        .subject = VEC3(0, 0, 0),
                        .up      = VEC3(0, 1, 0)
                };
                set_view(&cam);

                reset_dbuf();

                clock_t start = clock();

                render_model_smap(floor_model, &g_light);
                render_model_smap(main_model, &g_light);

#if 0
                draw_smap(g_light.shadow_map, smap_screen->pixels);
#endif

                render_model(floor_model, &cam);
                render_model(main_model, &cam);

                render_axes(&cam);

                g_light.pos = m4v3_mul(g_viewport, persp(m4v3_mul(cam.view_mat, g_light.pos)));

                if (!out_of_view(g_light.pos) && dbuf_read(g_light.pos.x, g_light.pos.y) < g_light.pos.z) {
                        point(g_light.pos.x, g_light.pos.y, RGBf(1.0, 1.0, 1.0));
                }

                double elapsed_ms = (double)(clock() - start) / CLOCKS_PER_SEC * 1000;
                printf(" %f fps [frame time: %.4fms]\t\r", 1000 / elapsed_ms, elapsed_ms);
                fflush(stdout);

                SDL_UpdateWindowSurface(g_window);
#if 0
                SDL_UpdateWindowSurface(smap_window);
#endif

                memset(g_frame_buffer, 0, g_window_width*g_window_height*sizeof(uint32_t));
                i += 0.01;
        }

_exit:
        arena_deinit(&render_arena);
        SDL_Quit();
        return 0;
}
