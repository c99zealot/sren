// @TODO start using a Scene object to organise everything instead of just rendering it separately
// @TODO make all the texture-y things Textures; depth buffer, shadow maps, etc. and write general sampling functions for 0,0 @ bottom left / 0,0 @ center
// @TODO handle asset loading failures with default assets
// @TODO phong lighting w/ normal maps
// @TODO serious OBJ parser
// @TODO make matrices structs so we can use =
// @TODO support more characters in render_text and render ? for unsupported characters, also prevent drawing glyphs outside of the framebuffer

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <stdarg.h>

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
uint32_t *g_framebuffer;
double *g_depth_buffer;

Light g_light = {
        .colour = VEC3(1.0, 1.0, 1.0),
        .subject = VEC3(0, 0, 0),
        .intensity = 1.0
};

const int g_window_width = 1000;
const int g_window_height = 1000;
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

enum {
        GLYPH_0,
        GLYPH_1,
        GLYPH_2,
        GLYPH_3,
        GLYPH_4,
        GLYPH_5,
        GLYPH_6,
        GLYPH_7,
        GLYPH_8,
        GLYPH_9,

        GLYPH_A = 'a',
        GLYPH_B = 'b',
        GLYPH_C = 'c',
        GLYPH_D = 'd',
        GLYPH_E = 'e',
        GLYPH_F = 'f',
        GLYPH_G = 'g',
        GLYPH_H = 'h',
        GLYPH_I = 'i',
        GLYPH_J = 'j',
        GLYPH_K = 'k',
        GLYPH_L = 'l',
        GLYPH_M = 'm',
        GLYPH_N = 'n',
        GLYPH_O = 'o',
        GLYPH_P = 'p',
        GLYPH_Q = 'q',
        GLYPH_R = 'r',
        GLYPH_S = 's',
        GLYPH_T = 't',
        GLYPH_U = 'u',
        GLYPH_V = 'v',
        GLYPH_W = 'w',
        GLYPH_X = 'x',
        GLYPH_Y = 'y',
        GLYPH_Z = 'z',

        GLYPH_DOT        = '.',
        GLYPH_SPACE      = ' ',
        GLYPH_MINUS      = '-',
        GLYPH_COMMA      = ',',
        GLYPH_RBRACKET   = ')',
        GLYPH_LBRACKET   = '(',
        GLYPH_QUESTION   = '?',
        GLYPH_UNDERSCORE = '_',
};

int g_glyph_seqs[][10] = {
        [GLYPH_0] = {9, 3, 5, 11, 9, 5,          -1},
        [GLYPH_1] = {6, 4, 10, 9, 11,            -1},
        [GLYPH_2] = {3, 5, 8, 6, 9, 11,          -1},
        [GLYPH_3] = {3, 5, 8, 6, 8, 11, 9,       -1},
        [GLYPH_4] = {11, 5, 6, 8,                -1},
        [GLYPH_5] = {5, 3, 6, 8, 11, 9,          -1},
        [GLYPH_6] = {5, 6, 8, 11, 9, 6,          -1},
        [GLYPH_7] = {3, 5, 9,                    -1},
        [GLYPH_8] = {8, 6, 3, 5, 11, 9, 6,       -1},
        [GLYPH_9] = {5, 6, 3, 5, 8, 6, 8, 11, 9, -1},

        [GLYPH_A] = {3, 5, 11, 9, 6, 8, -1},
        [GLYPH_B] = {3, 9, 11, 8, 6,    -1},
        [GLYPH_C] = {11, 9, 3, 5,       -1},
        [GLYPH_D] = {5, 11, 9, 6, 8,    -1},
        [GLYPH_E] = {11, 9, 3, 5, 6,    -1},
        [GLYPH_F] = {5, 3, 6, 8, 6, 9,  -1},
        [GLYPH_G] = {9, 11,5, 3, 6, 8,  -1},
        [GLYPH_H] = {3, 9, 6, 8, 11,    -1},
        [GLYPH_I] = {4, 10,             -1},
        [GLYPH_J] = {4, 5, 11, 9,       -1},
        [GLYPH_K] = {3, 9, 6, 11, 6, 5, -1},
        [GLYPH_L] = {3, 9, 11,          -1},
        [GLYPH_M] = {9, 3, 7, 5, 11,    -1},
        [GLYPH_N] = {9, 3, 11, 5,       -1},
        [GLYPH_O] = {9, 3, 5, 11, 9,    -1},
        [GLYPH_P] = {9, 3, 5, 8, 6,     -1},
        [GLYPH_Q] = {8, 10, 4, 3, 6, 7, -1},
        [GLYPH_R] = {9, 6, 5,           -1},
        [GLYPH_S] = {5, 3, 6, 8, 11, 9, -1},
        [GLYPH_T] = {7, 6, 3, 9, 11,    -1},
        [GLYPH_U] = {3, 9, 11, 5,       -1},
        [GLYPH_V] = {3, 10, 5,          -1},
        [GLYPH_W] = {3, 9, 7, 11, 5,    -1},
        [GLYPH_X] = {9, 5, 7, 3, 11,    -1},
        [GLYPH_Y] = {3, 7, 5, 7, 10,    -1},
        [GLYPH_Z] = {3, 5, 9, 11,       -1},

        [GLYPH_SPACE]    = {0,           -1},
        [GLYPH_MINUS]    = {6, 8,        -1},
        [GLYPH_COMMA]    = {7, 9,        -1},
        [GLYPH_QUESTION] = {0,           -1},
        [GLYPH_LBRACKET] = {2, 3, 6, 11, -1},
        [GLYPH_RBRACKET] = {0, 5, 8, 9,  -1},
        [GLYPH_DOT]      = {9, 9,        -1},
};

//
// point - draws a point at (x, y) in the frame buffer
//
static inline void point(int x, int y, int colour) {
        g_framebuffer[SCREEN(x, y)] = colour;
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
        if (a.x == b.x && a.y == b.y) {
                point(a.x, a.y, colour);
                return;
        }

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
// render_glyph - renders a glyph to the framebuffer
//
Vec3 render_glyph(char glyph_id, Vec3 pos, double unit) {
        Vec3 points[] = {
                pos,
                VEC3(pos.x + unit, pos.y, 0),
                VEC3(pos.x + 2*unit, pos.y, 0),

                VEC3(pos.x, pos.y - unit, 0),
                VEC3(pos.x + unit, pos.y - unit, 0),
                VEC3(pos.x + 2*unit, pos.y - unit, 0),

                VEC3(pos.x, pos.y - 2*unit, 0),
                VEC3(pos.x + unit, pos.y - 2*unit, 0),
                VEC3(pos.x + 2*unit, pos.y - 2*unit, 0),

                VEC3(pos.x, pos.y - 3*unit, 0),
                VEC3(pos.x + unit, pos.y - 3*unit, 0),
                VEC3(pos.x + 2*unit, pos.y - 3*unit, 0),
        };

        for (int i = 0; i < 12; ++i) {
                points[i] = m4v3_mul(g_viewport, points[i]);
        }

        Vec3 line_start = points[g_glyph_seqs[glyph_id][0]];
        for (int i = 1; g_glyph_seqs[glyph_id][i] != -1; ++i) {
                line(line_start, points[g_glyph_seqs[glyph_id][i]], 0x00FF00);
                line_start = points[g_glyph_seqs[glyph_id][i]];
        }

        pos.x += ((glyph_id == GLYPH_DOT) ? 1 : 3) * unit;
        return pos;
}

//
// render_text_int - called by render_text to render integers as text to the framebuffer
//
Vec3 render_text_int(Vec3 pos, double unit, int d) {
        if (d < 0) {
                pos = render_glyph(GLYPH_MINUS, pos, unit);
                d = -d;
        }

        int place_val = 1;
        while (place_val < d) {
                place_val *= 10;
        }

        if (place_val > d) {
                place_val /= 10;
        }

        while (place_val >= 1) {
                pos = render_glyph((d / place_val) % 10, pos, unit);
                place_val /= 10;
        }

        return pos;
}

//
// render_text_double - called by render_text to render doubles as text to the framebufer
//
Vec3 render_text_double(Vec3 pos, double unit, double f) {
        if (f < 0) {
                pos = render_glyph(GLYPH_MINUS, pos, unit);
                f = -f;
        }

        long whole = (long)f;

        if (whole == 0) {
                pos = render_glyph(GLYPH_0, pos, unit);
        }

        int place_val = 1;
        while (place_val < whole) {
                place_val *= 10;
        }

        if (place_val > whole) {
                place_val /= 10;
        }

        while (place_val >= 1) {
                pos = render_glyph((whole / place_val) % 10, pos, unit);
                place_val /= 10;
        }

        pos = render_glyph(GLYPH_DOT, pos, unit);

        int max_digits = 2;
        do {
                f *= 10;
                pos = render_glyph((long)f % 10, pos, unit);
        } while (f > (long)f && max_digits-- > 1);

        return pos;
}

//
// render_text_vector - called by render_text to render a vector as text to the framebuffer
//
Vec3 render_text_vector(Vec3 pos, double unit, Vec3 v) {
        pos = render_glyph(GLYPH_LBRACKET, pos, unit);
        pos = render_text_double(pos, unit, v.x);
        pos = render_glyph(GLYPH_COMMA, pos, unit);
        pos = render_text_double(pos, unit, v.y);
        pos = render_glyph(GLYPH_COMMA, pos, unit);
        pos = render_text_double(pos, unit, v.z);
        pos = render_glyph(GLYPH_RBRACKET, pos, unit);
}

//
// render_text - renders formatted text to the framebuffer
//
Vec3 render_text(Vec3 pos, double unit, const char *fmt, ...) {
        va_list arglist;
        va_start(arglist, fmt);

        Vec3 curs = pos;
        for (int i = 0; fmt[i] != '\0'; ++i) {
                if (fmt[i] == '\n') {
                        curs.x = pos.x;
                        curs.y -= 4*unit;
                } else if (fmt[i] == '%') {
                        ++i;
                        if (fmt[i] == 'd') {
                                curs = render_text_int(curs, unit, va_arg(arglist, int));
                        } else if (fmt[i] == 'f') {
                                curs = render_text_double(curs, unit, va_arg(arglist, double));
                        } else if (fmt[i] == 'v') {
                                curs = render_text_vector(curs, unit, va_arg(arglist, Vec3));
                        } else if (fmt[i] == 's') {
                                curs = render_text(curs, unit, va_arg(arglist, char*));
                        }
                } else {
                        curs = render_glyph(fmt[i], curs, unit);
                }
        }

        va_end(arglist);

        return curs;
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
void render_face(Face *f, Model *model, Camera *cam, Light *light) {
        Vec3 *v = model->obj->vertices;
        Vec3 *vt = model->obj->uvs;
        Vec3 *vn = model->obj->norms;

        Vec3 a = v[f->v0[VERTEX]];
        Vec3 b = v[f->v1[VERTEX]];
        Vec3 c = v[f->v2[VERTEX]];

        Vec3 a_lightview = m4v3_mul(g_light.view_mat, a);
        Vec3 b_lightview = m4v3_mul(g_light.view_mat, b);
        Vec3 c_lightview = m4v3_mul(g_light.view_mat, c);

        a = m4v3_mul(cam->view_mat, a);
        b = m4v3_mul(cam->view_mat, b);
        c = m4v3_mul(cam->view_mat, c);

        if (a.z > 0 || b.z > 0 || c.z > 0) {
                return;
        }

        a = persp(a);
        b = persp(b);
        c = persp(c);

        Vec3 face_norm = cross(vec3_sub(b, a), vec3_sub(c, a));
        if (face_norm.z <= 0) {
                return;
        }

        double a_recip_z = 1.0/a.z;
        double b_recip_z = 1.0/b.z;
        double c_recip_z = 1.0/c.z;

        Vec3 ta = vec3_scale(vt[f->v0[UV]], a_recip_z);
        Vec3 tb = vec3_scale(vt[f->v1[UV]], b_recip_z);
        Vec3 tc = vec3_scale(vt[f->v2[UV]], c_recip_z);

        Vec3 na = m4v3_mul(cam->inv_tr, vn[f->v0[NORM]]);
        Vec3 nb = m4v3_mul(cam->inv_tr, vn[f->v1[NORM]]);
        Vec3 nc = m4v3_mul(cam->inv_tr, vn[f->v2[NORM]]);

        Vec3 light_in_view = m4v3_mul(cam->view_mat, g_light.pos);

        Vec3 a_to_light = vec3_sub(light_in_view, a);
        Vec3 b_to_light = vec3_sub(light_in_view, b);
        Vec3 c_to_light = vec3_sub(light_in_view, c);

        Vec3 a_light_proj = vec3_scale(a_lightview, a_recip_z);
        Vec3 b_light_proj = vec3_scale(b_lightview, b_recip_z);
        Vec3 c_light_proj = vec3_scale(c_lightview, c_recip_z);

        double light_a = attenuate(vec3_dot(na, unit(a_to_light)), vec3_norm(a_to_light));
        double light_b = attenuate(vec3_dot(nb, unit(b_to_light)), vec3_norm(b_to_light));
        double light_c = attenuate(vec3_dot(nc, unit(c_to_light)), vec3_norm(c_to_light));

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

        #define INTERP(a, b, c)     (recip_area * (Ea*(a) + Eb*(b) + Ec*(c)))
        #define VEC_INTERP(a, b, c) (vec3_scale(vec3_add(vec3_scale(a, Ea), vec3_add(vec3_scale(b, Eb), \
                                        vec3_scale(c, Ec))), recip_area))

        for (int y = max_y; y >= min_y; --y) {
                for (int x = min_x; x <= max_x; ++x) {
                        if (Ea >= 0 && Eb >= 0 && Ec >= 0) {
                                double interp_recip_z = INTERP(a_recip_z, b_recip_z, c_recip_z);

                                double depth = -interp_recip_z;
                                if (depth <= dbuf_read(x, y)) {
                                        continue;
                                }

                                Vec3 uv = VEC_INTERP(ta, tb, tc);
                                uv = vec3_scale(uv, 1.0/interp_recip_z);

                                Vec3 interp_light_proj = VEC_INTERP(a_light_proj, b_light_proj, c_light_proj);
                                interp_light_proj = vec3_scale(interp_light_proj, 1.0/interp_recip_z);
                                interp_light_proj = m4v3_mul(g_viewport, persp(interp_light_proj));

                                double light_pixel = 0.1;
                                if (smap_read(light->shadow_map, interp_light_proj.x, interp_light_proj.y) <= interp_light_proj.z + 0.05) {
                                        light_pixel = INTERP(light_a, light_b, light_c);
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

        #undef INTERP
        #undef VEC_INTERP
}

//
// render_face_smap - renders a single face in a model to a shadow map
//
void render_face_smap(Face *f, Model *model, Light *light) {
        Vec3 *v = model->obj->vertices;

        Vec3 a = persp(m4v3_mul(light->view_mat, v[f->v0[VERTEX]]));
        Vec3 b = persp(m4v3_mul(light->view_mat, v[f->v1[VERTEX]]));
        Vec3 c = persp(m4v3_mul(light->view_mat, v[f->v2[VERTEX]]));

        if (a.z > 0 || b.z > 0 || c.z > 0) {
                return;
        }

        Vec3 face_norm = cross(vec3_sub(b, a), vec3_sub(c, a));
        if (face_norm.z <= 0) {
                return;
        }

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
        Face *f = obj->faces;

        for (int i = 0; i < obj->face_count; ++i) {
                render_face(&f[i], model, cam, &g_light);
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
                render_face_smap(&f[i], model, light);
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

#if 0
void draw_smap(Shadow_Map *smap, uint8_t *pixels) {
        for (int i = 0; i < smap->width*smap->height; ++i) {
                double d = (smap->data[i] + 4) / 8;
                uint32_t colour = RGBf(d, d, d);
                memcpy(pixels, &colour, 4);
                pixels += 4;
        }
}
#endif

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

        line(o, x, 0xff0000);
        line(o, y, 0x00ff00);
        line(o, z, 0x0000ff);
}

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
        g_framebuffer = g_screen->pixels;

        g_light.shadow_map = mk_smap(&render_arena, g_window_width, g_window_height);

        Camera cam = {
                .pos = VEC3(0, 0, 2),
                .subject = VEC3(0, 0, 1),
                .up = VEC3(0, 1, 0)
        };
        init_cam(&cam);

        SDL_WarpMouseInWindow(g_window, g_window_width/2, g_window_height/2);
        SDL_SetRelativeMouseMode(SDL_TRUE);

        const double mouse_sens = 0.01;

        double curs_dx = 0;
        double curs_dy = 0;

        Vec3 cam_vel = VEC3(0, 0, 0);

        int frames_drawn = 0;
        double fps = 0;
        double elapsed_ms = 0;
        int dig_counter = 0;
        int letter_counter = 0;

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

                cam_vel.x = keys[SDL_SCANCODE_A] * -0.02;
                cam_vel.y = keys[SDL_SCANCODE_LCTRL] * -0.02;
                cam_vel.z = keys[SDL_SCANCODE_W] * -0.02;

                cam_vel.x += keys[SDL_SCANCODE_D] * 0.02;
                cam_vel.y += keys[SDL_SCANCODE_LSHIFT] * 0.02;
                cam_vel.z += keys[SDL_SCANCODE_S] * 0.02;

                if (keys[SDL_SCANCODE_ESCAPE]) {
                        SDL_SetRelativeMouseMode(SDL_FALSE);
                }

                g_light.pos = VEC3(1.5*sin(i), 1, 1.5*cos(i));
                g_light.subject = VEC3(0, 0, 0);
                g_light.up = VEC3(0, 1, 0);
                set_light_view(&g_light);
                reset_smap(&g_light);

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

                reset_dbuf();

                clock_t start = clock();

                render_model_smap(floor_model, &g_light);
                render_model_smap(main_model, &g_light);

#if 0
                draw_smap(g_light.shadow_map, smap_screen->pixels);
#endif

                render_model(floor_model, &cam);
                render_model(main_model, &cam);

                if (frames_drawn % 64 == 0) {
                        elapsed_ms = (double)(clock() - start) / CLOCKS_PER_SEC * 1000;
                        fps = 1000 / elapsed_ms;
                }

                render_axes(&cam);

                g_light.pos = m4v3_mul(g_viewport, persp(m4v3_mul(cam.view_mat, g_light.pos)));

                if (!out_of_view(g_light.pos) && dbuf_read(g_light.pos.x, g_light.pos.y) < g_light.pos.z) {
                        point(g_light.pos.x, g_light.pos.y, RGBf(1.0, 1.0, 1.0));
                }

                render_text(VEC3(-0.95, 0.95, 0), 0.012, "-sren-\n%f fps\n%fms frame time\ncamera at %v", fps, elapsed_ms, cam.pos);

                SDL_UpdateWindowSurface(g_window);
#if 0
                SDL_UpdateWindowSurface(smap_window);
#endif

                memset(g_framebuffer, 0, g_window_width*g_window_height*sizeof(uint32_t));
                i += 0.01;

                ++frames_drawn;
        }

_exit:
        arena_deinit(&render_arena);
        SDL_Quit();
        return 0;
}
