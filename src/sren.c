#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>

#include "sren.h"
#include "arena.h"
#include "error_assets.h"

#define SCREEN(x, y) (g_fb_width * (g_fb_height/2 - y) - g_fb_width/2 + x)

static uint32_t *g_framebuffer;
static double *g_depthbuffer;
static int g_fb_width;
static int g_fb_height;
static int g_min_x;
static int g_min_y;
static int g_max_x;
static int g_max_y;

//
// signed_tri_area2 - calculates twice the signed area of the triangle abc
//
static inline double signed_tri_area2(Vec3 a, Vec3 b, Vec3 c) {
        return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
}

Mat4 g_viewport; // @XXX exposing this for now...

/*
 * @XXX shadow mapping @XXX
 */

Shadow_Map *mk_smap(Arena *arena, size_t w, size_t h) {
        Shadow_Map *smap = arena_alloc(arena, sizeof(Shadow_Map));
        smap->width = w;
        smap->height = h;
        smap->data = arena_alloc(arena, w*h*sizeof(double)*4);

        return smap;
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
double smap_read(Shadow_Map *smap, int x, int y) {
        size_t i = (smap->width * (smap->height/2 - y) - smap->width/2 + x);
        return (i >= smap->width*smap->height) ? 0 : smap->data[i];
}

//
// render_face_smap - renders a single face in a model to a shadow map
//
static void render_face_smap(Face *f, Model *model, Light *light) {
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
// reset_smap - resets a Light's shadow map to the maximum depth value
//
void reset_smap(Shadow_Map *smap) {
        for (int i = 0; i < smap->width*smap->height; ++i) {
                smap->data[i] = -1024.0;
        }
}

/*
 * @XXX depth buffer @XXX 
 */

//
// dbuf_write - writes a depth value at (x, y) in the depth buffer
//
static inline void dbuf_write(int x, int y, double depth) {
        g_depthbuffer[SCREEN(x, y)] = depth;
}

//
// dbuf_read - reads the depth value at (x, y) in the depth buffer
//
double dbuf_read(int x, int y) {
        return g_depthbuffer[SCREEN(x, y)];
}


//
// reset_dbuf - resets the depth buffer to the maximum depth value
//
void reset_dbuf(void) {
        for (int i = 0; i < g_fb_width*g_fb_height; ++i) {
                g_depthbuffer[i] = -1024.0;
        }
}

/*
 * @XXX rendering @XXX
 */

//
// init_renderer - initialises internal data, ready to render to the framebuffer
//
void init_renderer(Arena *arena, uint32_t *framebuffer, size_t fb_width, size_t fb_height) {
        g_fb_width = fb_width;
        g_fb_height = fb_height;

        g_min_x = -fb_width  / 2;
        g_min_y = -fb_height / 2;
        g_max_x =  fb_width  / 2 - 1;
        g_max_y =  fb_height / 2 - 1;

        g_viewport[0][0] = (fb_width - 1)/2;
        g_viewport[1][1] = fb_height/2 - 1;
        g_viewport[2][2] = 1.0;
        g_viewport[3][3] = 1.0;

        g_framebuffer = framebuffer;
        g_depthbuffer = arena_alloc(arena, fb_width*fb_height*sizeof(double));
}

//
// point - draws a point at (x, y) in the frame buffer
//
void point(int x, int y, Vec3 colour) {
        g_framebuffer[SCREEN(x, y)] = RGBf(colour.x, colour.y, colour.z);
}

//
// clear_screen - fills the framebuffer with an 8-bit repeating pattern
//
void clear_screen(char c) {
        memset(g_framebuffer, c, g_fb_width*g_fb_height*sizeof(uint32_t));
}

//
// line - draws a line in the frame buffer from (ax, ay) to (bx, by)
//
void line(Vec3 a, Vec3 b, Vec3 colour) {
        if (a.x == b.x && a.y == b.y) {
                point(a.x, a.y, colour);
                return;
        }

        #define SWAP(T, x, y) { \
                T _tmp = x;     \
                (x) = (y);      \
                (y) = _tmp;     \
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

        #undef SWAP

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
// sample texture - samples a texture at (x, y)
//
static inline Vec3 sample_texture(Texture *texture, double x, double y) {
        int u = x*(texture->width - 1);
        int v = (1.0 - y)*(texture->height - 1);

        size_t i = v*texture->width + u;

        double one_255th = 1.0/255;

        return (Vec3){
                .x = (double)texture->data[i * 3*sizeof(uint8_t)]     * one_255th,
                        .y = (double)texture->data[i * 3*sizeof(uint8_t) + 1] * one_255th,
                        .z = (double)texture->data[i * 3*sizeof(uint8_t) + 2] * one_255th
        };
}

//
// attenuate - attenuates a light value based on a given distance
//
static inline double attenuate(double intensity, double dist) {
        intensity = intensity < 0 ? 0 : intensity/(dist*dist);
        return intensity > 1 ? 1 : intensity;
}

//
// render_face - renders a single face in a model to the framebuffer, performs texture mapping & Gouraud shading
//
static void render_face(Face *f, Model *model, Camera *cam, Light *light) {
        Vec3 *v = model->obj->vertices;
        Vec3 *vt = model->obj->uvs;
        Vec3 *vn = model->obj->norms;

        Vec3 a = v[f->v0[VERTEX]];
        Vec3 b = v[f->v1[VERTEX]];
        Vec3 c = v[f->v2[VERTEX]];

        Vec3 a_lightview = m4v3_mul(light->view_mat, a);
        Vec3 b_lightview = m4v3_mul(light->view_mat, b);
        Vec3 c_lightview = m4v3_mul(light->view_mat, c);

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

        Vec3 light_in_view = m4v3_mul(cam->view_mat, light->pos);

        Vec3 a_to_light = vec3_sub(light_in_view, a);
        Vec3 b_to_light = vec3_sub(light_in_view, b);
        Vec3 c_to_light = vec3_sub(light_in_view, c);

        Vec3 a_light_proj = vec3_scale(a_lightview, a_recip_z);
        Vec3 b_light_proj = vec3_scale(b_lightview, b_recip_z);
        Vec3 c_light_proj = vec3_scale(c_lightview, c_recip_z);

        double light_a = light->intensity * attenuate(vec3_dot(na, unit(a_to_light)), vec3_norm(a_to_light));
        double light_b = light->intensity * attenuate(vec3_dot(nb, unit(b_to_light)), vec3_norm(b_to_light));
        double light_c = light->intensity * attenuate(vec3_dot(nc, unit(c_to_light)), vec3_norm(c_to_light));

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

                                double intensity = 0.1;
                                if (smap_read(light->shadow_map, interp_light_proj.x, interp_light_proj.y) <= interp_light_proj.z + 0.05) {
                                        intensity = INTERP(light_a, light_b, light_c);
                                }

                                Vec3 colour = sample_texture(model->texture, uv.x, uv.y);
                                colour.x *= light->colour.x;
                                colour.y *= light->colour.y;
                                colour.z *= light->colour.z;
                                colour = vec3_scale(colour, intensity);

                                point(x, y, colour);
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
// render_model - renders all the faces in a model
//
void render_model(Model *model, Camera *cam, Light *light) {
        Obj *obj = model->obj;
        Face *f = obj->faces;

        for (int i = 0; i < obj->face_count; ++i) {
                render_face(&f[i], model, cam, light);
        }
}

/*
 * @XXX text rendering
 */
enum {
        GLYPH_END  = -1,
        GLYPH_LIFT = -2,

        GLYPH_UNKNOWN = 128,

        GLYPH_0 = 0,
        GLYPH_1,
        GLYPH_2,
        GLYPH_3,
        GLYPH_4,
        GLYPH_5,
        GLYPH_6,
        GLYPH_7,
        GLYPH_8,
        GLYPH_9,

        GLYPH_ZERO  = '0',
        GLYPH_ONE   = '1',
        GLYPH_TWO   = '2',
        GLYPH_THREE = '3',
        GLYPH_FOUR  = '4',
        GLYPH_FIVE  = '5',
        GLYPH_SIX   = '6',
        GLYPH_SEVEN = '7',
        GLYPH_EIGHT = '8',
        GLYPH_NINE  = '9',

        GLYPH_A = 'A',
        GLYPH_B = 'B',
        GLYPH_C = 'C',
        GLYPH_D = 'D',
        GLYPH_E = 'E',
        GLYPH_F = 'F',
        GLYPH_G = 'G',
        GLYPH_H = 'H',
        GLYPH_I = 'I',
        GLYPH_J = 'J',
        GLYPH_K = 'K',
        GLYPH_L = 'L',
        GLYPH_M = 'M',
        GLYPH_N = 'N',
        GLYPH_O = 'O',
        GLYPH_P = 'P',
        GLYPH_Q = 'Q',
        GLYPH_R = 'R',
        GLYPH_S = 'S',
        GLYPH_T = 'T',
        GLYPH_U = 'U',
        GLYPH_V = 'V',
        GLYPH_W = 'W',
        GLYPH_X = 'X',
        GLYPH_Y = 'Y',
        GLYPH_Z = 'Z',

        GLYPH_a = 'a',
        GLYPH_b = 'b',
        GLYPH_c = 'c',
        GLYPH_d = 'd',
        GLYPH_e = 'e',
        GLYPH_f = 'f',
        GLYPH_g = 'g',
        GLYPH_h = 'h',
        GLYPH_i = 'i',
        GLYPH_j = 'j',
        GLYPH_k = 'k',
        GLYPH_l = 'l',
        GLYPH_m = 'm',
        GLYPH_n = 'n',
        GLYPH_o = 'o',
        GLYPH_p = 'p',
        GLYPH_q = 'q',
        GLYPH_r = 'r',
        GLYPH_s = 's',
        GLYPH_t = 't',
        GLYPH_u = 'u',
        GLYPH_v = 'v',
        GLYPH_w = 'w',
        GLYPH_x = 'x',
        GLYPH_y = 'y',
        GLYPH_z = 'z',

        GLYPH_DOT         = '.',
        GLYPH_SPACE       = ' ',
        GLYPH_COMMA       = ',',
        GLYPH_LBRACK      = '[',
        GLYPH_RBRACK      = ']',
        GLYPH_LPAREN      = '(',
        GLYPH_RPAREN      = ')',
        GLYPH_LBRACE      = '{',
        GLYPH_RBRACE      = '}',
        GLYPH_QUESTION    = '?',
        GLYPH_UNDERSCORE  = '_',
        GLYPH_EXCLAMATION = '!',
        GLYPH_SQUOTE      = '\'',
        GLYPH_DQUOTE      = '"',
        GLYPH_DOLLAR      = '$',
        GLYPH_PERCENT     = '%',
        GLYPH_CARET       = '^',
        GLYPH_AMPERSAND   = '&',
        GLYPH_ASTERISK    = '*',
        GLYPH_MINUS       = '-',
        GLYPH_PLUS        = '+',
        GLYPH_AT          = '@',
        GLYPH_TILDE       = '~',
        GLYPH_HASH        = '#',
        GLYPH_EQUALS      = '=',
        GLYPH_LT          = '<',
        GLYPH_GT          = '>',
        GLYPH_COLON       = ':',
        GLYPH_SEMICOLON   = ';',
        GLYPH_FSLASH      = '/',
        GLYPH_BSLASH      = '\\',
        GLYPH_PIPE        = '|',
        GLYPH_BACKTICK    = '`',
};

static int16_t g_glyph_seqs[][16] = {
        [GLYPH_UNKNOWN] = {0, 2, 11, 9, 0, 11, 9, 2, GLYPH_END},

        [GLYPH_0] = {9, 3, 5, 11, 9, 5,          GLYPH_END},
        [GLYPH_1] = {6, 4, 10, 9, 11,            GLYPH_END},
        [GLYPH_2] = {3, 5, 8, 6, 9, 11,          GLYPH_END},
        [GLYPH_3] = {3, 5, 8, 6, 8, 11, 9,       GLYPH_END},
        [GLYPH_4] = {11, 5, 6, 8,                GLYPH_END},
        [GLYPH_5] = {5, 3, 6, 8, 11, 9,          GLYPH_END},
        [GLYPH_6] = {5, 6, 8, 11, 9, 6,          GLYPH_END},
        [GLYPH_7] = {3, 5, 9,                    GLYPH_END},
        [GLYPH_8] = {8, 6, 3, 5, 11, 9, 6,       GLYPH_END},
        [GLYPH_9] = {5, 6, 3, 5, 8, 6, 8, 11, 9, GLYPH_END},

        [GLYPH_ZERO]  = {9, 3, 5, 11, 9, 5,          GLYPH_END},
        [GLYPH_ONE]   = {6, 4, 10, 9, 11,            GLYPH_END},
        [GLYPH_TWO]   = {3, 5, 8, 6, 9, 11,          GLYPH_END},
        [GLYPH_THREE] = {3, 5, 8, 6, 8, 11, 9,       GLYPH_END},
        [GLYPH_FOUR]  = {11, 5, 6, 8,                GLYPH_END},
        [GLYPH_FIVE]  = {5, 3, 6, 8, 11, 9,          GLYPH_END},
        [GLYPH_SIX]   = {5, 6, 8, 11, 9, 6,          GLYPH_END},
        [GLYPH_SEVEN] = {3, 5, 9,                    GLYPH_END},
        [GLYPH_EIGHT] = {8, 6, 3, 5, 11, 9, 6,       GLYPH_END},
        [GLYPH_NINE]  = {5, 6, 3, 5, 8, 6, 8, 11, 9, GLYPH_END},

        [GLYPH_A] = {6, 8, GLYPH_LIFT, 9, 3, 1, 5, 11,           GLYPH_END},
        [GLYPH_B] = {3, 9, 10, 8, 4, 3, 0, 1, 5, 4,              GLYPH_END},
        [GLYPH_C] = {11, 10, 6, 3, 1, 2,                         GLYPH_END},
        [GLYPH_D] = {0, 1, 5, 8, 10, 9, 0,                       GLYPH_END},
        [GLYPH_E] = {2, 0, 9, 11, GLYPH_LIFT, 3, 5,              GLYPH_END},
        [GLYPH_F] = {2, 0, 3, 5, 3, 9,                           GLYPH_END},
        [GLYPH_G] = {2, 0, 9, 11, 8, 7,                          GLYPH_END},
        [GLYPH_H] = {0, 9, GLYPH_LIFT, 2, 11, GLYPH_LIFT, 3, 5,  GLYPH_END},
        [GLYPH_I] = {0, 2, GLYPH_LIFT, 9, 11, GLYPH_LIFT, 1, 10, GLYPH_END},
        [GLYPH_J] = {0, 2, GLYPH_LIFT, 1, 10, 9,                 GLYPH_END},
        [GLYPH_K] = {0, 9, GLYPH_LIFT, 2, 3, 11,                 GLYPH_END},
        [GLYPH_L] = {0, 9, 11,                                   GLYPH_END},
        [GLYPH_M] = {9, 0, 7, 2, 11,                             GLYPH_END},
        [GLYPH_N] = {9, 0, 11, 2,                                GLYPH_END},
        [GLYPH_O] = {9, 0, 2, 11, 9,                             GLYPH_END},
        [GLYPH_P] = {9, 0, 2, 8, 6,                              GLYPH_END},
        [GLYPH_Q] = {3, 1, 5, 8, 10, 6, 3, GLYPH_LIFT, 7, 11,    GLYPH_END},
        [GLYPH_R] = {9, 0, 1, 5, 7, 11, 7, 6,                    GLYPH_END},
        [GLYPH_S] = {2, 3, 8, 9,                                 GLYPH_END},
        [GLYPH_T] = {1, 10, GLYPH_LIFT, 0, 2,                    GLYPH_END},
        [GLYPH_U] = {0, 9, 11, 2,                                GLYPH_END},
        [GLYPH_V] = {0, 10, 2,                                   GLYPH_END},
        [GLYPH_W] = {0, 9, 7, 11, 2,                             GLYPH_END},
        [GLYPH_X] = {9, 2, 7, 0, 11,                             GLYPH_END},
        [GLYPH_Y] = {0, 7, 2, 7, 10,                             GLYPH_END},
        [GLYPH_Z] = {0, 2, 9, 11,                                GLYPH_END},

        [GLYPH_a] = {3, 5, 11, 9, 6, 8,             GLYPH_END},
        [GLYPH_b] = {3, 9, 11, 8, 6,                GLYPH_END},
        [GLYPH_c] = {11, 9, 3, 5,                   GLYPH_END},
        [GLYPH_d] = {5, 11, 9, 6, 8,                GLYPH_END},
        [GLYPH_e] = {11, 9, 3, 5, 6,                GLYPH_END},
        [GLYPH_f] = {5, 4, 10, 9, GLYPH_LIFT, 6, 8, GLYPH_END},
        [GLYPH_g] = {9, 11,5, 3, 6, 8,              GLYPH_END},
        [GLYPH_h] = {3, 9, 6, 8, 11,                GLYPH_END},
        [GLYPH_i] = {4, 10,                         GLYPH_END},
        [GLYPH_j] = {4, 10, 9,                      GLYPH_END},
        [GLYPH_k] = {3, 9, 6, 11, 6, 5,             GLYPH_END},
        [GLYPH_l] = {4, 10, 11,                     GLYPH_END},
        [GLYPH_m] = {9, 3, 7, 5, 11,                GLYPH_END},
        [GLYPH_n] = {9, 3, 11, 5,                   GLYPH_END},
        [GLYPH_o] = {9, 3, 5, 11, 9,                GLYPH_END},
        [GLYPH_p] = {9, 3, 5, 8, 6,                 GLYPH_END},
        [GLYPH_q] = {8, 10, 4, 3, 6, 7,             GLYPH_END},
        [GLYPH_r] = {9, 6, 5,                       GLYPH_END},
        [GLYPH_s] = {5, 6, 8, 9,                    GLYPH_END},
        [GLYPH_t] = {7, 6, 3, 9, 11,                GLYPH_END},
        [GLYPH_u] = {3, 9, 11, 5,                   GLYPH_END},
        [GLYPH_v] = {3, 10, 5,                      GLYPH_END},
        [GLYPH_w] = {3, 9, 7, 11, 5,                GLYPH_END},
        [GLYPH_x] = {9, 5, 7, 3, 11,                GLYPH_END},
        [GLYPH_y] = {3, 7, 5, 7, 10,                GLYPH_END},
        [GLYPH_z] = {3, 5, 9, 11,                   GLYPH_END},

        [GLYPH_PERCENT]     = {0, 1, 4, 3, 0, GLYPH_LIFT, 9, 2, GLYPH_LIFT, 11, 10, 7, 8, 11, GLYPH_END},
        [GLYPH_ASTERISK]    = {3, 11, GLYPH_LIFT, 9, 5, GLYPH_LIFT, 6, 8, GLYPH_LIFT, 4, 10,  GLYPH_END},
        [GLYPH_HASH]        = {3, 5, GLYPH_LIFT, 6, 8, GLYPH_LIFT, 9, 1, GLYPH_LIFT, 10, 2,   GLYPH_END},
        [GLYPH_QUESTION]    = {10, 10, GLYPH_LIFT, 7, 4, 5, 2, 0,                             GLYPH_END},
        [GLYPH_DOLLAR]      = {9, 8, 3, 2, GLYPH_LIFT, 1, 10,                                 GLYPH_END},
        [GLYPH_EXCLAMATION] = {10, 10, GLYPH_LIFT, 7, 1,                                      GLYPH_END},
        [GLYPH_COLON]       = {4, 4, GLYPH_LIFT, 10, 10,                                      GLYPH_END},
        [GLYPH_PLUS]        = {6, 8, GLYPH_LIFT, 4, 10,                                       GLYPH_END},
        [GLYPH_AT]          = {11, 9, 0, 2, 8, 7, 4, 5,                                       GLYPH_END},
        [GLYPH_DQUOTE]      = {0, 3, GLYPH_LIFT, 1, 4,                                        GLYPH_END},
        [GLYPH_EQUALS]      = {3, 5, GLYPH_LIFT, 6, 8,                                        GLYPH_END},
        [GLYPH_SEMICOLON]   = {4, 4, GLYPH_LIFT, 7, 9,                                        GLYPH_END},
        [GLYPH_AMPERSAND]   = {11, 3, 1, 5, 6, 9, 8,                                          GLYPH_END},
        [GLYPH_RBRACK]      = {1, 2, 11, 10,                                                  GLYPH_END},
        [GLYPH_LPAREN]      = {2, 4, 7, 11,                                                   GLYPH_END},
        [GLYPH_LBRACE]      = {2, 4, 7, 11,                                                   GLYPH_END},
        [GLYPH_LBRACK]      = {1, 0, 9, 10,                                                   GLYPH_END},
        [GLYPH_RPAREN]      = {0, 4, 7, 9,                                                    GLYPH_END},
        [GLYPH_TILDE]       = {6, 4, 8, 5,                                                    GLYPH_END},
        [GLYPH_RBRACE]      = {0, 4, 7, 9,                                                    GLYPH_END},
        [GLYPH_LT]          = {5, 6, 11,                                                      GLYPH_END},
        [GLYPH_GT]          = {3, 8, 9,                                                       GLYPH_END},
        [GLYPH_CARET]       = {3, 1, 5,                                                       GLYPH_END},
        [GLYPH_UNDERSCORE]  = {9, 11,                                                         GLYPH_END},
        [GLYPH_BSLASH]      = {0, 10,                                                         GLYPH_END},
        [GLYPH_PIPE]        = {1, 10,                                                         GLYPH_END},
        [GLYPH_MINUS]       = {6, 8,                                                          GLYPH_END},
        [GLYPH_COMMA]       = {7, 9,                                                          GLYPH_END},
        [GLYPH_SQUOTE]      = {1, 4,                                                          GLYPH_END},
        [GLYPH_DOT]         = {9, 9,                                                          GLYPH_END},
        [GLYPH_FSLASH]      = {9, 1,                                                          GLYPH_END},
        [GLYPH_BACKTICK]    = {0, 4,                                                          GLYPH_END},
        [GLYPH_SPACE]       = {0,                                                             GLYPH_END},
};

//
// is_valid_glyph - checks whether a glyph is valid
//
static inline int is_valid_glyph(uint8_t glyph_id) {
        return glyph_id <= 9 || glyph_id == GLYPH_SPACE || isprint(glyph_id);
}

//
// render_glyph - renders a glyph to the framebuffer
//
static Vec3 render_glyph(Vec3 pos, double unit, Vec3 colour, uint8_t glyph_id) {
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

        if (!is_valid_glyph(glyph_id)) {
                glyph_id = GLYPH_UNKNOWN;
        }

        Vec3 line_start = points[g_glyph_seqs[glyph_id][0]];
        for (int i = 1; g_glyph_seqs[glyph_id][i] != GLYPH_END; ++i) {
                int point_index = g_glyph_seqs[glyph_id][i];
                if (point_index == GLYPH_LIFT) {
                        line_start = points[g_glyph_seqs[glyph_id][i++ + 1]];
                        continue;
                }
                line(line_start, points[point_index], colour);
                line_start = points[point_index];
        }

        pos.x += ((glyph_id == GLYPH_DOT) ? 1 : 3) * unit;
        return pos;
}

//
// render_text_int - called by render_text to render integers as text to the framebuffer
//
static Vec3 render_text_int(Vec3 pos, double unit, Vec3 colour, int d) {
        if (d == 0) {
                return render_glyph(pos, unit, colour, GLYPH_0);
        }

        if (d < 0) {
                pos = render_glyph(pos, unit, colour, GLYPH_MINUS);
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
                pos = render_glyph(pos, unit, colour, (d / place_val) % 10);
                place_val /= 10;
        }

        return pos;
}

//
// render_text_double - called by render_text to render doubles as text to the framebufer
//
static Vec3 render_text_double(Vec3 pos, double unit, Vec3 colour, double f) {
        if (f < 0) {
                pos = render_glyph(pos, unit, colour, GLYPH_MINUS);
                f = -f;
        }

        long whole = (long)f;

        if (whole == 0) {
                pos = render_glyph(pos, unit, colour, GLYPH_0);
        }

        int place_val = 1;
        while (place_val < whole) {
                place_val *= 10;
        }

        if (place_val > whole) {
                place_val /= 10;
        }

        while (place_val >= 1) {
                pos = render_glyph(pos, unit, colour, (whole / place_val) % 10);
                place_val /= 10;
        }

        pos = render_glyph(pos, unit, colour, GLYPH_DOT);

        int max_digits = 2;
        do {
                f *= 10;
                pos = render_glyph(pos, unit, colour, (long)f % 10);
        } while (f > (long)f && max_digits-- > 1);

        return pos;
}

//
// render_text_vector - called by render_text to render a vector as text to the framebuffer
//
static Vec3 render_text_vector(Vec3 pos, double unit, Vec3 colour, const Vec3 *v) {
        pos = render_glyph(pos, unit, colour, GLYPH_LPAREN);
        pos = render_text_double(pos, unit, colour, v->x);
        pos = render_glyph(pos, unit, colour, GLYPH_COMMA);
        pos = render_text_double(pos, unit, colour, v->y);
        pos = render_glyph(pos, unit, colour, GLYPH_COMMA);
        pos = render_text_double(pos, unit, colour, v->z);
        pos = render_glyph(pos, unit, colour, GLYPH_RPAREN);
        return pos;
}

//
// render_text - renders formatted text to the framebuffer
//
Vec3 render_text(Vec3 pos, double unit, Vec3 colour, const char *fmt, ...) {
        va_list arglist;
        va_start(arglist, fmt);

        if (pos.y >= 1.0 || pos.x <= -1.0)  {
                goto _exit;
        }

        Vec3 curs = pos;
        for (int i = 0; fmt[i] != '\0'; ++i) {
                if (curs.y - 4*unit <= -1.0) {
                        goto _exit;
                }

                if (fmt[i] == '\n' || curs.x + 3*unit >= 1.0) {
                        curs.x = pos.x;
                        curs.y -= 4*unit;
                } else if (fmt[i] == '%') {
                        ++i;
                        if (fmt[i] == 'd') {
                                curs = render_text_int(curs, unit, colour, va_arg(arglist, int));
                        } else if (fmt[i] == 'f') {
                                curs = render_text_double(curs, unit, colour, va_arg(arglist, double));
                        } else if (fmt[i] == 'v') {
                                curs = render_text_vector(curs, unit, colour, va_arg(arglist, Vec3*));
                        } else if (fmt[i] == 's') {
                                curs = render_text(curs, unit, colour, va_arg(arglist, char*));
                        } else if (fmt[i] == 'c') {
                                curs = render_glyph(curs, unit, colour, va_arg(arglist, int));
                        } else if (fmt[i] == '%') {
                                curs = render_glyph(curs, unit, colour, GLYPH_PERCENT);
                        }
                } else {
                        curs = render_glyph(curs, unit, colour, fmt[i]);
                }
        }

_exit:
        va_end(arglist);
        return curs;
}

/*
 * @XXX asset loading @XXX
 */

//
// parse_double - takes a pointer to a string storing a textual floating point number and returns that number as a double
//
static double parse_double(char **src) {
        char *p = *src;
        double result = 0.0;

        int sign;
        if (*p == '-') {
                ++p;
                sign = -1;
        } else {
                sign = 1;
        }

        while (*p >= '0' && *p <= '9') {
                result *= 10;
                result += *p - '0';
                ++p;
        }

        double sig = 0.1;
        if (*p++ == '.') {
                while (*p >= '0' && *p <= '9') {
                        result += (double)(*p - '0') * sig;
                        sig *= 0.1;
                        ++p;
                }
        }

        *src = p;
        return result * sign;
}

//
// parse_int - takes a pointer to a string storing a textual integer and returns that integer as an int
//
static int parse_int(char **src) {
        char *p = *src;
        int result = 0;

        while (*p >= '0' && *p <= '9') {
                result *= 10;
                result += *p - '0';
                ++p;
        }

        *src = p;
        return result;
}

//
// is_space - returns true if a character is a whitespace character, else false
//
static int is_space(unsigned c) {
        return (c == ' ')      ||
               (c - '\a' <= 6);
}

//
// load_obj - loads an OBJ file from disk and parses geometry data into vertex arrays etc.
//
static Obj *load_obj(Arena *arena, char *filename) {
        FILE *obj_file = fopen(filename, "rb");
        if (!obj_file) {
                return &g_error_obj;
        }

        fseek(obj_file, 0, SEEK_END);
        size_t obj_file_len = ftell(obj_file);
        rewind(obj_file);
        char *src = arena_alloc(arena, obj_file_len + 1);
        fread(src, 1, obj_file_len, obj_file);
        fclose(obj_file);

        Obj *obj = arena_alloc(arena, sizeof(Obj));

        size_t vert_buf_size = 1 << 12;
        obj->vertices = arena_alloc(arena, vert_buf_size);

        #define SKIP_SPACE    {while (is_space(*src)) ++src;}
        #define SKIP_TO_SPACE {while (!is_space(*src)) ++src;}

        int i = 0;
        for (; *src == 'v' && src[1] == ' '; ++i) {
                src += 2;
                if (i * sizeof(Vec3) >= vert_buf_size) {
                        arena_commit_at(arena, obj->vertices, vert_buf_size *= 2);
                }

                SKIP_SPACE;
                obj->vertices[i].x = parse_double(&src);
                SKIP_SPACE;
                obj->vertices[i].y = parse_double(&src);
                SKIP_SPACE;
                obj->vertices[i].z = parse_double(&src);
                SKIP_SPACE;
        }
        obj->vertex_count = i;

        while (!(*src == 'v' && src[1] == 't')) {
                ++src;
        }

        size_t uv_buf_size = 1 << 12;
        obj->uvs = arena_alloc(arena, uv_buf_size);

        i = 0;
        for (; *src == 'v' && src[1] == 't'; ++i) {
                src += 2;
                if (i * sizeof(Vec3) >= uv_buf_size) {
                        arena_commit_at(arena, obj->uvs, uv_buf_size *= 2);
                }

                SKIP_SPACE;
                obj->uvs[i].x = parse_double(&src);
                SKIP_SPACE;
                obj->uvs[i].y = parse_double(&src);
                SKIP_SPACE;
                obj->uvs[i].z = parse_double(&src);
                SKIP_SPACE;
        }
        obj->uv_count = i;

        while (!(*src == 'v' && src[1] == 'n')) {
                ++src;
        }

        size_t norm_buf_size = 1 << 12;
        obj->norms = arena_alloc(arena, norm_buf_size);

        i = 0;
        for (; *src == 'v' && src[1] == 'n'; ++i) {
                src += 2;
                if (i * sizeof(Vec3) >= norm_buf_size) {
                        arena_commit_at(arena, obj->norms, norm_buf_size *= 2);
                }

                SKIP_SPACE;
                obj->norms[i].x = parse_double(&src);
                SKIP_SPACE;
                obj->norms[i].y = parse_double(&src);
                SKIP_SPACE;
                obj->norms[i].z = parse_double(&src);
                SKIP_SPACE;
        }
        obj->norm_count = i;

        while (*src != 'f') {
                ++src;
        }

        size_t face_buf_size = 1 << 12;
        obj->faces = arena_alloc(arena, face_buf_size);

        i = 0;
        for (; *src == 'f'; ++i) {
                ++src;
                if (i * sizeof(Face) >= face_buf_size) {
                        arena_commit_at(arena, obj->faces, face_buf_size *= 2);
                }
                
                SKIP_SPACE;
                obj->faces[i].v0[VERTEX] = parse_int(&src) - 1;
                ++src;
                obj->faces[i].v0[UV] = parse_int(&src) - 1;
                ++src;
                obj->faces[i].v0[NORM] = parse_int(&src) - 1;

                SKIP_SPACE;
                obj->faces[i].v1[VERTEX] = parse_int(&src) - 1;
                ++src;
                obj->faces[i].v1[UV] = parse_int(&src) - 1;
                ++src;
                obj->faces[i].v1[NORM] = parse_int(&src) - 1;


                SKIP_SPACE;
                obj->faces[i].v2[VERTEX] = parse_int(&src) - 1;
                ++src;
                obj->faces[i].v2[UV] = parse_int(&src) - 1;
                ++src;
                obj->faces[i].v2[NORM] = parse_int(&src) - 1;

                SKIP_SPACE;
        }
        obj->face_count = i;

        return obj;
}

//
// load_texture - loads a texture from disk
//
static Texture *load_texture(Arena *arena, char *filename, size_t w, size_t h) {
        Texture *texture = arena_alloc(arena, sizeof(Texture));
        texture->width = w;
        texture->height = h;

        FILE *tex_file = fopen(filename, "rb");
        if (!tex_file) {
                return &g_error_tex;
        }

        fseek(tex_file, 0, SEEK_END);
        size_t tex_file_len = ftell(tex_file);
        rewind(tex_file);

        texture->data = arena_alloc(arena, tex_file_len);
        fread(texture->data, 1, tex_file_len, tex_file);
        fclose(tex_file);

        return texture;
}

//
// load_model - loads a model by loading specified OBJ, texture, and normal map for the model, falls back to error assets on failure
//
Model *load_model(
        Arena *arena,
        char *obj_filename,
        char *tm_filename,
        char *nm_filename,
        size_t tm_w,
        size_t tm_h,
        size_t nm_w,
        size_t nm_h
) {
        Model *model = arena_alloc(arena, sizeof(Model));

        model->obj = load_obj(arena, obj_filename);
        model->texture = load_texture(arena, tm_filename, tm_w, tm_h);
        model->norm_map = load_texture(arena, nm_filename, nm_w, nm_h);

        return model;
}

/*
 * @XXX misc. @XXX
 */

Scene *mkscene(Arena *arena) {
        // @TODO might want to add the ability to add multiple reserved regions to an arena (linked list) so each
        // growable array (models, lights, etc.) can be allocated at the same time and just grow without interfering with
        // each other, I think this would be best implemented with a separate alloc function like arena_alloc_group
        // and then just have allocate_at (I feel like this also needs to be smoothed out...) look at the pointer and
        // address the appropriate region, rather than just bumping the end pointer
        return arena_alloc(arena, sizeof(Scene));
}

void add_model(Model *model) {
        
}

void add_light(Light *light) {

}
