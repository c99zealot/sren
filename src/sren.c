// @TODO do calculations in fixed-point when rasterising, precision errors are causing smearing in large faces.
//       rendering models in order of descending distance from the camera seems to fix this for now.

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include <stdbool.h>
#include <float.h>

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

Mat4 g_viewport; // @XXX exposing this for now...

Arena g_arena;

/*
 * @XXX shadow mapping @XXX
 */

//
// reset_smap - resets a Light's shadow map to the maximum depth value
//
static void reset_smap(Shadow_Map *smap) {
        for (int i = 0; i < smap->width*smap->height; ++i) {
                smap->depths[i] = -1024.0;
        }
}

//
// mk_smap - allocates and initialises a Shadow_Map
//
static Shadow_Map *mk_smap(size_t w, size_t h, int lifetime) {
        Shadow_Map *smap = arena_alloc(&g_arena, sizeof(Shadow_Map));
        *smap = (Shadow_Map){
                .width = w,
                .height = h,
                .min_x = -w/2,
                .min_y = -h/2,
                .max_x = w/2 - 1,
                .max_y = h/2 - 1,
                .depths = arena_alloc(&g_arena, w*h*sizeof(double)),
                .lifetime = lifetime,
                .age = 0,
        };

        reset_smap(smap);
        return smap;
}

/*
 * @XXX shadow mapping @XXX
 */

//
// reset_fmap - resets a Light's filter map to transparent white
//
static void reset_fmap(Filter_Map *fmap) {
        memset(fmap->colours, 0xFF, fmap->width*fmap->height*sizeof(uint32_t));
}

//
// mk_fmap - allocates and initialises a Shadow_Map
//
static Filter_Map *mk_fmap(size_t w, size_t h, int lifetime) {
        Filter_Map *fmap = arena_alloc(&g_arena, sizeof(Filter_Map));
        *fmap = (Filter_Map){
                .width = w,
                .height = h,
                .min_x = -w/2,
                .min_y = -h/2,
                .max_x = w/2 - 1,
                .max_y = h/2 - 1,
                .colours = arena_alloc(&g_arena, w*h*sizeof(uint32_t)),
                .lifetime = lifetime,
                .age = 0,
        };

        reset_fmap(fmap);
        return fmap;
}

//
// smap_write - writes a depth value at (x, y) in a shadow map
//
static inline void smap_write(Shadow_Map *smap, int x, int y, double depth) {
        size_t i = (smap->width * (smap->height/2 - y) - smap->width/2 + x);
        smap->depths[i] = depth;
}

//
// smap_read - reads the depth value at (x, y) in a shadow map
//
double smap_read(Shadow_Map *smap, int x, int y) {
        size_t i = (smap->width * (smap->height/2 - y) - smap->width/2 + x);
        return (i >= smap->width*smap->height) ? -1024.0 : smap->depths[i];
}

//
// fmap_write - writes a colour value at (x, y) in a filter map
//
static inline void fmap_write(Filter_Map *fmap, int x, int y, Vec4 colour) {
        size_t i = (fmap->width * (fmap->height/2 - y) - fmap->width/2 + x);
        fmap->colours[i] = vec4_to_rgba(colour);
}

//
// fmap_read - reads the colour value at (x, y) in a shadow map
//
Vec4 fmap_read(Filter_Map *fmap, int x, int y) {
        size_t i = (fmap->width * (fmap->height/2 - y) - fmap->width/2 + x);
        return (i >= fmap->width*fmap->height) ? (Vec4){0} : rgba_to_vec4(fmap->colours[i]);
}

//
// signed_tri_area2 - calculates twice the signed area of the triangle abc
//
static inline double signed_tri_area2(Vec3 a, Vec3 b, Vec3 c) {
        return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
}

// @XXX render_face_smap needs this, remove when things are better organised
static inline Vec4 sample_texture(Texture *texture, double x, double y);

//
// render_face_smap - renders a single face in a model to a shadow map
//
static void render_face_smap(Face *f, Model *model, Light *light) {
        Vec3 *v = model->obj->vertices;
        Vec3 *vt = model->obj->uvs;

        Vec3 a = m4v3_mul(light->view_mat, v[f->v0[VERTEX]]);
        Vec3 b = m4v3_mul(light->view_mat, v[f->v1[VERTEX]]);
        Vec3 c = m4v3_mul(light->view_mat, v[f->v2[VERTEX]]);

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

        a = m4v3_mul(light->viewport, a);
        b = m4v3_mul(light->viewport, b);
        c = m4v3_mul(light->viewport, c);

        double recip_area = 1.0/signed_tri_area2(a, b, c);

        int min_x = MIN3(a.x, b.x, c.x);
        int max_x = MAX3(a.x, b.x, c.x);
        int min_y = MIN3(a.y, b.y, c.y);
        int max_y = MAX3(a.y, b.y, c.y);

        min_x = MAX(min_x, light->shadow_map->min_x);
        min_y = MAX(min_y, light->shadow_map->min_y);
        max_x = MIN(max_x, light->shadow_map->max_x);
        max_y = MIN(max_y, light->shadow_map->max_y);

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

        #define INTERP(r, s, t)     (recip_area * (Ea*(r) + Eb*(s) + Ec*(t)))
        #define VEC_INTERP(r, s, t) (vec3_scale(vec3_add(vec3_scale(r, Ea), vec3_add(vec3_scale(s, Eb), \
                                        vec3_scale(t, Ec))), recip_area))

        for (int y = max_y; y >= min_y; --y) {
                for (int x = min_x; x <= max_x; ++x) {
                        if (Ea >= 0 && Eb >= 0 && Ec >= 0) {
                                double interp_recip_z = INTERP(a_recip_z, b_recip_z, c_recip_z);
                                Vec3 uv = vec3_scale(VEC_INTERP(ta, tb, tc), 1.0/interp_recip_z);
                                Vec4 texel = sample_texture(model->texture, uv.x, uv.y);

                                double depth = INTERP(a.z, b.z, c.z);
                                if (texel.w == 1 && depth > smap_read(light->shadow_map, x, y)) {
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

        #undef INTERP
        #undef VEC_INTERP
}

//
// render_face_fmap - renders a single face in a model with transparency to a filter map
//
static void render_face_fmap(Face *f, Model *model, Light *light) {
        Vec3 *v = model->obj->vertices;
        Vec3 *vt = model->obj->uvs;

        Vec3 a = m4v3_mul(light->view_mat, v[f->v0[VERTEX]]);
        Vec3 b = m4v3_mul(light->view_mat, v[f->v1[VERTEX]]);
        Vec3 c = m4v3_mul(light->view_mat, v[f->v2[VERTEX]]);

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

        a = m4v3_mul(light->viewport, a);
        b = m4v3_mul(light->viewport, b);
        c = m4v3_mul(light->viewport, c);

        double recip_area = 1.0/signed_tri_area2(a, b, c);

        int min_x = MIN3(a.x, b.x, c.x);
        int max_x = MAX3(a.x, b.x, c.x);
        int min_y = MIN3(a.y, b.y, c.y);
        int max_y = MAX3(a.y, b.y, c.y);

        min_x = MAX(min_x, light->filter_map->min_x);
        min_y = MAX(min_y, light->filter_map->min_y);
        max_x = MIN(max_x, light->filter_map->max_x);
        max_y = MIN(max_y, light->filter_map->max_y);

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

        #define INTERP(r, s, t)     (recip_area * (Ea*(r) + Eb*(s) + Ec*(t)))
        #define VEC_INTERP(r, s, t) (vec3_scale(vec3_add(vec3_scale(r, Ea), vec3_add(vec3_scale(s, Eb), \
                                        vec3_scale(t, Ec))), recip_area))

        for (int y = max_y; y >= min_y; --y) {
                for (int x = min_x; x <= max_x; ++x) {
                        if (Ea >= 0 && Eb >= 0 && Ec >= 0) {
                                double interp_recip_z = INTERP(a_recip_z, b_recip_z, c_recip_z);
                                Vec3 uv = vec3_scale(VEC_INTERP(ta, tb, tc), 1.0/interp_recip_z);
                                Vec4 texel = sample_texture(model->texture, uv.x, uv.y);

                                double depth = INTERP(a.z, b.z, c.z);
                                if (depth > smap_read(light->shadow_map, x, y)) {
                                        fmap_write(light->filter_map, x, y, texel);
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

        #undef INTERP
        #undef VEC_INTERP
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
// render_model_fmap - renders a model to a filter map
//
void render_model_fmap(Model *model, Light *light) {
        Obj *obj = model->obj;
        Vec3 *v = obj->vertices;
        Face *f = obj->faces;

        for (int i = 0; i < obj->face_count; ++i) {
                render_face_fmap(&f[i], model, light);
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
static void reset_dbuf(void) {
        memset(g_depthbuffer, 0, g_fb_width*g_fb_height*sizeof(double));
}

/*
 * @XXX rendering @XXX
 */

//
// init_renderer - initialises internal data for rendering
//
void init_renderer(uint32_t *framebuffer, size_t fb_width, size_t fb_height) {
        arena_init(&g_arena, ARENA_RESERVE_DEFAULT, 1 << 24, 1);

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
        g_depthbuffer = arena_alloc(&g_arena, fb_width*fb_height*sizeof(double));
}

//
// deinit_renderer - deinitialises a renderer, releasing allocations etc.
//
void deinit_renderer(void) {
        arena_deinit(&g_arena);
}

//
// init_light - initialises a Light
//
void init_light(Light *light, size_t map_w, size_t map_h, int smap_lt, int fmap_lt) {
        light->shadow_map = mk_smap(map_w, map_h, smap_lt);
        light->filter_map = mk_fmap(map_w, map_h, fmap_lt);

        memset(light->viewport, 0, sizeof(Mat4));
        light->viewport[0][0] = (map_w - 1)/2;
        light->viewport[1][1] = map_h/2 - 1;
        light->viewport[2][2] = 1.0;
        light->viewport[3][3] = 1.0;
}

//
// point - draws a point at (x, y) in the frame buffer
//
void point(int x, int y, Vec4 colour) {
        g_framebuffer[SCREEN(x, y)] = RGBf(colour.x, colour.y, colour.z);
}

//
// get_pixel - reads an RGBA value from the frame buffer and returns it as a Vec4
//
Vec4 get_pixel(int x, int y) {
        return xrgb_to_vec4(g_framebuffer[SCREEN(x, y)]);
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
void line(Vec3 a, Vec3 b, Vec4 colour) {
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
static inline Vec4 sample_texture(Texture *texture, double x, double y) {
        size_t u = x*(texture->width - 1);
        size_t v = (1.0 - y)*(texture->height - 1);

        size_t i = 4*sizeof(uint8_t) * (v*texture->width + u);

        Vec4 result = {
                .x = texture->data[i],
                .y = texture->data[i + 1],
                .z = texture->data[i + 2],
                .w = texture->data[i + 3]
        };

        return vec4_scale(result, 1.0/255);
}

//
// out_of_view - checks whether a point resides outside of the screen
//
int out_of_view(Vec3 v) {
        return (v.x > g_max_x) ||
               (v.x < g_min_x) ||
               (v.y > g_max_y) ||
               (v.y < g_min_y);
}

//
// alpha_blend - alpha blend two colour vectors
//
static inline Vec4 alpha_blend(Vec4 a, Vec4 b) {
        return vec4_add(
                vec4_scale(a, a.w),
                vec4_scale(b, 1.0 - a.w)
        );
}

//
// render_image_fragment - renders a rectangular fragment of a Texture
//
void render_image_fragment(Vec3 pos, Vec4 colour, Texture *texture, Vec3 uv, double frag_w, double frag_h, double scale) {
        pos = m4v3_mul(g_viewport, pos);
        if (out_of_view(pos)) {
                return;
        }

        size_t target_w = texture->width * frag_w * scale;
        size_t target_h = texture->height * frag_h * scale;

        target_w = (pos.x + target_w > g_max_x) ? g_max_x - pos.x : target_w;
        target_h = (pos.y + target_h > g_max_y) ? g_max_y - pos.y : target_h;

        double dx = 1.0/target_w;
        double dy = 1.0/target_h;

        for (double y = 0; y <= 1.0; y += dy) {
                for (double x = 0; x <= 1.0; x += dx) {
                        size_t target_x = round(pos.x + x*target_w);
                        size_t target_y = round(pos.y + y*target_h);

                        Vec4 texel = vec4_mul(sample_texture(texture, uv.x + x*frag_w, uv.y + y*frag_h), colour);
                        Vec4 pixel = get_pixel(target_x, target_y);

                        point(target_x, target_y, alpha_blend(texel, pixel));
                }
        }
}

//
// render_glyph - renders a glyph from a fontset texture
//
Vec3 render_glyph(Vec3 pos, Vec4 colour, double scale, Texture *fontset, char c) {
        double glyph_width = 77.0/fontset->width;

        if (!isprint(c)) {
                render_image_fragment(pos, colour, fontset, VEC3(('~' - ' ' + 1) * glyph_width, 0, 0), glyph_width, 1.0, scale);
                pos.x += (77.0 * scale)/g_max_x;
                return pos;
        }

        Vec3 uv = VEC3((c - ' ') * glyph_width, 0, 0);
        render_image_fragment(pos, colour, fontset, uv, glyph_width, 1.0, scale);

        pos.x += (77.0 * scale)/g_max_x;
        return pos;
}

//
// vec4_clamp - clamps a 4D vector's components to (-inf, 1]
//
static inline Vec4 vec4_clamp(Vec4 v) {
        v.x = MIN(1, v.x);
        v.y = MIN(1, v.y);
        v.z = MIN(1, v.z);
        v.w = MIN(1, v.w);

        return v;
}

//
// reflect - reflects a unit vector r about another unit vector n
//
static inline Vec3 reflect(Vec3 r, Vec3 n) {
        return vec3_sub(vec3_scale(n, 2*vec3_dot(n, r)), r);
}

//
// render_face - renders a single, textured & lit face from a model
//
static void render_face(Face *f, Model *model, Camera *cam, Light *light) {
        Vec3 *v = model->obj->vertices;
        Vec3 *vt = model->obj->uvs;
        Vec3 *vn = model->obj->normals;

        Material *mat = model->material;

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

        Vec3 light_camview = m4v3_mul(cam->view_mat, light->pos);
        Vec3 a_to_light = vec3_sub(light_camview, a);
        Vec3 b_to_light = vec3_sub(light_camview, b);
        Vec3 c_to_light = vec3_sub(light_camview, c);

        double a_recip_z = 1.0/a.z;
        double b_recip_z = 1.0/b.z;
        double c_recip_z = 1.0/c.z;

        Vec3 a_to_cam = vec3_scale(a, -a_recip_z);
        Vec3 b_to_cam = vec3_scale(b, -b_recip_z);
        Vec3 c_to_cam = vec3_scale(c, -c_recip_z);

        a = persp(a);
        b = persp(b);
        c = persp(c);

        Vec3 face_norm = cross(vec3_sub(b, a), vec3_sub(c, a));
        if (face_norm.z <= 0) {
                return;
        }

        Vec3 ta = vec3_scale(vt[f->v0[UV]], a_recip_z);
        Vec3 tb = vec3_scale(vt[f->v1[UV]], b_recip_z);
        Vec3 tc = vec3_scale(vt[f->v2[UV]], c_recip_z);

        Vec3 na = vec3_scale(m4v3_mul(cam->inv_tr, vn[f->v0[NORM]]), a_recip_z);
        Vec3 nb = vec3_scale(m4v3_mul(cam->inv_tr, vn[f->v1[NORM]]), b_recip_z);
        Vec3 nc = vec3_scale(m4v3_mul(cam->inv_tr, vn[f->v2[NORM]]), c_recip_z);

        a_to_light = vec3_scale(a_to_light, a_recip_z);
        b_to_light = vec3_scale(b_to_light, b_recip_z);
        c_to_light = vec3_scale(c_to_light, c_recip_z);

        a_lightview = vec3_scale(a_lightview, a_recip_z);
        b_lightview = vec3_scale(b_lightview, b_recip_z);
        c_lightview = vec3_scale(c_lightview, c_recip_z);

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

        #define INTERP(r, s, t)     (recip_area * (Ea*(r) + Eb*(s) + Ec*(t)))
        #define VEC_INTERP(r, s, t) (vec3_scale(vec3_add(vec3_scale(r, Ea), vec3_add(vec3_scale(s, Eb), \
                                        vec3_scale(t, Ec))), recip_area))

        double LaKa = light->ambient * mat->ambient;
        double LdKd = light->diffuse * mat->diffuse;
        double LsKs = light->specular * mat->specular;

        for (int y = max_y; y >= min_y; --y) {
                for (int x = min_x; x <= max_x; ++x) {
                        if (Ea >= 0 && Eb >= 0 && Ec >= 0) {
                                double recip_z = INTERP(a_recip_z, b_recip_z, c_recip_z);
                                if (recip_z > dbuf_read(x, y)) {
                                        continue;
                                }

                                Vec3 lightview = VEC_INTERP(a_lightview, b_lightview, c_lightview);
                                lightview = vec3_scale(lightview, 1.0/recip_z);
                                lightview = m4v3_mul(light->viewport, persp(lightview));

                                Vec3 uv = vec3_scale(VEC_INTERP(ta, tb, tc), 1.0/recip_z);
                                Vec4 lit_colour = vec4_mul(sample_texture(model->texture, uv.x, uv.y), light->colour);

                                if (smap_read(light->shadow_map, lightview.x, lightview.y) - lightview.z <= 0.05) {
                                        Vec3 to_cam = vec3_scale(VEC_INTERP(a_to_cam, b_to_cam, c_to_cam), 1.0/recip_z);
                                        Vec3 to_light = vec3_scale(
                                                VEC_INTERP(a_to_light, b_to_light, c_to_light),
                                                1.0/recip_z
                                        );

                                        Vec3 light_dir = unit(to_light);
                                        Vec3 cam_dir = unit(to_cam);

                                        double dist = vec3_norm(to_light);
                                        double atten = 1.0/(1.0 + light->dropoff*dist*dist);
                                        Vec3 norm = unit(vec3_scale(VEC_INTERP(na, nb, nc), 1.0/recip_z));

                                        double ambient = LaKa;
                                        double diffuse = atten * LdKd * MAX(0, vec3_dot(norm, light_dir));
                                        double specular = atten * LsKs * pow(MAX(0, vec3_dot(reflect(light_dir, norm), cam_dir)), mat->shininess); // @TODO LUT for powers


                                        lit_colour = vec4_scale(lit_colour, ambient + diffuse);
                                        lit_colour = vec4_add(lit_colour, vec4_scale(light->colour, specular));

                                        Vec4 filtered_colour = fmap_read(light->filter_map, lightview.x, lightview.y);
                                        lit_colour = filtered_colour = vec4_mul(filtered_colour, lit_colour);
                                        //lit_colour = alpha_blend(filtered_colour, lit_colour);
                                } else {
                                        lit_colour = vec4_scale(lit_colour, LaKa);
                                }

                                point(x, y, lit_colour);
                                dbuf_write(x, y, recip_z);
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
// render_alpha_face - renders a single, textured & lit translucent face from a model
// @TODO lift inner loops of this and render_face to two separate functions (check that the per-face branch involved isn't too bad)
//
static void render_alpha_face(Face *f, Model *model, Camera *cam, Light *light) {
        Vec3 *v = model->obj->vertices;
        Vec3 *vt = model->obj->uvs;
        Vec3 *vn = model->obj->normals;

        Material *mat = model->material;

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

        Vec3 light_camview = m4v3_mul(cam->view_mat, light->pos);
        Vec3 a_to_light = vec3_sub(light_camview, a);
        Vec3 b_to_light = vec3_sub(light_camview, b);
        Vec3 c_to_light = vec3_sub(light_camview, c);

        double a_recip_z = 1.0/a.z;
        double b_recip_z = 1.0/b.z;
        double c_recip_z = 1.0/c.z;

        Vec3 a_to_cam = vec3_scale(a, -a_recip_z);
        Vec3 b_to_cam = vec3_scale(b, -b_recip_z);
        Vec3 c_to_cam = vec3_scale(c, -c_recip_z);

        a = persp(a);
        b = persp(b);
        c = persp(c);

        Vec3 face_norm = cross(vec3_sub(b, a), vec3_sub(c, a));
        if (face_norm.z <= 0) {
                return;
        }

        Vec3 ta = vec3_scale(vt[f->v0[UV]], a_recip_z);
        Vec3 tb = vec3_scale(vt[f->v1[UV]], b_recip_z);
        Vec3 tc = vec3_scale(vt[f->v2[UV]], c_recip_z);

        Vec3 na = vec3_scale(m4v3_mul(cam->inv_tr, vn[f->v0[NORM]]), a_recip_z);
        Vec3 nb = vec3_scale(m4v3_mul(cam->inv_tr, vn[f->v1[NORM]]), b_recip_z);
        Vec3 nc = vec3_scale(m4v3_mul(cam->inv_tr, vn[f->v2[NORM]]), c_recip_z);

        a_to_light = vec3_scale(a_to_light, a_recip_z);
        b_to_light = vec3_scale(b_to_light, b_recip_z);
        c_to_light = vec3_scale(c_to_light, c_recip_z);

        a_lightview = vec3_scale(a_lightview, a_recip_z);
        b_lightview = vec3_scale(b_lightview, b_recip_z);
        c_lightview = vec3_scale(c_lightview, c_recip_z);

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

        #define INTERP(r, s, t)     (recip_area * (Ea*(r) + Eb*(s) + Ec*(t)))
        #define VEC_INTERP(r, s, t) (vec3_scale(vec3_add(vec3_scale(r, Ea), vec3_add(vec3_scale(s, Eb), \
                                        vec3_scale(t, Ec))), recip_area))

        double LaKa = light->ambient * mat->ambient;
        double LdKd = light->diffuse * mat->diffuse;
        double LsKs = light->specular * mat->specular;

        for (int y = max_y; y >= min_y; --y) {
                for (int x = min_x; x <= max_x; ++x) {
                        if (Ea >= 0 && Eb >= 0 && Ec >= 0) {
                                double recip_z = INTERP(a_recip_z, b_recip_z, c_recip_z);
                                if (recip_z > dbuf_read(x, y)) {
                                        continue;
                                }

                                Vec3 lightview = VEC_INTERP(a_lightview, b_lightview, c_lightview);
                                lightview = vec3_scale(lightview, 1.0/recip_z);
                                lightview = m4v3_mul(light->viewport, persp(lightview));

                                Vec3 uv = vec3_scale(VEC_INTERP(ta, tb, tc), 1.0/recip_z);
                                Vec4 lit_colour = vec4_mul(sample_texture(model->texture, uv.x, uv.y), light->colour);

                                Vec4 pixel = get_pixel(x, y);

                                if (smap_read(light->shadow_map, lightview.x, lightview.y) <= lightview.z + 0.05) {
                                        Vec3 to_cam = vec3_scale(VEC_INTERP(a_to_cam, b_to_cam, c_to_cam), 1.0/recip_z);
                                        Vec3 to_light = vec3_scale(
                                                VEC_INTERP(a_to_light, b_to_light, c_to_light),
                                                1.0/recip_z
                                        );

                                        Vec3 light_dir = unit(to_light);
                                        Vec3 cam_dir = unit(to_cam);

                                        double dist = vec3_norm(to_light);
                                        double atten = 1.0/(1.0 + light->dropoff*dist*dist);
                                        Vec3 norm = unit(vec3_scale(VEC_INTERP(na, nb, nc), 1.0/recip_z));

                                        double ambient = LaKa;
                                        double diffuse = atten * LdKd * MAX(0, vec3_dot(norm, light_dir));
                                        double specular = atten * LsKs * pow(MAX(0, vec3_dot(reflect(light_dir, norm), cam_dir)), mat->shininess); // @TODO LUT for powers

                                        double old_alpha = lit_colour.w;
                                        lit_colour = vec4_scale(lit_colour, ambient + diffuse);
                                        lit_colour = vec4_add(lit_colour, vec4_scale(light->colour, specular));
                                        lit_colour.w = old_alpha;

                                        lit_colour = vec4_clamp(alpha_blend(lit_colour, pixel));
                                } else {
                                        lit_colour = vec4_scale(alpha_blend(lit_colour, pixel), LaKa);
                                }

                                point(x, y, lit_colour);
                                dbuf_write(x, y, recip_z);
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
// fog - screen-space fog effect
//
void fog(double thickness, Vec4 colour) {
        for (int y = g_min_y; y <= g_max_y; ++y) {
                for (int x = g_min_x; x <= g_max_x; ++x) {
                        colour.w = thickness * (1 - MIN(1, -dbuf_read(x, y)));
                        point(x, y, alpha_blend(colour, get_pixel(x, y)));
                }
        }
}

//
// render_model - renders a model
//
void render_model(Model *model, Camera *cam, Light *light) {
        Obj *obj = model->obj;
        Face *f = obj->faces;

        for (int i = 0; i < obj->face_count; ++i) {
                render_face(&f[i], model, cam, light);
        }
}

//
// render_alpha_model - renders a model with transparency
//
void render_alpha_model(Model *model, Camera *cam, Light *light) {
        Obj *obj = model->obj;
        Face *f = obj->faces;

        for (int i = 0; i < obj->face_count; ++i) {
                render_alpha_face(&f[i], model, cam, light);
        }
}

//
// render_scene - renders all the Models in a Scene, lit by all the Lights in that Scene,
//                from the perspective of a Camera
//
void render_scene(Scene *scene, Camera *cam) {
        reset_dbuf();

        int smap_lifetime = scene->lights[0]->shadow_map->lifetime;
        int smap_age = ++scene->lights[0]->shadow_map->age;
        if (smap_age >= smap_lifetime) {
                reset_smap(scene->lights[0]->shadow_map);

                for (int i = 0; i < scene->model_count; ++i) {
                        render_model_smap(scene->models[i], scene->lights[0]);
                }

                for (int i = 0; i < scene->alpha_model_count; ++i) {
                        render_model_smap(scene->alpha_models[i], scene->lights[0]);
                }

                scene->lights[0]->shadow_map->age = 0;
        }

        int fmap_lifetime = scene->lights[0]->filter_map->lifetime;
        int fmap_age = ++scene->lights[0]->filter_map->age;
        if (fmap_age >= fmap_lifetime) {
                reset_fmap(scene->lights[0]->filter_map);

                for (int i = 0; i < scene->alpha_model_count; ++i) {
                        render_model_fmap(scene->alpha_models[i], scene->lights[0]);
                }

                scene->lights[0]->filter_map->age = 0;
        }

        for (int i = 0; i < scene->model_count; ++i) {
                render_model(scene->models[i], cam, scene->lights[0]);
        }

        for (int i = 0; i < scene->alpha_model_count; ++i) {
                render_alpha_model(scene->alpha_models[i], cam, scene->lights[0]);
        }

#if 0
        Filter_Map *fmap = scene->lights[0]->filter_map;
        for (int i = 0; i < fmap->width*fmap->height; ++i) {
                g_framebuffer[i] = vec4_to_xrgb(rgba_to_vec4(fmap->colours[i]));
        }
#endif
}

//
// render_text_int - called by render_text to render integers as text
//
static Vec3 render_text_int(Vec3 pos, Vec4 colour, double scale, Texture *fontset, int d) {
        if (d == 0) {
                return render_glyph(pos, colour, scale, fontset, '0');
        }

        if (d < 0) {
                pos = render_glyph(pos, colour, scale, fontset, '-');
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
                pos = render_glyph(pos, colour, scale, fontset, '0' + (d / place_val) % 10);
                place_val /= 10;
        }

        return pos;
}

//
// render_text_double - called by render_text to render doubles as text
//
static Vec3 render_text_double(Vec3 pos, Vec4 colour, double scale, Texture *fontset, double f) {
        if (f < 0) {
                pos = render_glyph(pos, colour, scale, fontset, '-');
                f = -f;
        }

        long whole = (long)f;

        if (whole == 0) {
                pos = render_glyph(pos, colour, scale, fontset, '0');
        }

        int place_val = 1;
        while (place_val < whole) {
                place_val *= 10;
        }

        if (place_val > whole) {
                place_val /= 10;
        }

        while (place_val >= 1) {
                pos = render_glyph(pos, colour, scale, fontset, '0' + (whole / place_val) % 10);
                place_val /= 10;
        }

        pos = render_glyph(pos, colour, scale, fontset, '.');

        int max_digits = 2;
        do {
                f *= 10;
                pos = render_glyph(pos, colour, scale, fontset, '0' + (long)f % 10);
        } while (f > (long)f && max_digits-- > 1);

        return pos;
}

//
// render_text_vector - called by render_text to render a vector as text
//
static Vec3 render_text_vector(Vec3 pos, Vec4 colour, double scale, Texture *fontset, const Vec3 *v) {
        pos = render_glyph(pos, colour, scale, fontset, '(');
        pos = render_text_double(pos, colour, scale, fontset, v->x);
        pos = render_glyph(pos, colour, scale, fontset, ',');
        pos = render_glyph(pos, colour, scale, fontset, ' ');
        pos = render_text_double(pos, colour, scale, fontset, v->y);
        pos = render_glyph(pos, colour, scale, fontset, ',');
        pos = render_glyph(pos, colour, scale, fontset, ' ');
        pos = render_text_double(pos, colour, scale, fontset, v->z);
        pos = render_glyph(pos, colour, scale, fontset, ')');
        return pos;
}

//
// render_text - renders formatted text
//
Vec3 render_text(Vec3 pos, Vec4 colour, double scale, Texture *fontset, const char *fmt, ...) {
        va_list arglist;
        va_start(arglist, fmt);

        if (pos.y >= 1.0 || pos.x <= -1.0)  {
                goto _exit;
        }

        Vec3 curs = pos;
        for (int i = 0; fmt[i] != '\0'; ++i) {
                if (fmt[i] == '\n') {
                        curs.y -= (128 * scale)/g_max_y;
                        curs.x = pos.x;
                } else if (fmt[i] == '%') {
                        ++i;
                        if (fmt[i] == 'd') {
                                curs = render_text_int(curs, colour, scale, fontset, va_arg(arglist, int));
                        } else if (fmt[i] == 'f') {
                                curs = render_text_double(curs, colour, scale, fontset, va_arg(arglist, double));
                        } else if (fmt[i] == 'v') {
                                curs = render_text_vector(curs, colour, scale, fontset, va_arg(arglist, Vec3*));
                        } else if (fmt[i] == 's') { // @TODO render_text_string
                                curs = render_text(curs, colour, scale, fontset, va_arg(arglist, char*));
                        } else if (fmt[i] == 'c') {
                                curs = render_glyph(curs, colour, scale, fontset, va_arg(arglist, int));
                        } else if (fmt[i] == '%') {
                                curs = render_glyph(curs, colour, scale, fontset, '%');
                        }
                } else {
                        curs = render_glyph(curs, colour, scale, fontset, fmt[i]);
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
// load_obj - loads an OBJ file and parses geometry data into vertex arrays etc.
//
static Obj *load_obj(char *filename) {
        FILE *obj_file = fopen(filename, "rb");
        if (!obj_file) {
                return &g_error_obj;
        }

        fseek(obj_file, 0, SEEK_END);
        size_t obj_file_len = ftell(obj_file);
        rewind(obj_file);
        char *src = arena_alloc(&g_arena, obj_file_len + 1);
        fread(src, 1, obj_file_len, obj_file);
        fclose(obj_file);

        Obj *obj = arena_alloc(&g_arena, sizeof(Obj));

        size_t vert_buf_size = 1 << 12;
        obj->vertices = arena_alloc(&g_arena, vert_buf_size);

        #define SKIP_SPACE    {while (is_space(*src)) ++src;}
        #define SKIP_TO_SPACE {while (!is_space(*src)) ++src;}

        int i = 0;
        for (; *src == 'v' && src[1] == ' '; ++i) {
                src += 2;
                if (i * sizeof(Vec3) >= vert_buf_size) {
                        arena_commit_at(&g_arena, obj->vertices, vert_buf_size *= 2);
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
        obj->uvs = arena_alloc(&g_arena, uv_buf_size);

        i = 0;
        for (; *src == 'v' && src[1] == 't'; ++i) {
                src += 2;
                if (i * sizeof(Vec3) >= uv_buf_size) {
                        arena_commit_at(&g_arena, obj->uvs, uv_buf_size *= 2);
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
        obj->normals = arena_alloc(&g_arena, norm_buf_size);

        i = 0;
        for (; *src == 'v' && src[1] == 'n'; ++i) {
                src += 2;
                if (i * sizeof(Vec3) >= norm_buf_size) {
                        arena_commit_at(&g_arena, obj->normals, norm_buf_size *= 2);
                }

                SKIP_SPACE;
                obj->normals[i].x = parse_double(&src);
                SKIP_SPACE;
                obj->normals[i].y = parse_double(&src);
                SKIP_SPACE;
                obj->normals[i].z = parse_double(&src);
                SKIP_SPACE;
        }
        obj->norm_count = i;

        while (*src != 'f') {
                ++src;
        }

        size_t face_buf_size = 1 << 12;
        obj->faces = arena_alloc(&g_arena, face_buf_size);

        i = 0;
        for (; *src == 'f'; ++i) {
                ++src;
                if (i * sizeof(Face) >= face_buf_size) {
                        arena_commit_at(&g_arena, obj->faces, face_buf_size *= 2);
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
// load_texture - loads a texture file
//
Texture *load_texture(char *filename, size_t w, size_t h) {
        Texture *texture = arena_alloc(&g_arena, sizeof(Texture));
        texture->width = w;
        texture->height = h;

        FILE *tex_file = fopen(filename, "rb");
        if (!tex_file) {
                return &g_error_tex;
        }

        fseek(tex_file, 0, SEEK_END);
        size_t tex_file_len = ftell(tex_file);
        rewind(tex_file);

        texture->data = arena_alloc(&g_arena, tex_file_len);
        fread(texture->data, 1, tex_file_len, tex_file);
        fclose(tex_file);

        return texture;
}

//
// load_model - populates a Model by loading specified OBJ, texture, and normal map files for the model,
//              falls back to error assets on failure
//
Model *load_model(
        char *obj_filename,
        char *tm_filename,
        char *nm_filename,
        Material *mat,
        size_t tm_w,
        size_t tm_h,
        size_t nm_w,
        size_t nm_h
) {
        Model *model = arena_alloc(&g_arena, sizeof(Model));

        model->obj = load_obj(obj_filename);
        model->texture = load_texture(tm_filename, tm_w, tm_h);
        model->norm_map = load_texture(nm_filename, nm_w, nm_h);
        model->material = (mat == NULL) ? &g_error_material : mat;

        return model;
}

/*
 * @XXX misc. @XXX
 */

size_t get_mem_usage(void) {
        return arena_get_usage(&g_arena);
}

//
// init_scene - initialises a Scene by allocating lights and models arrays
//
void init_scene(Scene *scene, size_t model_limit, size_t alpha_model_limit, size_t light_limit) {
        scene->model_count = scene->alpha_model_count = scene->light_count = 0;

        scene->alpha_model_limit = alpha_model_limit;
        scene->model_limit = model_limit;
        scene->light_limit = light_limit;

        scene->alpha_models = arena_alloc(&g_arena, alpha_model_limit * sizeof(Model*));
        scene->models = arena_alloc(&g_arena, model_limit * sizeof(Model*));
        scene->lights = arena_alloc(&g_arena, light_limit * sizeof(Light*));
}

//
// add_model - adds a Model to a Scene
//
void add_model(Scene *scene, Model *model) {
        if (scene->model_count == scene->model_limit) {
                return;
        }

        if (model == NULL) {
                model = load_model("", "", NULL, NULL, 0, 0, 0, 0);
        }

        scene->models[scene->model_count] = model;
        ++scene->model_count;
}

//
// add_model - adds a Model to be rendered with transparency to a Scene
//
void add_alpha_model(Scene *scene, Model *model) {
        if (scene->alpha_model_count == scene->alpha_model_limit) {
                return;
        }

        if (model == NULL) {
                model = load_model("", "", NULL, NULL, 0, 0, 0, 0);
        }

        scene->alpha_models[scene->alpha_model_count] = model;
        ++scene->alpha_model_count;
}

//
// add_light - adds a Light to a Scene
//
void add_light(Scene *scene, Light *light) {
        if (scene->light_count == scene->light_limit) {
                return;
        }

        if (light == NULL) {
                light = &g_error_light;
        }

        scene->lights[scene->light_count] = light;
        ++scene->light_count;
}
