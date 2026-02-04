#ifndef SREN_H_INCLUDED
        #define SREN_H_INCLUDED 1

        #include "arena.h"

        #include <stdint.h>
        #include <string.h>
        #include <stddef.h>
        #include <math.h>

        #define VEC3(x, y, z) ((Vec3){(x), (y), (z)})
        #define VEC4(x, y, z, w) ((Vec4){(x), (y), (z), (w)})

        #define MIN(x, y) ((x) < (y) ? (x) : (y))
        #define MAX(x, y) ((x) > (y) ? (x) : (y))
        #define MAX3(x, y, z) (MAX(MAX(x, y), z))
        #define MIN3(x, y, z) (MIN(MIN(x, y), z))
        #define CLAMP(x, y, z) (MAX(x, MIN(y, z)))
        
        #define RGB(x, y, z)  (((x) << 16) | ((y) << 8) | (z))
        #define RGBf(x, y, z) (((uint8_t)(x*255) << 16) | ((uint8_t)(y*255) << 8) | (uint8_t)(z*255))

        #define fmt_vec3 "(%f, %f, %f)"
        #define fmt_vec4 "(%f, %f, %f, %f)"
        #define fmt_mat3 "%f, %f, %f\n%f, %f, %f\n%f, %f, %f"
        #define fmt_mat4 "%f, %f, %f, %f\n%f, %f, %f, %f\n%f, %f, %f, %f\n%f, %f, %f, %f"

        enum {
                VERTEX = 0,
                UV     = 1,
                NORM   = 2,
        };

        typedef struct {
                double x, y, z;
        } Vec3;

        typedef struct {
                double x, y, z, w;
        } Vec4;

        typedef double Mat3[3][3];
        typedef double Mat4[4][4];

        typedef size_t Vertex_Data[3];

        typedef struct {
                Vertex_Data v0;
                Vertex_Data v1;
                Vertex_Data v2;
        } Face;

        typedef struct {
                size_t width;
                size_t height;

                uint8_t *data;
        } Texture;

        typedef struct {
                double *data;

                size_t width;
                size_t height;

                int max_x;
                int max_y;
                int min_x;
                int min_y;
        } Shadow_Map;

        typedef struct {
                Mat4 view_mat;
                Mat4 inv_tr;

                Vec3 pos;
                Vec3 subject;
                Vec3 up;
                Vec3 l;
                Vec3 m;
                Vec3 n;
        } Camera;

        typedef struct {
                Mat4 view_mat;
                Mat4 viewport;

                Vec3 pos;
                Vec3 subject;
                Vec3 up;

                Vec4 colour;

                Shadow_Map *shadow_map;

                double ambient;
                double diffuse;
                double specular;

                double dropoff;
        } Light;

        typedef struct {
                Vec3 *vertices;
                Vec3 *uvs;
                Vec3 *normals;

                Face *faces;

                size_t vertex_count;
                size_t uv_count;
                size_t norm_count;
                size_t face_count;
        } Obj;

        typedef struct {
                double ambient;
                double diffuse;
                double specular;
                int shininess;
        } Material;

        typedef struct {
                Obj *obj;

                Texture *texture;
                Texture *norm_map;

                Material *material;
        } Model;

        typedef struct {
                Model **models;

                Light **lights;

                size_t model_count;
                size_t light_count;
        } Scene;

        extern Mat4 g_viewport; // @XXX exposing this for now...

        //
        // cvec3s_to_mat3 - constructs a 3x3 matrix from three 3D column vectors
        //
        static inline void cvec3s_to_mat3(Mat3 m, Vec3 *v0, Vec3 *v1, Vec3 *v2) {
                m[0][0] = v0->x;
                m[0][1] = v1->x;
                m[0][2] = v2->x;

                m[1][0] = v0->y;
                m[1][1] = v1->y;
                m[1][2] = v2->y;

                m[2][0] = v0->z;
                m[2][1] = v1->z;
                m[2][2] = v2->z;
        }

        //
        // rvec3s_to_mat3 - constructs a 3x3 matrix from three 3D row vectors
        //
        static inline void rvec3s_to_mat3(Mat3 m, Vec3 *v0, Vec3 *v1, Vec3 *v2) {
                memcpy(*m, v0, sizeof(Vec3));
                memcpy(*m + 3, v1, sizeof(Vec3));
                memcpy(*m + 6, v2, sizeof(Vec3));
        }

        //
        // cvec3s_to_hmat4 - constructs a 4x4 matrix from three homogenised 3D column vectors
        //
        static inline void cvec3s_to_hmat4(Mat4 m, Vec3 *v0, Vec3 *v1, Vec3 *v2) {
                m[0][0] = v0->x;
                m[0][1] = v1->x;
                m[0][2] = v2->x;
                m[0][3] = 0.0;

                m[1][0] = v0->y;
                m[1][1] = v1->y;
                m[1][2] = v2->y;
                m[1][3] = 0.0;

                m[2][0] = v0->z;
                m[2][1] = v1->z;
                m[2][2] = v2->z;
                m[2][3] = 0.0;

                m[3][0] = m[3][1] = m[3][2] = 0.0;
                m[3][3] = 1.0;
        }

        //
        // rvec3s_to_hmat4 - constructs a 4x4 matrix from three homogenised 3D row vectors
        //
        static inline void rvec3s_to_hmat4(Mat4 m, Vec3 *v0, Vec3 *v1, Vec3 *v2) {
                memcpy(*m, v0, sizeof(Vec3));
                memcpy(*m + 4, v1, sizeof(Vec3));
                memcpy(*m + 8, v2, sizeof(Vec3));

                m[0][3] = 0.0;
                m[1][3] = 0.0;
                m[2][3] = 0.0;

                m[3][0] = 0.0;
                m[3][1] = 0.0;
                m[3][2] = 0.0;

                m[3][3] = 1.0;
        }

        //
        // mv3_mul - computes the product of a 3x3 matrix and a 3D vector
        //
        static inline Vec3 mv3_mul(Mat3 m, Vec3 v) {
                Mat3 _m;
                memcpy(_m, m, sizeof(Mat3));

                _m[0][0] *= v.x;
                _m[0][1] *= v.y;
                _m[0][2] *= v.z;

                _m[1][0] *= v.x;
                _m[1][1] *= v.y;
                _m[1][2] *= v.z;

                _m[2][0] *= v.x;
                _m[2][1] *= v.y;
                _m[2][2] *= v.z;

                v.x = _m[0][0] + _m[0][1] + _m[0][2];
                v.y = _m[1][0] + _m[1][1] + _m[1][2];
                v.z = _m[2][0] + _m[2][1] + _m[2][2];

                return v;
        }

        //
        // m4v3_mul - computes the product of a 4x4 matrix and a homogenised 3D vector
        //
        static inline Vec3 m4v3_mul(Mat4 m, Vec3 v) {
                Mat4 _m;
                memcpy(_m, m, sizeof(Mat4));

                Vec4 _v = {.x = v.x, .y = v.y, .z = v.z, .w = 1.0};

                _m[0][0] *= _v.x;
                _m[0][1] *= _v.y;
                _m[0][2] *= _v.z;
                _m[0][3] *= _v.w;

                _m[1][0] *= _v.x;
                _m[1][1] *= _v.y;
                _m[1][2] *= _v.z;
                _m[1][3] *= _v.w;

                _m[2][0] *= _v.x;
                _m[2][1] *= _v.y;
                _m[2][2] *= _v.z;
                _m[2][3] *= _v.w;

                v.x = _m[0][0] + _m[0][1] + _m[0][2] + _m[0][3];
                v.y = _m[1][0] + _m[1][1] + _m[1][2] + _m[1][3];
                v.z = _m[2][0] + _m[2][1] + _m[2][2] + _m[2][3];

                return v;
        }

        //
        // mv4_mul - computes the product of a 4x4 matrix and a 4D vector
        //
        static inline Vec4 mv4_mul(Mat4 m, Vec4 v) {
                Mat4 _m;
                memcpy(_m, m, sizeof(Mat4));

                _m[0][0] *= v.x;
                _m[0][1] *= v.y;
                _m[0][2] *= v.z;
                _m[0][3] *= v.w;

                _m[1][0] *= v.x;
                _m[1][1] *= v.y;
                _m[1][2] *= v.z;
                _m[1][3] *= v.w;

                _m[2][0] *= v.x;
                _m[2][1] *= v.y;
                _m[2][2] *= v.z;
                _m[2][3] *= v.w;

                v.x = _m[0][0] + _m[0][1] + _m[0][2] + _m[0][3];
                v.y = _m[1][0] + _m[1][1] + _m[1][2] + _m[1][3];
                v.z = _m[2][0] + _m[2][1] + _m[2][2] + _m[2][3];
                v.w = _m[2][0] + _m[2][1] + _m[2][2] + _m[3][3];

                return v;
        }

        //
        // matmul4 - computes the product of two 4x4 matrices
        //
        static inline void matmul4(Mat4 C, Mat4 A, Mat4 B) {
                for (int i = 0; i < 4; ++i) {
                        double a0 = A[i][0], a1 = A[i][1], a2 = A[i][2], a3 = A[i][3];

                        C[i][0] = a0 * B[0][0] + a1 * B[1][0] + a2 * B[2][0] + a3 * B[3][0];
                        C[i][1] = a0 * B[0][1] + a1 * B[1][1] + a2 * B[2][1] + a3 * B[3][1];
                        C[i][2] = a0 * B[0][2] + a1 * B[1][2] + a2 * B[2][2] + a3 * B[3][2];
                        C[i][3] = a0 * B[0][3] + a1 * B[1][3] + a2 * B[2][3] + a3 * B[3][3];
                }
        }

        //
        // m4_inverse_transpose - computes the inverse transpose of a 4x4 matrix
        //
        static inline void m4_inverse_transpose(Mat4 out, Mat4 m) {
                Mat4 cof;
                for (int r = 0; r < 4; r++) {
                        for (int c = 0; c < 4; c++) {

                                double sub[3][3];
                                int sr = 0, sc;

                                for (int i = 0; i < 4; i++) {
                                        if (i == r) continue;
                                        sc = 0;
                                        for (int j = 0; j < 4; j++) {
                                                if (j == c) continue;
                                                sub[sr][sc] = m[i][j];
                                                sc++;
                                        }
                                        sr++;
                                }

                                double sign = ((r + c) & 1) ? -1.0 : 1.0;

                                cof[r][c] = sign * (
                                        sub[0][0] * (sub[1][1] * sub[2][2] - sub[1][2] * sub[2][1]) -
                                        sub[0][1] * (sub[1][0] * sub[2][2] - sub[1][2] * sub[2][0]) +
                                        sub[0][2] * (sub[1][0] * sub[2][1] - sub[1][1] * sub[2][0])
                                );
                        }
                }

                double det = 0.0;
                for (int j = 0; j < 4; j++) {
                        det += m[0][j] * cof[0][j];
                }

                double recip_det = 1.0/det;
                for (int r = 0; r < 4; r++) {
                        for (int c = 0; c < 4; c++) {
                                out[r][c] = cof[r][c] * recip_det;
                        }
                }
        }

        //
        // vec4_mul - computes the Hadamard product of two 4D vectors
        //
        static inline Vec4 vec4_mul(Vec4 v0, Vec4 v1) {
                return VEC4(v0.x*v1.x, v0.y*v1.y, v0.z*v1.z, v0.w*v1.w);
        }

        //
        // vec3_mul - computes the Hadamard product of two 3D vectors
        //
        static inline Vec3 vec3_mul(Vec3 v0, Vec3 v1) {
                return VEC3(v0.x*v1.x, v0.y*v1.y, v0.z*v1.z);
        }

        //
        // vec3_dot - computes the dot product of two 3D vectors
        //
        static inline double vec3_dot(Vec3 v0, Vec3 v1) {
                return v0.x*v1.x + v0.y*v1.y + v0.z*v1.z;
        }

        //
        // vec3_norm - computes the norm of a 3D vector
        //
        static inline double vec3_norm(Vec3 v) {
                return sqrt(vec3_dot(v, v));
        }

        //
        // unit - returns a unit vector pointing in the direction of v
        //
        static inline Vec3 unit(Vec3 v) {
                double c = 1.0/vec3_norm(v);
                v.x *= c;
                v.y *= c;
                v.z *= c;
                return v;
        }

        //
        // cross - computes the cross product of two vectors
        //
        static inline Vec3 cross(Vec3 v0, Vec3 v1) {
                return (Vec3){
                        .x = v0.y*v1.z - v0.z*v1.y,
                        .y = v0.z*v1.x - v0.x*v1.z,
                        .z = v0.x*v1.y - v0.y*v1.x
                };
        }

        //
        // vec3_add - subtracts a 3D vector from another
        //
        static inline Vec3 vec3_sub(Vec3 v0, Vec3 v1) {
                return (Vec3){v0.x - v1.x, v0.y - v1.y, v0.z - v1.z};
        }

        //
        // vec3_add - adds two 3D vectors
        //
        static inline Vec3 vec3_add(Vec3 v0, Vec3 v1) {
                return (Vec3){v0.x + v1.x, v0.y + v1.y, v0.z + v1.z};
        }

        //
        // vec3_scale - scales a 3D vector v by a scalar s
        //
        static inline Vec3 vec3_scale(Vec3 v, double s) {
                return (Vec3){s*v.x, s*v.y, s*v.z};
        }

        //
        // vec4_add - adds two 4D vectors
        //
        static inline Vec4 vec4_add(Vec4 v0, Vec4 v1) {
                return (Vec4){v0.x + v1.x, v0.y + v1.y, v0.z + v1.z, v0.w + v1.w};
        }

        //
        // Vec4_scale - scales a 3D vector v by a scalar s
        //
        static inline Vec4 vec4_scale(Vec4 v, double s) {
                return (Vec4){s*v.x, s*v.y, s*v.z, s*v.w};
        }

        //
        // xrgb_to_vec4 - converts an XRGB8888 value to a Vec4
        //
        static inline Vec4 xrgb_to_vec4(uint32_t xrgb) {
                uint8_t r = ((uint8_t*)&xrgb)[2];
                uint8_t g = ((uint8_t*)&xrgb)[1];
                uint8_t b = ((uint8_t*)&xrgb)[0];

                return vec4_scale(VEC4(r, g, b, 1), 1.0/255);
        }

        //
        // xrgb_to_vec4 - converts a Vec4 to an XRGB8888 value
        //
        static inline uint32_t vec4_to_xrgb(Vec4 xrgb) {
                uint32_t result;
                uint8_t *p = (uint8_t*)&result;

                p[2] = (uint8_t)(xrgb.x*255);
                p[1] = (uint8_t)(xrgb.y*255);
                p[0] = (uint8_t)(xrgb.z*255);

                return result;
        }

        //
        // persp - basic perspective projection
        //
        static inline Vec3 persp(Vec3 v) {
                double recip_z = -1.0/v.z;
                return VEC3(v.x * recip_z, v.y * recip_z, v.z);
        }

        //
        // init_cam - initialises a Camera
        //
        static inline void init_cam(Camera *cam) {
                cam->n = unit(vec3_sub(cam->pos, cam->subject));
                cam->l = unit(cross(cam->up, cam->n));
                cam->m = cross(cam->n, cam->l);

                Mat4 view = {
                        cam->l.x,  cam->l.y,  cam->l.z, -vec3_dot(cam->l, cam->pos),
                        cam->m.x,  cam->m.y,  cam->m.z, -vec3_dot(cam->m, cam->pos),
                        cam->n.x,  cam->n.y,  cam->n.z, -vec3_dot(cam->n, cam->pos),
                        0.0,       0.0,       0.0,       1.0
                };

                memcpy(cam->view_mat, view, sizeof(Mat4));
                m4_inverse_transpose(cam->inv_tr, cam->view_mat);
        }

        //
        // move_cam - rotates and displaces a Camera
        //
        static inline void move_cam(Camera *cam, Mat3 rot, Vec3 disp) {
                Mat4 view_inv = {
                        cam->l.x, cam->m.x, cam->n.x, cam->pos.x,
                        cam->l.y, cam->m.y, cam->n.y, cam->pos.y,
                        cam->l.z, cam->m.z, cam->n.z, cam->pos.z,
                        0.0,      0.0,      0.0,      1.0
                };

                cam->subject = m4v3_mul(cam->view_mat, cam->subject);
                cam->subject = mv3_mul(rot, cam->subject);
                cam->subject = m4v3_mul(view_inv, cam->subject);

                cam->n = unit(vec3_sub(cam->pos, cam->subject));
                cam->l = unit(cross(cam->up, cam->n));
                cam->m = cross(cam->n, cam->l);

                disp = vec3_add(
                        vec3_scale(cam->l, disp.x),
                        vec3_add(
                                vec3_scale(cam->m, disp.y),
                                vec3_scale(cam->n, disp.z)
                        )
                );
                cam->pos = vec3_add(cam->pos, disp);
                cam->subject = vec3_add(cam->subject, disp);

                Mat4 view = {
                        cam->l.x, cam->l.y, cam->l.z, -vec3_dot(cam->l, cam->pos),
                        cam->m.x, cam->m.y, cam->m.z, -vec3_dot(cam->m, cam->pos),
                        cam->n.x, cam->n.y, cam->n.z, -vec3_dot(cam->n, cam->pos),
                        0.0,      0.0,      0.0,       1.0
                };

                memcpy(cam->view_mat, view, sizeof(Mat4));
                m4_inverse_transpose(cam->inv_tr, cam->view_mat);
        }

        //
        // move_light_to - moves a Light to a position in world-space and updates its view matrix
        //
        static inline void move_light_to(Light *light, Vec3 pos) {
                light->pos = pos;

                Vec3 n = unit(vec3_sub(light->pos, light->subject));
                Vec3 l = unit(cross(light->up, n));
                Vec3 m = cross(n, l);

                Mat4 view = {
                        l.x,  l.y,  l.z, -vec3_dot(l, light->pos),
                        m.x,  m.y,  m.z, -vec3_dot(m, light->pos),
                        n.x,  n.y,  n.z, -vec3_dot(n, light->pos),
                        0.0,  0.0,  0.0,  1.0
                };

                memcpy(light->view_mat, view, sizeof(Mat4));
        }

        extern void clear_screen(char c);
        extern void point(int x, int y, Vec4 colour);
        extern void init_renderer(uint32_t *framebuffer, size_t fb_width, size_t fb_height);
        extern void deinit_renderer(void);
        extern void line(Vec3 a, Vec3 b, Vec4 colour);
        extern void render_model(Model *model, Camera *cam, Light *light);
        extern double dbuf_read(int x, int y);
        extern void reset_dbuf(void);
        extern Vec3 render_text(Vec3 pos, Vec4 colour, double scale, Texture *fontset, const char *fmt, ...);
        extern void render_image_fragment(Vec3 pos, Vec4 colour, Texture *texture, Vec3 uv, double frag_w, double frag_h, double scale);
        extern Vec3 render_glyph(Vec3 pos, Vec4 colour, double scale, Texture *fontset, char c);
        extern void init_light(Light *light, size_t smap_width, size_t smap_height);
        extern int out_of_view(Vec3 v);
        extern void fog(double thickness, Vec4 colour);
        extern double smap_read(Shadow_Map *shadow_map, int x, int y);
        extern void render_model_smap(Model *model, Light *light);
        extern void reset_smap(Shadow_Map *smap);
        extern size_t get_mem_usage(void);
        extern Texture *load_texture(char *filename, size_t w, size_t h);
        extern Model *load_model(
                char *obj_filename,
                char *tm_filename,
                char *nm_filename,
                Material *mat,
                size_t tm_w,
                size_t tm_h,
                size_t nm_w,
                size_t nm_h
        );
#endif
