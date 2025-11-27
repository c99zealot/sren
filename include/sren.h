#ifndef SREN_H_INCLUDED
        #define SREN_H_INCLUDED 1

        #include "arena.h"

        #include <stdint.h>
        #include <stddef.h>

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
                Vec3 pos;
                Vec3 subject;
                Vec3 up;
                Mat4 view_mat;
                Mat4 inv_tr;
        } Camera;

        typedef struct {
                Vec3 pos;
                Vec3 colour;
                double intensity;
        } Light;

        typedef struct {
                Vec3 *vertices;
                Vec3 *uvs;
                Vec3 *norms;
                Face *faces;
                size_t vertex_count;
                size_t uv_count;
                size_t norm_count;
                size_t face_count;
        } Obj;

        typedef struct {
                size_t width;
                size_t height;
                uint8_t *data;
        } Texture;

        typedef struct {
                Obj *obj;
                Texture *texture;
                Texture *norm_map;
        } Model;

        typedef struct {
                Model **models;
                Light **lights;

                size_t model_count;
                size_t light_count;
        } Scene;

        inline void cvec3s_to_mat3(Mat3 m, Vec3 *v0, Vec3 *v1, Vec3 *v2) {
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

        inline void rvec3s_to_mat3(Mat3 m, Vec3 *v0, Vec3 *v1, Vec3 *v2) {
                memcpy(*m, v0, sizeof(Vec3));
                memcpy(*m + 3, v1, sizeof(Vec3));
                memcpy(*m + 6, v2, sizeof(Vec3));
        }

        inline void cvec3s_to_hmat4(Mat4 m, Vec3 *v0, Vec3 *v1, Vec3 *v2) {
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

        inline void rvec3s_to_hmat4(Mat4 m, Vec3 *v0, Vec3 *v1, Vec3 *v2) {
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

        inline Vec3 mv3_mul(Mat3 m, Vec3 v) {
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

        inline Vec3 m4v3_mul(Mat4 m, Vec3 v) {
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

        inline Vec4 mv4_mul(Mat4 m, Vec4 v) {
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

        inline void matmul4(Mat4 C, Mat4 A, Mat4 B) {
                for (int i = 0; i < 4; ++i) {
                        double a0 = A[i][0], a1 = A[i][1], a2 = A[i][2], a3 = A[i][3];

                        C[i][0] = a0 * B[0][0] + a1 * B[1][0] + a2 * B[2][0] + a3 * B[3][0];
                        C[i][1] = a0 * B[0][1] + a1 * B[1][1] + a2 * B[2][1] + a3 * B[3][1];
                        C[i][2] = a0 * B[0][2] + a1 * B[1][2] + a2 * B[2][2] + a3 * B[3][2];
                        C[i][3] = a0 * B[0][3] + a1 * B[1][3] + a2 * B[2][3] + a3 * B[3][3];
                }
        }

        inline void m4_inverse_transpose(Mat4 out, Mat4 m) {
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

        inline double vec3_dot(Vec3 v0, Vec3 v1) {
                return v0.x*v1.x + v0.y*v1.y + v0.z*v1.z;
        }

        inline double vec3_norm(Vec3 v) {
                return sqrt(vec3_dot(v, v));
        }

        inline Vec3 unit(Vec3 v) {
                double c = 1.0/vec3_norm(v);
                v.x *= c;
                v.y *= c;
                return v;
        }

        inline Vec3 cross(Vec3 v0, Vec3 v1) {
                return (Vec3){
                        .x = v0.y*v1.z - v0.z*v1.y,
                        .y = v0.z*v1.x - v0.x*v1.z,
                        .z = v0.x*v1.y - v0.y*v1.x
                };
        }

        inline Vec3 vec3_sub(Vec3 v0, Vec3 v1) {
                return (Vec3){v0.x - v1.x, v0.y - v1.y, v0.z - v1.z};
        }

        inline Vec3 vec3_add(Vec3 v0, Vec3 v1) {
                return (Vec3){v0.x + v1.x, v0.y + v1.y, v0.z + v1.z};
        }

        inline Vec3 vec3_scale(Vec3 v, double s) {
                return (Vec3){s*v.x, s*v.y, s*v.z};
        }

        inline double signed_tri_area2(Vec3 a, Vec3 b, Vec3 c) {
                return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
        }

        inline double attenuate(double intensity, double dist) {
                intensity = intensity < 0 ? 0 : intensity/dist;
                return intensity > 1 ? 1 : intensity;
        }

        inline void set_view(Camera *cam) {
                Vec3 n = unit(vec3_sub(cam->pos, cam->subject));
                Vec3 l = unit(cross(cam->up, n));
                Vec3 m = unit(cross(n, l));

                rvec3s_to_hmat4(cam->view_mat, &l, &m, &n);
                Mat4 C = {
                        1.0, 0.0, 0.0, -cam->pos.x,
                        0.0, 1.0, 0.0, -cam->pos.y,
                        0.0, 0.0, 1.0, -cam->pos.z,
                        0.0, 0.0, 0.0, 1.0
                };

                matmul4(cam->view_mat, cam->view_mat, C);
                m4_inverse_transpose(cam->inv_tr, cam->view_mat);
        }

        static inline Vec3 sample_texture(Texture *texture, double x, double y) {
                int u = x*(texture->width - 1);
                int v = y*(texture->height - 1);

                int v_flipped = texture->height - 1 - v;
                size_t i = v_flipped*texture->width + u;

                return (Vec3){
                        .x = (double)texture->data[i * 3*sizeof(uint8_t)]/255,
                        .y = (double)texture->data[i * 3*sizeof(uint8_t) + 1]/255,
                        .z = (double)texture->data[i * 3*sizeof(uint8_t) + 2]/255
                };
        }

        extern Texture *load_texture(Arena *arena, char *filename, size_t w, size_t h);
        extern Obj *load_obj(Arena *arena, char *filename);
        extern Model *load_model(
                Arena *arena,
                char *obj_filename,
                char *tm_filename,
                char *nm_filename,
                size_t tm_w,
                size_t tm_h,
                size_t nm_w,
                size_t nm_h
        );
#endif
