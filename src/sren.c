#include <stdio.h>
#include <string.h>
#include <math.h>

#include "sren.h"
#include "arena.h"

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

static int is_space(unsigned c) {
        return (c == ' ')      ||
               (c - '\a' <= 6);
}

Texture *load_texture(Arena *arena, char *filename, size_t w, size_t h) {
        Texture *texture = arena_alloc(arena, sizeof(Obj));
        texture->width = w;
        texture->height = h;

        FILE *tex_file = fopen(filename, "rb");
        if (!tex_file) {
                return NULL;
        }

        fseek(tex_file, 0, SEEK_END);
        size_t tex_file_len = ftell(tex_file);
        rewind(tex_file);

        texture->data = arena_alloc(arena, tex_file_len);
        fread(texture->data, 1, tex_file_len, tex_file);
        fclose(tex_file);

        return texture;
}

Obj *load_obj(Arena *arena, char *filename) {
        FILE *obj_file = fopen(filename, "rb");
        if (!obj_file) {
                return NULL;
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

        if (obj_filename != NULL) {
                model->obj = load_obj(arena, obj_filename);
        }

        if (tm_filename != NULL) {
                model->texture = load_texture(arena, tm_filename, tm_w, tm_h);
        }

        if (nm_filename != NULL) {
                model->norm_map = load_texture(arena, nm_filename, nm_w, nm_h);
        }

        return model;
}

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
