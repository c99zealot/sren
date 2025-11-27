#ifndef ARENA_H_INCLUDED
        #define ARENA_H_INCLUDED 1

        #include <stddef.h>
        
        typedef struct {
                void *base;
                size_t unused;
                size_t uncommitted;
                size_t unreserved;
                size_t alignment;
        } Arena;

        enum {ARENA_RESERVE_DEFAULT = 0};

        extern void *arena_init(Arena *arena, size_t init_reserve, size_t init_commit, size_t alignment);
        extern void arena_deinit(Arena *arena);
        extern void *arena_alloc(Arena *arena, size_t nobjs);
        extern void arena_clear(Arena *arena);
        extern void arena_commit_at(Arena *arena, void *ptr, size_t nobjs);
        extern void arena_debump(Arena *arena, size_t nobjs);
#endif
