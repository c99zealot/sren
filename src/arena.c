// @TODO make this work with LeakSanitizer

#include <stdint.h>
#include <string.h>

#include <sys/mman.h>

#ifdef ASAN
        #include <sanitizer/asan_interface.h>
#endif

#include "arena.h"
#include "ASSERT.h"

#define VOID_PTR_ADD(p, off) ((void*)((uintptr_t)(p) + (off)))
#define VOID_PTR_DIFF(p, q) ((uintptr_t)(p) - (uintptr_t)(q))

#define DEFAULT_RESERVE (1L << 36) // 64 GiB
enum {
        PAGE_SIZE = 4096,
        DEFAULT_ALIGNMENT = 8,
};

//
// alignto - aligns n up to alignment
//
static inline size_t alignto(size_t n, size_t alignment) {
        return (n + alignment - 1) & ~(alignment - 1);
}

//
// arena_init - initialises an Arena by reserving and committing memory,
//              returns the lowest address in the reserved range, asserts on error
//
void *arena_init(Arena *arena, size_t init_reserve, size_t init_commit, size_t alignment) {
        init_reserve = (init_reserve == ARENA_RESERVE_DEFAULT) ? DEFAULT_RESERVE :
                alignto(init_reserve * alignment, PAGE_SIZE);
        init_commit  = alignto(init_commit * alignment, PAGE_SIZE);

        arena->base = mmap(NULL, init_reserve, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        ASSERTF(arena->base != MAP_FAILED, "failed to reserve memory for arena (0x%X)\n", arena);

        arena->unreserved = init_reserve;
        arena->unused = arena->uncommitted = init_commit;
        arena->alignment = (alignment == 0) ? DEFAULT_ALIGNMENT : alignment;

        int err = mprotect(arena->base, init_commit, PROT_READ | PROT_WRITE);
        ASSERTF(err != -1, "failed to commit memory (0x%X) for arena (0x%X)\n", arena->base, arena);

#ifdef ASAN
        __asan_poison_memory_region((uint8_t*)arena->base + arena->unused, arena->uncommitted - arena->unused);
#endif

        return arena->base;
}

//
// arena_deinit - deinitialises an Arena by unmapping its memory and setting all its members to zero
//
void arena_deinit(Arena *arena) {
#ifdef ASAN
        __asan_unpoison_memory_region(arena->base, arena->uncommitted);
#endif

        int err = munmap(arena->base, arena->unreserved);
        ASSERTF(err != -1, "failed to unmap memory (0x%X) for arena (0x%X)\n", arena->base, arena);
        memset(arena, 0, sizeof(Arena));
}

//
// arena_clear - clears an Arena
//
void arena_clear(Arena *arena) {
        int err = mprotect(arena->base, arena->unused, PROT_NONE);
        ASSERTF(err != -1, "failed to change protection of memory (0x%X) for arena (0x%X)\n", arena->base, arena);
        arena->uncommitted = arena->unused = 0;
}

//
// arena_alloc - returns a pointer to an available region within the committed memory of an Arena
//               capable of storing nobjs objects of size arena->alignment bytes, pointer is aligned
//               to arena->alignment
//
void *arena_alloc(Arena *arena, size_t nobjs) {
        size_t bytes = nobjs * arena->alignment;
        size_t alloc_size = alignto(bytes, PAGE_SIZE);

        ASSERTF(arena->uncommitted + bytes < arena->unreserved, "can't allocate object; arena (0x%X) is full\n", arena);

        if (arena->unused + alloc_size > arena->uncommitted) {
                int err = mprotect((uint8_t*)arena->base + arena->uncommitted, alloc_size, PROT_READ | PROT_WRITE);
                ASSERTF(err != -1, "failed to commit (%zu bytes @ base+0x%X) for arena (0x%X)\n", alloc_size,
                        arena->uncommitted, arena);
                arena->uncommitted += alloc_size;
        }

#ifdef ASAN
        __asan_unpoison_memory_region((uint8_t*)arena->base + arena->unused, alloc_size);
#endif

        void *alloc_base = (uint8_t*)arena->base + arena->unused;
        arena->unused += bytes;

        return alloc_base;
}

//
// arena_commit_at - takes an address which was previously returned by arena_alloc and commits pages within the range
//                   [ptr, ptr + nobjs * arena->alignment]
//
void arena_commit_at(Arena *arena, void *ptr, size_t nobjs) {
        size_t bytes = nobjs * arena->alignment;
        size_t alloc_size = alignto(bytes, PAGE_SIZE);
        size_t alloc_offset = (uintptr_t)ptr - (uintptr_t)arena->base;
        size_t next_page_offset = alignto(alloc_offset, PAGE_SIZE);

        ASSERTF(alloc_offset + bytes < arena->unreserved, "can't allocate object; arena (0x%X) is full\n", arena);

        if (next_page_offset + alloc_size >= arena->uncommitted) {
                int err = mprotect((uint8_t*)arena->base + next_page_offset, alloc_size, PROT_READ | PROT_WRITE);
                ASSERTF(err != -1, "failed to commit (%zu bytes @ base+0x%X) for arena (0x%X)\n", alloc_size,
                        next_page_offset, arena);
                arena->uncommitted = next_page_offset + alloc_size;
        }

#ifdef ASAN
        __asan_unpoison_memory_region((uint8_t*)arena->base + next_page_offset, alloc_size);
#endif

        if (arena->unused < alloc_offset + bytes) {
                arena->unused = alloc_offset + bytes;
        }
}

//
// arena_debump - reduces the portion of committed memory which is in use by an Arena by freeing enough
//                space at the end of that portion to store nobjs objects of size arena->alignment
//
void arena_debump(Arena *arena, size_t nobjs) {
        if (arena->unused > 0) {
                arena->unused -= nobjs * arena->alignment;
        }

#ifdef ASAN
        __asan_poison_memory_region((uint8_t*)arena->base + arena->unused, nobjs * arena->alignment);
#endif
}

//
// arena_get_usage - returns the number of bytes allocated in an Arena by arena_alloc
//
size_t arena_get_usage(Arena *arena) {
        return arena->unused;
}
