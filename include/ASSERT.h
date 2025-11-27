#ifndef ASSERT_H_INCLUDED
        #define ASSERT_H_INCLUDED 1

        #include <stdlib.h>
        #include <stdio.h>

        #ifdef NOP_ASSERTIONS
                #define ASSERT(x, msg) ((void)0)
        #else
                #define EXPAND(x) #x
                #define DELAY_EXPAND(x) EXPAND(x)

                #define ASSERT(x, msg) (x ? (void)0 : (                          \
                        fprintf(stderr,                                          \
                                "-- Failed assertion in " __FILE__               \
                                " --\n%s() {\n\t...\n" DELAY_EXPAND(__LINE__)    \
                                " | \tASSERT(" #x ", \"" msg "\");\n\t...\n}\n", \
                                __func__                                         \
                        ), fputs(msg "\n", stderr),                              \
                           exit(EXIT_FAILURE)))
        #endif

        #ifdef NOP_ASSERTIONS
                #define ASSERTF(x, ...) ((void)0)
        #else
                #undef EXPAND
                #undef DELAY_EXPAND

                #define EXPAND(x) #x
                #define DELAY_EXPAND(x) EXPAND(x)

                #define ASSERTF(x, ...) (x ? (void)0 : (                      \
                        fprintf(stderr,                                       \
                                "-- Failed assertion in " __FILE__            \
                                " --\n%s() {\n\t...\n" DELAY_EXPAND(__LINE__) \
                                " | \tASSERT(" #x ", ...);\n\t...\n}\n",      \
                                __func__                                      \
                        ), fprintf(stderr, __VA_ARGS__),                      \
                           putchar('\n'),                                     \
                           exit(EXIT_FAILURE)))
        #endif
#endif
