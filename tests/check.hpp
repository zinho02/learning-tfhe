// check.hpp -- Minimal test "framework": a single assertion macro.
//
// CHECK(cond) does nothing if cond holds. Otherwise it prints the condition,
// file and line, and exits with status 1, which makes "make test" stop at the
// first failing test.
#pragma once
#include <cstdio>
#include <cstdlib>

#define CHECK(cond) \
    do { \
        if (!(cond)) { \
            std::printf("FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__); \
            std::exit(1); \
        } \
    } while (0)
