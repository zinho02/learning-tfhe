// poly.hpp -- Polynomials modulo X^N + 1.
//
// A polynomial of degree less than N is the vector of its N coefficients,
// where p[j] is the coefficient of X^j. Since X^N = -1, every term of degree
// N or more wraps around to the low degrees with its sign flipped (negacyclic
// product). Exponents of X therefore live modulo 2N, since X^(2N) = 1.
#pragma once
#include <vector>
#include "torus.hpp"

namespace tfhe {

using IntPoly = std::vector<int32_t>;
using TorusPoly = std::vector<Torus>;

// r = a * b mod X^N + 1 by the schoolbook method, O(N^2). Too slow for the
// scheme itself; it is the reference against which the FFT is tested.
inline TorusPoly mulNaive(const IntPoly& a, const TorusPoly& b) {
    int N = static_cast<int>(a.size());
    TorusPoly r(N, 0);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            Torus p = a[i] * b[j];
            if (i + j < N)  r[i + j] += p;
            else            r[i + j - N] -= p;
        }
    }
    return r;
}

// r = X^e * p mod X^N + 1, for 0 <= e < 2N. No multiplication is needed: it is
// a rotation of the coefficients by e positions, negating those that wrap
// around. For e >= N, X^e = -X^(e-N), i.e. the same rotation with all signs
// flipped.
inline TorusPoly mulByXe(const TorusPoly& p, int e) {
    int N = static_cast<int>(p.size());
    TorusPoly r(N);
    if (e < N) {
        for (int j = 0; j < e; ++j) r[j] = -p[j - e + N];
        for (int j = e; j < N; ++j) r[j] = p[j - e];
    } else {
        e -= N;
        for (int j = 0; j < e; ++j) r[j] = p[j - e + N];
        for (int j = e; j < N; ++j) r[j] = -p[j - e];
    }
    return r;
}

// r <- r + p and r <- r - p, coefficient-wise.
inline void polyAddTo(TorusPoly& r, const TorusPoly& p) {
    for (size_t j = 0; j < r.size(); ++j) r[j] += p[j];
}

inline void polySubTo(TorusPoly& r, const TorusPoly& p) {
    for (size_t j = 0; j < r.size(); ++j) r[j] -= p[j];
}

}  // namespace tfhe
