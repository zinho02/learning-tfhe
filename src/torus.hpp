// torus.hpp -- The discretized torus T = R/Z and the sources of randomness.
//
// A torus element is stored as a 32-bit unsigned integer t that represents the
// real number t / 2^32 (mod 1). Unsigned arithmetic in C++ is defined modulo
// 2^32, so addition, subtraction and multiplication by an integer on uint32_t
// are exactly the torus operations, with the reduction mod 1 for free.
#pragma once
#include <cmath>
#include <cstdint>
#include <random>

namespace tfhe {

using Torus = uint32_t;

// 2^32 as a double: the scale factor between reals and torus integers.
#define MAX_INT 4294967296.0

// Real -> torus: reduces d modulo 1 and rounds it to the nearest multiple of
// 2^-32. A result of exactly 2^32 wraps to 0, which is the same torus point.
inline Torus dtot(double d) {
    double frac = d - std::floor(d);
    return static_cast<Torus>(std::llround(frac * MAX_INT));
}

// Torus -> real: returns the representative in [-1/2, 1/2), obtained by
// reading the 32 bits as a signed integer. This is the natural way to look at
// errors, since a small negative noise shows up as -0.001 and not as 0.999.
inline double ttod(Torus t) {
    return static_cast<int32_t>(t) / MAX_INT;
}

// Encodes mu in Z_M as the torus point mu / M (computed in 64 bits).
inline Torus encode(int32_t mu, int32_t M) {
    return static_cast<Torus>(mu * (int64_t{1} << 32) / M);
}

// Decodes a phase to the closest message in {0, ..., M-1}. Adding half an
// interval before the integer division turns truncation into rounding, and
// the final "% M" handles the wrap-around from M back to 0. Any error smaller
// than 1/(2M) is removed.
inline int32_t decode(Torus phase, int32_t M) {
    uint64_t interv = (uint64_t{1} << 32) / M;
    return static_cast<int32_t>(((phase + interv / 2) / interv) % M);
}

// Global pseudorandom generator, seeded from std::random_device.
// NOTE: mt19937_64 is NOT cryptographically secure. It is adequate for
// studying and measuring the scheme, but a real deployment needs a CSPRNG.
inline std::mt19937_64& rng() {
    static std::mt19937_64 gen{std::random_device{}()};
    return gen;
}

// Fixes the seed; the tests use it to be reproducible.
inline void setSeed(uint64_t s) { rng().seed(s); }

// Uniform torus element (used for the masks of ciphertexts).
inline Torus uniformTorus() { return static_cast<Torus>(rng()()); }

// Uniform bit (used for the secret keys).
inline int32_t uniformBit() { return static_cast<int32_t>(rng()() & 1); }

// Centered Gaussian noise with standard deviation stdev, given in torus units
// (e.g. stdev = 2^-15), mapped to the torus.
inline Torus gaussianTorus(double stdev) {
    std::normal_distribution<double> dist(0.0, stdev);
    return dtot(dist(rng()));
}

}  // namespace tfhe
