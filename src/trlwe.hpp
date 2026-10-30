// trlwe.hpp -- TRLWE samples (ring LWE over the torus), with k = 1.
//
// With a binary polynomial key s' in B_N[X], an encryption of mu in T_N[X] is
//     c = (a, b),   b = a * s' + mu + e   (mod X^N + 1),
// and its phase b - a * s' = mu + e is a polynomial. Besides being additive,
// a TRLWE sample can be rotated: (X^e a, X^e b) encrypts X^e * mu with the
// noise only rotated, never amplified. Blind rotation is built on this.
#pragma once
#include "fft.hpp"

namespace tfhe {

// Secret key: a polynomial with N binary coefficients.
struct TrlweKey {
    IntPoly s;
};

// TRLWE ciphertext: the pair of polynomials (a, b).
struct TrlweSample {
    TorusPoly a, b;
    explicit TrlweSample(int N = 0) : a(N, 0), b(N, 0) {}
};

inline TrlweKey trlweKeyGen(int N) {
    TrlweKey key;
    key.s.resize(N);
    for (auto& sj : key.s) sj = uniformBit();
    return key;
}

// Encrypts the polynomial mu, with independent Gaussian noise of standard
// deviation stdev on every coefficient.
inline TrlweSample trlweEncrypt(const TorusPoly& mu, double stdev,
                                const TrlweKey& key, const NegFFT& fft) {
    int N = fft.N();
    TrlweSample c(N);
    for (auto& aj : c.a) aj = uniformTorus();
    c.b = mulFFT(fft, key.s, c.a);
    for (int j = 0; j < N; ++j) c.b[j] += mu[j] + gaussianTorus(stdev);
    return c;
}

// Phase b - a * s' = mu + e, without rounding.
inline TorusPoly trlwePhase(const TrlweSample& c, const TrlweKey& key,
                            const NegFFT& fft) {
    TorusPoly phase = c.b;
    polySubTo(phase, mulFFT(fft, key.s, c.a));
    return phase;
}

}  // namespace tfhe
