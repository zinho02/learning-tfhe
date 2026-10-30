// lwe.hpp -- LWE samples over the torus.
//
// With a binary secret key s in B^n, an encryption of mu in T is the pair
//     c = (a, b),   a uniform in T^n,   b = <a, s> + mu + e,
// where e is a small Gaussian noise. The phase phi(c) = b - <a, s> = mu + e is
// linear in c, which makes the scheme additively homomorphic.
#pragma once
#include <vector>
#include "torus.hpp"

namespace tfhe {

// Secret key: n bits, stored as int32_t so that a[i] * s[i] is a plain
// "integer times torus" product.
struct LweKey {
    std::vector<int32_t> s;
};

// LWE ciphertext: mask a (n torus elements) and body b.
struct LweSample {
    std::vector<Torus> a;
    Torus b = 0;
    explicit LweSample(int n = 0) : a(n,0) {}
};

// Samples a uniform binary key of dimension n.
inline LweKey lweKeyGen(int n) {
    LweKey key;
    key.s.resize(n);
    for (auto& si : key.s) si = uniformBit();
    return key;
}

// Encrypts the torus point mu with Gaussian noise of standard deviation stdev.
// The body starts as mu + e and accumulates <a, s> while the mask is sampled.
inline LweSample lweEncrypt(Torus mu, double stdev, const LweKey& key) {
    int n = static_cast<int>(key.s.size());
    LweSample c(n);
    c.b = mu + gaussianTorus(stdev);
    for (int i = 0; i < n; ++i) {
        c.a[i] = uniformTorus();
        c.b += c.a[i] * key.s[i];
    }
    return c;
}

// Phase b - <a, s> = mu + e, i.e. decryption without the final rounding.
// Callers either decode it or subtract the expected message to measure noise.
inline Torus lwePhase(const LweSample& c, const LweKey& key) {
    Torus phase = c.b;
    for (size_t i = 0; i < key.s.size(); ++i) phase -= c.a[i] * key.s[i];
    return phase;
}

// Trivial (noiseless) encryption (0, mu). Its phase is mu under any key, so it
// is how public constants are added to ciphertexts.
inline LweSample lweTrivial(Torus mu, int n) {
    LweSample c(n);
    c.b = mu;
    return c;
}

// r <- r + k * c for a small integer k. Messages add up (scaled by k) and so
// do the noises; the noise standard deviation of c is scaled by |k|.
inline void lweAddMulTo(LweSample& r, int32_t k, const LweSample& c) {
    for (size_t i = 0; i < r.a.size(); ++i) r.a[i] += k * c.a[i];
    r.b += k * c.b;
}

// r <- r + c and r <- r - c.
inline void lweAddTo(LweSample& r, const LweSample& c) { lweAddMulTo(r, 1, c); }
inline void lweSubTo(LweSample& r, const LweSample& c) { lweAddMulTo(r, -1, c); }

}  // namespace tfhe
