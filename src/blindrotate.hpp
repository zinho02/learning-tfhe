// blindrotate.hpp -- Blind rotation and sample extraction.
//
// Given an LWE sample c = (a, b) under s, blind rotation produces a TRLWE
// encryption of X^(-phi) * v, where phi = bbar - sum_i abar_i s_i (mod 2N) is
// the phase of c rounded to Z_2N and v is a public test polynomial. Since
//     X^(-phi) = X^(-bbar) * prod_i X^(abar_i s_i),
// each factor is a choice controlled by the secret bit s_i, i.e. a CMux with
// the bootstrapping key bk[i]. The constant coefficient of X^(-phi) * v is
// v_phi for 0 <= phi < N and -v_(phi-N) for N <= phi < 2N, so v acts as a
// lookup table indexed by the encrypted phase.
#pragma once
#include "keys.hpp"

namespace tfhe {

// Rounds a torus element to the nearest multiple of 1/(2N) and returns the
// numerator in {0, ..., 2N-1}: (t / 2^32) * 2N computed in 64 bits, plus 2^31
// for rounding, then the 32 low bits are dropped. "% twoN" maps 2N to 0.
inline int modSwitch(Torus t, int twoN) {
    return static_cast<int>(((uint64_t{t} * twoN + (uint64_t{1} << 31)) >> 32) % twoN);
}

// Starts from the trivial sample (0, X^(-bbar) * v) and applies, for each i,
//     acc <- CMux(bk[i], X^abar_i * acc, acc) = acc + bk[i] [.] ((X^abar_i - 1) * acc).
// X^(-bbar) is written as X^(2N - bbar) because mulByXe takes exponents in
// [0, 2N). Steps with abar_i = 0 are skipped, since X^0 = 1. The loop runs n
// times and each step costs one external product: this is where almost all
// of the gate time is spent.
inline TrlweSample blindRotate(const TorusPoly& testv, const LweSample& c, const EvalKey& ek) {
    int N = ek.p.N, twoN = 2 * N;
    int bbar = modSwitch(c.b, twoN);

    TrlweSample acc(N);
    acc.b = mulByXe(testv, (twoN - bbar) % twoN);

    for (int i = 0; i < ek.p.n; ++i) {
        int abar = modSwitch(c.a[i], twoN);
        if (abar == 0) continue;

        TrlweSample diff(N);
        diff.a = mulByXe(acc.a, abar);
        diff.b = mulByXe(acc.b, abar);
        polySubTo(diff.a, acc.a);
        polySubTo(diff.b, acc.b);
        TrlweSample prod = externalProduct(ek.bk[i], diff, ek.g, ek.fft);
        polyAddTo(acc.a, prod.a);
        polyAddTo(acc.b, prod.b);
    }
    return acc;
}

// Extracts an LWE sample of the constant coefficient of a TRLWE message. The
// constant coefficient of b - a*s' is b_0 - a_0 s'_0 - sum_{j>=1} a_(N-j) s'_j,
// so the result is (a', b_0) with a' = (a_0, -a_(N-1), ..., -a_1), of
// dimension N, decrypted by the LWE key formed by the coefficients of s'.
// It is a rearrangement of coefficients and adds no noise.
inline LweSample sampleExtract(const TrlweSample& c) {
    int N = static_cast<int>(c.a.size());
    LweSample r(N);
    r.a[0] = c.a[0];
    for (int j = 1; j < N; ++j) r.a[j] = -c.a[N - j];
    r.b = c.b[0];
    return r;
}

}  // namespace tfhe
