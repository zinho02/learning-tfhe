// bootstrap.hpp -- Key switching, gate bootstrapping and bit encoding.
//
// bootstrap(c) = KeySwitch(SampleExtract(BlindRotate(v, c))). The output noise
// depends only on the evaluation key, not on the noise of c: the input noise
// only has to be small enough for the phase to fall in the correct half of
// the torus. This is what allows circuits of arbitrary depth.
#pragma once
#include "blindrotate.hpp"

namespace tfhe {

// Switches an LWE sample from key s' (dimension N) to key s (dimension n),
// preserving its message. Each a_i is rounded to t digits in base 2^basebit
// (adding half of the least significant digit, as in the gadget
// decomposition), and for every nonzero digit v at level j the tabulated
// encryption of v * s'_i / base^(j+1) is subtracted from (0, b). There are no
// multiplications: at most N * t subtractions of samples of size n + 1.
inline LweSample keySwitch(const LweSample& c, const KeySwitchKey& K) {
    LweSample r = lweTrivial(c.b, K.n);

    const uint32_t round = 1u << (32 - (1 + K.basebit * K.t));
    const uint32_t mask = K.base - 1;
    for (int i = 0; i < K.N; ++i) {
        uint32_t ai = c.a[i] + round;
        for (int j = 0; j < K.t; ++j) {
            uint32_t v = (ai >> (32 - (j + 1) * K.basebit)) & mask;
            if (v != 0) lweSubTo(r, K.at(i, j, v));
        }
    }
    return r;
}

// Bootstrapping without the final key switching. With the test polynomial
// v = mu * (1 + X + ... + X^(N-1)) the output is an LWE sample of dimension N,
// under s', of +mu if the phase of c lies in (0, 1/2) and -mu if it lies in
// (-1/2, 0). Kept separate because the MUX gate adds two of these outputs
// before a single key switching.
inline LweSample bootstrapNoKS(Torus mu, const LweSample& c, const EvalKey& ek) {
    TorusPoly testv(ek.p.N, mu);
    return sampleExtract(blindRotate(testv, c, ek));
}

// Full bootstrapping: the output is again a sample of dimension n under s.
inline LweSample bootstrap(Torus mu, const LweSample& c, const EvalKey& ek) {
    return keySwitch(bootstrapNoKS(mu, c, ek), ek.ksk);
}

// Bit encoding: true is +1/8 and false is -1/8.
const Torus MU = encode(1, 8);

inline LweSample encryptBit(bool m, const SecretKey& sk) {
    return lweEncrypt(m ? MU : -MU, sk.p.lweStdev, sk.lwe);
}

// Decrypts by looking only at the sign of the phase, so the noise tolerance
// is 1/8 on each side (the distance from +-1/8 to 0 and to 1/2).
inline bool decryptBit(const LweSample& c, const SecretKey& sk) {
    return static_cast<int32_t>(lwePhase(c, sk.lwe)) > 0;
}

}  // namespace tfhe
