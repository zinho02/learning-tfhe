// keys.hpp -- Secret key (client side) and evaluation key (server side).
//
// The evaluation key has two parts:
//   * bootstrapping key: for each bit s_i of the LWE key, a TRGSW encryption
//     of s_i under the TRLWE key s' (n samples, stored in the FFT domain);
//   * key-switching key: LWE encryptions under s of v * s'_i / base^(j+1), for
//     every coefficient i of s', every digit level j and every digit value v.
// Each key encrypts the other one, so security relies on the usual circular
// security assumption of bootstrapped FHE schemes.
#pragma once
#include <vector>
#include "lwe.hpp"
#include "params.hpp"
#include "trgsw.hpp"

namespace tfhe {

// What the client keeps: s encrypts and decrypts the user's bits, and s' is
// used only inside bootstrapping.
struct SecretKey {
    Params p;
    LweKey lwe;
    TrlweKey trlwe;
};

// Table of N * t * base LWE samples stored in a single vector, with the
// entries of the same (i, j) contiguous in memory.
struct KeySwitchKey {
    int n, N, t, basebit, base;
    std::vector<LweSample> ks;

    // Encryption of v * s'_i / base^(j+1).
    const LweSample& at(int i, int j, int v) const {
        return ks[(i * t + j) * base + v];
    }
};

// Everything the server needs to evaluate gates; nothing here is secret. The
// Gadget and NegFFT objects are not keys, but depend only on the parameters
// and are bundled here so that every gate takes a single argument.
struct EvalKey {
    Params p;
    Gadget g;
    NegFFT fft;
    std::vector<TrgswFFT> bk;
    KeySwitchKey ksk;
    explicit EvalKey(const Params& p_) : p(p_), g(p_.ell, p_.Bgbit), fft(p_.N) {}
};

// Samples the two secret keys.
inline SecretKey secretKeyGen(const Params& p) {
    SecretKey sk;
    sk.p = p;
    sk.lwe = lweKeyGen(p.n);
    sk.trlwe = trlweKeyGen(p.N);
    return sk;
}

// Derives the evaluation key from the secret key. The push_back order of the
// key-switching table (i, then j, then v) must match KeySwitchKey::at. Entries
// with v = 0 are never used; a trivial encryption of 0 keeps indexing simple.
inline EvalKey evalKeyGen(const SecretKey& sk) {
    const Params& p = sk.p;
    EvalKey ek(p);

    ek.bk.reserve(p.n);
    for (int i = 0; i < p.n; ++i) {
        TrgswSample C = trgswEncrypt(sk.lwe.s[i], p.bkStdev, sk.trlwe, ek.g, ek.fft);
        ek.bk.push_back(toFFT(C, ek.fft));
    }

    KeySwitchKey& K = ek.ksk;
    K.n = p.n;
    K.N = p.N;
    K.t = p.ksT;
    K.basebit = p.ksBasebit;
    K.base = 1 << p.ksBasebit;
    K.ks.reserve(static_cast<size_t>(K.N) * K.t * K.base);
    for (int i = 0; i < K.N; ++i) {
        for (int j = 0; j < K.t; ++j) {
            Torus unit = 1u << (32 - (j + 1) * K.basebit);
            for (int v = 0; v < K.base; ++v) {
                Torus mu = v * sk.trlwe.s[i] * unit;
                K.ks.push_back(v == 0 ? lweTrivial(0, K.n)
                                    : lweEncrypt(mu, p.lweStdev, sk.lwe));
            }
        }
    }
    return ek;
}

}  // namespace tfhe
