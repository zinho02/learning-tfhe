// test07_keys.cpp -- The evaluation key contains what it should.
#include <cmath>
#include <cstdio>
#include "../src/keys.hpp"
#include "check.hpp"
using namespace tfhe;

int main() {
    setSeed(7);
    Params p;
    SecretKey sk = secretKeyGen(p);
    EvalKey ek = evalKeyGen(sk);

    // 1) sizes
    CHECK(static_cast<int>(sk.lwe.s.size()) == p.n);
    CHECK(static_cast<int>(sk.trlwe.s.size()) == p.N);
    CHECK(static_cast<int>(ek.bk.size()) == p.n);
    CHECK(static_cast<int>(ek.bk[0].a.size()) == 2 * p.ell);
    CHECK(static_cast<int>(ek.ksk.ks.size()) == p.N * p.ksT * (1 << p.ksBasebit));

    // 2) bk[i] encrypts s_i: used as a CMux selector it must pick d1 when
    //    s_i = 1 and d0 when s_i = 0
    TorusPoly mu0(p.N, encode(1, 8)), mu1(p.N, encode(5, 8));
    TrlweSample d0 = trlweEncrypt(mu0, p.bkStdev, sk.trlwe, ek.fft);
    TrlweSample d1 = trlweEncrypt(mu1, p.bkStdev, sk.trlwe, ek.fft);
    for (int i = 0; i < p.n; ++i) {
        TrlweSample r = cmux(ek.bk[i], d1, d0, ek.g, ek.fft);
        TorusPoly phase = trlwePhase(r, sk.trlwe, ek.fft);
        CHECK(decode(phase[0], 8) == (sk.lwe.s[i] ? 5 : 1));
    }

    // 3) entry (i, j, v) of the key-switching table encrypts, under s,
    //    v * s'_i / base^(j+1) with LWE-sized noise
    const KeySwitchKey& K = ek.ksk;
    double maxErr = 0;
    for (int i = 0; i < K.N; ++i) {
        for (int j = 0; j < K.t; ++j) {
            for (int v = 0; v < K.base; ++v) {
                Torus mu = v * sk.trlwe.s[i] * (1u << (32 - (j + 1) * K.basebit));
                double err = std::fabs(ttod(lwePhase(K.at(i, j, v), sk.lwe) - mu));
                if (err > maxErr) maxErr = err;
            }
        }
    }
    std::printf("largest noise in the key-switching table: %.2e (standard deviation %.2e)\n",
                maxErr, p.lweStdev);
    CHECK(maxErr < 8 * p.lweStdev);

    std::printf("test07_keys: OK\n");
    return 0;
}
