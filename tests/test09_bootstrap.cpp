// test09_bootstrap.cpp -- Key switching and the complete gate bootstrapping.
#include <cmath>
#include <cstdio>
#include "../src/bootstrap.hpp"
#include "check.hpp"
using namespace tfhe;

int main() {
    setSeed(9);
    Params p;
    SecretKey sk = secretKeyGen(p);
    EvalKey ek = evalKeyGen(sk);
    LweKey extracted{sk.trlwe.s};

    // 1) encrypt and decrypt bits
    for (int rep = 0; rep < 100; ++rep) {
        bool bit = rep % 2;
        CHECK(decryptBit(encryptBit(bit, sk), sk) == bit);
    }

    // 2) key switching from s' (dimension N) to s (dimension n) keeps the message
    for (int mu = 0; mu < 8; ++mu) {
        LweSample big = lweEncrypt(encode(mu, 8), p.bkStdev, extracted);
        LweSample small = keySwitch(big, ek.ksk);
        CHECK(static_cast<int>(big.a.size()) == p.N);
        CHECK(static_cast<int>(small.a.size()) == p.n);
        CHECK(decode(lwePhase(small, sk.lwe), 8) == mu);
    }

    // 3) bootstrapping without key switching: output +-1/8 under s'
    for (int rep = 0; rep < 20; ++rep) {
        bool bit = rep % 2;
        LweSample r = bootstrapNoKS(MU, encryptBit(bit, sk), ek);
        CHECK(decode(lwePhase(r, extracted), 8) == (bit ? 1 : 7));
    }

    // 4) full bootstrapping keeps the bit and can be chained
    for (int rep = 0; rep < 20; ++rep) {
        bool bit = rep % 2;
        LweSample c = encryptBit(bit, sk);
        for (int k = 0; k < 3; ++k) c = bootstrap(MU, c, ek);
        CHECK(decryptBit(c, sk) == bit);
    }

    // 5) bootstrapping removes noise: inputs with standard deviation 0.01
    //    (more than 300 times the usual) come out with the usual output noise
    const double noisy = 0.01;
    const int reps = 50;
    double in2 = 0, out2 = 0;
    for (int rep = 0; rep < reps; ++rep) {
        LweSample c = lweEncrypt(MU, noisy, sk.lwe);
        LweSample r = bootstrap(MU, c, ek);
        CHECK(decryptBit(r, sk) == true);
        double eIn = ttod(lwePhase(c, sk.lwe) - MU);
        double eOut = ttod(lwePhase(r, sk.lwe) - MU);
        in2 += eIn * eIn;
        out2 += eOut * eOut;
    }
    double sdIn = std::sqrt(in2 / reps), sdOut = std::sqrt(out2 / reps);
    std::printf("noise before bootstrapping: %.2e after: %.2e\n", sdIn, sdOut);
    CHECK(sdOut < sdIn / 2);

    std::printf("test09_bootstrap: OK\n");
    return 0;
}
