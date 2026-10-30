// test02_lwe.cpp -- LWE encryption, decryption and homomorphic linear operations.
#include <cmath>
#include <cstdio>
#include "../src/lwe.hpp"
#include "check.hpp"
using namespace tfhe;

int main() {
    setSeed(2);
    const int n = 630;
    const double stdev = std::ldexp(1.0, -15);
    LweKey key = lweKeyGen(n);

    // 1) encrypt and decrypt messages of Z_8
    for (int rep = 0; rep < 1000; ++rep) {
        int mu = rep % 8;
        LweSample c = lweEncrypt(encode(mu, 8), stdev, key);
        CHECK(decode(lwePhase(c, key), 8) == mu);
    }

    // 2) the sum of ciphertexts encrypts the sum of messages
    LweSample c1 = lweEncrypt(encode(3, 8), stdev, key);
    LweSample c2 = lweEncrypt(encode(6, 8), stdev, key);
    LweSample sum = c1;
    lweAddTo(sum, c2);
    CHECK(decode(lwePhase(sum, key), 8) == (3 + 6) % 8);

    // 3) subtraction and multiplication by an integer
    LweSample dif = c1;
    lweSubTo(dif, c2);
    CHECK(decode(lwePhase(dif, key), 8) == (3 - 6 + 8) % 8);
    LweSample tri = lweTrivial(0, n);
    lweAddMulTo(tri, 3, c1);
    CHECK(decode(lwePhase(tri, key), 8) == (3 * 3) % 8);

    // 4) trivial samples decrypt under any key and add public constants
    LweSample plus2 = lweTrivial(encode(2, 8), n);
    CHECK(decode(lwePhase(plus2, lweKeyGen(n)), 8) == 2);
    lweAddTo(plus2, c1);
    CHECK(decode(lwePhase(plus2, key), 8) == 3 + 2);

    // 5) fresh noise has the requested standard deviation, and the noise of a
    //    sum of two independent samples is sqrt(2) times larger
    const int reps = 20000;
    double fresh2 = 0, added2 = 0;
    for (int rep = 0; rep < reps; ++rep) {
        LweSample z1 = lweEncrypt(0, stdev, key), z2 = lweEncrypt(0, stdev, key);
        double e1 = ttod(lwePhase(z1, key));
        lweAddTo(z1, z2);
        double e12 = ttod(lwePhase(z1, key));
        fresh2 += e1 * e1;
        added2 += e12 * e12;
    }

    double sdFresh = std::sqrt(fresh2 / reps), sdAdded = std::sqrt(added2 / reps);
    std::printf("noise of a fresh sample: %.3e (requested %.3e)\n", sdFresh, stdev);
    std::printf("noise of the sum of two: %.3e (ratio %.3f; sqrt(2) = 1.414)\n",
                sdAdded, sdAdded / sdFresh);
    CHECK(std::fabs(sdFresh / stdev - 1.0) < 0.05);
    CHECK(std::fabs(sdAdded / sdFresh - std::sqrt(2.0)) < 0.05);

    std::printf("test02_lwe: OK\n");
    return 0;
}
