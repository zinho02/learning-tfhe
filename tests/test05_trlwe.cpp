// test05_trlwe.cpp -- TRLWE encryption, addition, rotation by X^e and trivial samples.
#include <cstdio>
#include <vector>
#include "../src/trlwe.hpp"
#include "check.hpp"
using namespace tfhe;

int main() {
    setSeed(5);
    const int N = 1024;
    const double stdev = std::ldexp(1.0, -25);
    NegFFT fft(N);
    TrlweKey key = trlweKeyGen(N);

    // a message polynomial with coefficients in Z_8
    std::vector<int> msg(N);
    TorusPoly mu(N);
    for (int j = 0; j < N; ++j) {
        msg[j] = static_cast<int>(rng()() % 8);
        mu[j] = encode(msg[j], 8);
    }

    // 1) every coefficient of the phase decodes to the message
    TrlweSample c = trlweEncrypt(mu, stdev, key, fft);
    TorusPoly phase = trlwePhase(c, key, fft);
    for (int j = 0; j < N; ++j) CHECK(decode(phase[j], 8) == msg[j]);

    // 2) the sum of ciphertexts encrypts the sum of messages
    TrlweSample c2 = trlweEncrypt(mu, stdev, key, fft);
    polyAddTo(c2.a, c.a);
    polyAddTo(c2.b, c.b);
    phase = trlwePhase(c2, key, fft);
    for (int j = 0; j < N; ++j) CHECK(decode(phase[j], 8) == (2 * msg[j]) % 8);

    // 3) multiplying both components by X^e rotates the message
    const int e = 100;
    TrlweSample rot(N);
    rot.a = mulByXe(c.a, e);
    rot.b = mulByXe(c.b, e);
    TorusPoly rotMu = mulByXe(mu, e);
    phase = trlwePhase(rot, key, fft);
    for (int j = 0; j < N; ++j) CHECK(decode(phase[j], 8) == decode(rotMu[j], 8));

    // 4) the trivial sample (0, mu) has phase exactly mu
    TrlweSample triv(N);
    triv.b = mu;
    CHECK(trlwePhase(triv, key, fft) == mu);

    std::printf("test05_trlwe: OK\n");
    return 0;
}
