// test08_blindrotate.cpp -- Modulus switching, sample extraction and blind rotation.
#include <cstdio>
#include "../src/blindrotate.hpp"
#include "check.hpp"
using namespace tfhe;

int main() {
    setSeed(8);
    Params p;
    SecretKey sk = secretKeyGen(p);
    EvalKey ek = evalKeyGen(sk);
    const int N = p.N, twoN = 2 * N;
    // s' seen as an LWE key of dimension N: decrypts extracted samples
    LweKey extracted{sk.trlwe.s};
    const Torus mu = encode(1, 8);

    // 1) rounding to Z_2N, including the wrap-around of 0.9999 to 1 = 0
    CHECK(modSwitch(dtot(0.25), twoN) == twoN / 4);
    CHECK(modSwitch(dtot(0.5), twoN) == N);
    CHECK(modSwitch(dtot(0.9999), twoN) == 0);

    // 2) extraction turns the constant coefficient of a TRLWE into an LWE
    TorusPoly msg(N);
    for (auto& c : msg) c = encode(static_cast<int>(rng()() % 8), 8);
    TrlweSample d = trlweEncrypt(msg, p.bkStdev, sk.trlwe, ek.fft);
    LweSample e = sampleExtract(d);
    CHECK(static_cast<int>(e.a.size()) == N);
    CHECK(decode(lwePhase(e, extracted), 8) == decode(msg[0], 8));

    // 3) blind rotation of a trivial sample (0, b) with b = 3/2N (2N = 2^11):
    //    all abar_i are zero, so the result is exactly X^-3 * X^10 * mu = X^7 * mu
    TorusPoly v(N, 0);
    v[10] = mu;
    LweSample triv = lweTrivial(Torus(3) << 21, p.n);
    TorusPoly ph = trlwePhase(blindRotate(v, triv, ek), sk.trlwe, ek.fft);
    for (int j = 0; j < N; ++j) CHECK(ph[j] == (j == 7 ? mu : 0));

    // 4) gate-bootstrapping test polynomial v = mu (1 + X + ... + X^(N-1)):
    //    positive phase gives +1/8 (decodes to 1), negative gives -1/8 = 7/8
    TorusPoly ones(N, mu);
    for (int rep = 0; rep < 20; ++rep) {
        bool positive = rep % 2;
        LweSample c = lweEncrypt(positive ? mu : -mu, p.lweStdev, sk.lwe);
        LweSample r = sampleExtract(blindRotate(ones, c, ek));
        CHECK(decode(lwePhase(r, extracted), 8) == (positive ? 1 : 7));
    }

    // 5) programmable bootstrapping: coefficient j stores (4j/N)/8, so the
    //    input (2m + 1)/16, m in {0, 1, 2, 3}, comes out as m/8
    TorusPoly table(N);
    for (int j = 0; j < N; ++j) table[j] = encode(j * 4 / N, 8);
    for (int m = 0; m < 4; ++m) {
        LweSample c = lweEncrypt(encode(2 * m + 1, 16), p.lweStdev, sk.lwe);
        LweSample r = sampleExtract(blindRotate(table, c, ek));
        CHECK(decode(lwePhase(r, extracted), 8) == m);
    }

    std::printf("test08_blindrotate: OK\n");
    return 0;
}
