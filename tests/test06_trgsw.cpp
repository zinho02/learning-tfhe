// test06_trgsw.cpp -- Gadget decomposition, external product and CMux.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "../src/trgsw.hpp"
#include "check.hpp"
using namespace tfhe;

int main() {
    setSeed(6);
    const int N = 1024;
    const double stdev = std::ldexp(1.0, -25);
    NegFFT fft(N);
    Gadget g(3, 7);
    TrlweKey key = trlweKeyGen(N);

    // 1) digits lie in [-Bg/2, Bg/2) and sum_i d_i / Bg^(i+1) reconstructs each
    //    coefficient up to 1 / (2 Bg^ell) = 2^-22, i.e. 2^10 units of 2^-32
    TorusPoly p(N);
    for (auto& c : p) c = uniformTorus();
    std::vector<IntPoly> dec;
    decompose(dec, p, g);
    for (int j = 0; j < N; ++j) {
        Torus rec = 0;
        for (int i = 0; i < g.ell; ++i) {
            CHECK(dec[i][j] >= -64 && dec[i][j] < 64);
            rec += dec[i][j] * g.h(i);
        }
        int32_t err = static_cast<int32_t>(p[j] - rec);
        CHECK(std::abs(err) <= (1 << (32 - g.ell * g.Bgbit - 1)));
    }

    // a TRLWE encryption of a polynomial with coefficients in Z_8
    std::vector<int> msg(N);
    TorusPoly mu(N);
    for (int j = 0; j < N; ++j) {
        msg[j] = static_cast<int>(rng()() % 8);
        mu[j] = encode(msg[j], 8);
    }
    TrlweSample d = trlweEncrypt(mu, stdev, key, fft);

    // 2) TRGSW(m) [.] TRLWE(mu) = TRLWE(m * mu), and the output noise
    for (int m = 0; m <= 1; ++m) {
        TrgswFFT C = toFFT(trgswEncrypt(m, stdev, key, g, fft), fft);
        TorusPoly phase = trlwePhase(externalProduct(C, d, g, fft), key, fft);
        double err2 = 0;
        for (int j = 0; j < N; ++j) {
            CHECK(decode(phase[j], 8) == m * msg[j]);
            double e = ttod(phase[j] - m * mu[j]);
            err2 += e * e;
        }
        std::printf("external product with m = %d: output noise %.2e\n",
                    m, std::sqrt(err2 / N));
    }

    // 3) CMux selects d1 when the bit is 1 and d0 when it is 0
    TorusPoly mu0(N, encode(1, 8)), mu1(N, encode(5, 8));
    TrlweSample d0 = trlweEncrypt(mu0, stdev, key, fft);
    TrlweSample d1 = trlweEncrypt(mu1, stdev, key, fft);
    for (int bit = 0; bit <= 1; ++bit) {
        TrgswFFT C = toFFT(trgswEncrypt(bit, stdev, key, g, fft), fft);
        TorusPoly phase = trlwePhase(cmux(C, d1, d0, g, fft), key, fft);
        for (int j = 0; j < N; ++j) CHECK(decode(phase[j], 8) == (bit ? 5 : 1));
    }

    std::printf("test06_trgsw: OK\n");
    return 0;
}
