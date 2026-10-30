// test04_fft.cpp -- Negacyclic FFT against the schoolbook product: exactness and speed.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include "../src/fft.hpp"
#include "check.hpp"
using namespace tfhe;

int main() {
    setSeed(4);
    const int N = 1024;
    NegFFT fft(N);

    // 1) forward followed by backward returns the original polynomial exactly
    TorusPoly p(N), back;
    for (auto& c : p) c = uniformTorus();
    FreqPoly f;
    fft.forward(f, p);
    fft.backward(back, f);
    CHECK(back == p);

    // 2) FFT product against the schoolbook product, on inputs shaped like the
    //    scheme's: digits in [-64, 64) times arbitrary torus coefficients.
    //    The bootstrapping key noise (2^-25) is 128 units of 2^-32.
    int32_t maxErr = 0;
    IntPoly a(N);
    for (int rep = 0; rep < 20; ++rep) {
        for (auto& c : a) c = static_cast<int32_t>(rng()() % 128) - 64;
        for (auto& c : p) c = uniformTorus();
        TorusPoly r1 = mulNaive(a, p), r2 = mulFFT(fft, a, p);
        for (int j = 0; j < N; ++j) {
            int32_t err = std::abs(static_cast<int32_t>(r1[j] - r2[j]));
            if (err > maxErr) maxErr = err;
        }
    }
    std::printf("largest FFT error: %d units of 2^-32\n", maxErr);
    CHECK(maxErr < 16);

    // 3) timing of both methods; the volatile sink keeps the compiler from
    //    removing loops whose results are otherwise unused
    using clk = std::chrono::steady_clock;
    volatile Torus sink = 0;
    auto t0 = clk::now();
    for (int rep = 0; rep < 50; ++rep) sink = sink + mulNaive(a, p)[0];
    auto t1 = clk::now();
    for (int rep = 0; rep < 5000; ++rep) sink = sink + mulFFT(fft, a, p)[0];
    auto t2 = clk::now();
    double usNaive = std::chrono::duration<double, std::micro>(t1 - t0).count() / 50;
    double usFFT = std::chrono::duration<double, std::micro>(t2 - t1).count() / 5000;
    std::printf("schoolbook product: %.1f us FFT product: %.1f us (%.0fx)\n",
                usNaive, usFFT, usNaive / usFFT);
    std::printf("test04_fft: OK\n");
    return 0;
}
