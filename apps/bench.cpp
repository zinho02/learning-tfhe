// bench.cpp -- Measurements of time, memory, noise and decryption error rate.
//
// Usage: ./build/bench [number of gates per measurement]
//
// The exit status is the number of decryption errors (0 means every gate was
// correct).
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "../src/gates.hpp"

using namespace tfhe;

using clk = std::chrono::steady_clock;

// Milliseconds between two instants.
static double ms(clk::time_point a, clk::time_point b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
}

// Theoretical upper bound on the noise variance after blind rotation (k = 1):
// n external products, each adding 2 ell N beta^2 var_BK of fresh noise plus
// (1 + N) eps^2 of decomposition error, with beta = Bg/2 and eps = 1/(2 Bg^ell).
static double varBlindRotate(const Params& p) {
    double beta = std::ldexp(1.0, p.Bgbit - 1);
    double eps = std::ldexp(1.0, -(p.ell * p.Bgbit + 1));
    return p.n * 2.0 * p.ell * p.N * beta * beta * p.bkStdev * p.bkStdev +
            p.n * (1.0 + p.N) * eps * eps;
}

// Theoretical upper bound on the variance added by key switching: up to t N
// table samples, plus the error of rounding each a_i to t digits.
static double varKeySwitch(const Params& p) {
    return 1.0 * p.ksT * p.N * p.lweStdev * p.lweStdev +
            p.N * std::ldexp(1.0, -2 * (p.ksBasebit * p.ksT + 1));
}

int main(int argc, char** argv) {
    int reps = argc > 1 ? std::atoi(argv[1]) : 200;
    Params p;

    // key generation time
    auto t0 = clk::now();
    SecretKey sk = secretKeyGen(p);
    auto t1 = clk::now();
    EvalKey ek = evalKeyGen(sk);
    auto t2 = clk::now();
    std::printf("secret key: %.2f ms evaluation key: %.2f ms\n",
                ms(t0, t1), ms(t1, t2));

    // memory, computed from the parameters: n TRGSW samples of 2 ell rows with
    // 2 polynomials of N/2 complex numbers; N t base LWE samples of n + 1
    // torus elements
    double mb = 1024.0 * 1024.0;
    double bkBytes = 1.0 * p.n * 2 * p.ell * 2 * (p.N / 2) * sizeof(cplx);
    double kskBytes = 1.0 * p.N * p.ksT * (1 << p.ksBasebit) * (p.n + 1) * sizeof(Torus);
    std::printf("bootstrapping key: %.1f MB key-switching key: %.1f MB\n",
                bkBytes / mb, kskBytes / mb);
    std::printf("ciphertext of one bit: %zu bytes\n", (p.n + 1) * sizeof(Torus));

    // time per gate, and the two halves of the bootstrapping separately
    LweSample x = encryptBit(1, sk), y = encryptBit(0, sk), r, u;
    auto a0 = clk::now();
    for (int i = 0; i < reps; ++i) r = gateNAND(x, y, ek);
    auto a1 = clk::now();
    std::printf("NAND: %.2f ms per gate (average of %d)\n", ms(a0, a1) / reps, reps);

    a0 = clk::now();
    for (int i = 0; i < reps; ++i) u = bootstrapNoKS(MU, x, ek);
    a1 = clk::now();
    for (int i = 0; i < reps; ++i) r = keySwitch(u, ek.ksk);
    auto a2 = clk::now();
    std::printf(" blind rotation + extraction: %.2f ms key switching: %.2f ms\n",
                ms(a0, a1) / reps, ms(a1, a2) / reps);

    // noise (phase minus expected message) before and after key switching, on
    // NAND inputs built by hand from random bits, and decryption error rate
    LweKey extracted{sk.trlwe.s};
    double sBR = 0, sBR2 = 0, sKS = 0, sKS2 = 0;
    int errors = 0;
    for (int i = 0; i < reps; ++i) {
        bool b1 = uniformBit(), b2 = uniformBit();
        LweSample t = lweTrivial(MU, p.n);
        lweSubTo(t, encryptBit(b1, sk));
        lweSubTo(t, encryptBit(b2, sk));
        Torus expected = !(b1 && b2) ? MU : -MU;

        LweSample v = bootstrapNoKS(MU, t, ek);
        double eBR = ttod(lwePhase(v, extracted) - expected);
        LweSample w = keySwitch(v, ek.ksk);
        double eKS = ttod(lwePhase(w, sk.lwe) - expected);
        sBR += eBR;
        sBR2 += eBR * eBR;
        sKS += eKS;
        sKS2 += eKS * eKS;
        if (decryptBit(w, sk) != !(b1 && b2)) ++errors;
    }
    double mBR = sBR / reps, mKS = sKS / reps;
    double vBR = varBlindRotate(p), vKS = varKeySwitch(p);
    std::printf("error after blind rotation: mean %+.2e variance %.2e (bound %.2e)\n",
                mBR, sBR2 / reps - mBR * mBR, vBR);
    std::printf("error after key switching: mean %+.2e variance %.2e (bound %.2e)\n",
                mKS, sKS2 / reps - mKS * mKS, vBR + vKS);
    std::printf("decryption errors: %d out of %d gates\n", errors, reps);

    return errors;
}
