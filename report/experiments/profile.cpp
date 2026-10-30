// profile.cpp -- Cost breakdown of one gate bootstrapping (Section 5 of the report).
//
// Usage: ./profile [repetitions of each primitive]
//
// Times every primitive that blind rotation calls, multiplies each one by the
// number of times it runs per gate, and writes the breakdown to stdout as
// "name microseconds" lines (the format of ../data/profile.dat).
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include "../../src/gates.hpp"

using namespace tfhe;
using clk = std::chrono::steady_clock;

// Average time of f() in microseconds.
template <class F>
static double usPer(int reps, F f) {
    auto t0 = clk::now();
    for (int i = 0; i < reps; ++i) f();
    return std::chrono::duration<double, std::micro>(clk::now() - t0).count() / reps;
}

int main(int argc, char** argv) {
    int reps = argc > 1 ? std::atoi(argv[1]) : 20000;
    setSeed(42);
    Params p;
    SecretKey sk = secretKeyGen(p);
    EvalKey ek = evalKeyGen(sk);
    const int N = p.N;

    TorusPoly pa(N), pb(N), out;
    IntPoly digits(N);
    for (auto& c : pa) c = uniformTorus();
    for (auto& c : pb) c = uniformTorus();
    for (auto& c : digits) c = static_cast<int32_t>(rng()() % 128) - 64;
    FreqPoly f, g, scratch, acc(N / 2);
    ek.fft.forward(g, pb);
    std::vector<IntPoly> dec;
    TrlweSample t(N), d(N);
    t.a = pa;
    t.b = pb;
    d.a = pb;
    d.b = pa;

    double fwd = usPer(reps, [&] { ek.fft.forward(f, digits); });
    double copy = usPer(reps, [&] { scratch = g; });
    double bwd = usPer(reps, [&] { scratch = g; ek.fft.backward(out, scratch); }) - copy;
    double mac = usPer(reps, [&] { mulAcc(acc, f, g); });
    // The same product, but reading each time a row of a different TRGSW
    // sample of the bootstrapping key, as blind rotation does: the 59 MiB key
    // does not fit in cache and is streamed from memory at every gate.
    int row = 0, sample = 0;
    double macKey = usPer(reps, [&] {
        mulAcc(acc, f, ek.bk[sample].a[row]);
        if (++row == 2 * p.ell) row = 0, sample = (sample + 1) % p.n;
    });
    double dcp = usPer(reps, [&] { decompose(dec, pa, ek.g); });
    double rot = usPer(reps, [&] {
        TrlweSample diff(N);
        diff.a = mulByXe(t.a, 77);
        diff.b = mulByXe(t.b, 77);
        polySubTo(diff.a, t.a);
        polySubTo(diff.b, t.b);
        polyAddTo(d.a, diff.a);
        polyAddTo(d.b, diff.b);
    });
    TrlweSample prod;
    double ext = usPer(reps / 10, [&] { prod = externalProduct(ek.bk[0], t, ek.g, ek.fft); });
    sample = 0;
    double extKey = usPer(reps / 10, [&] {
        prod = externalProduct(ek.bk[sample], t, ek.g, ek.fft);
        sample = (sample + 1) % p.n;
    });

    LweSample c = encryptBit(true, sk), y = encryptBit(false, sk), u, r;
    int steps = 0;
    for (int i = 0; i < p.n; ++i) steps += modSwitch(c.a[i], 2 * N) != 0;
    TorusPoly testv(N, MU);
    TrlweSample accBR;
    int brReps = reps / 400 > 10 ? reps / 400 : 10;
    double br = usPer(brReps, [&] { accBR = blindRotate(testv, c, ek); });
    u = sampleExtract(accBR);
    double ks = usPer(brReps * 4, [&] { r = keySwitch(u, ek.ksk); });
    double gate = usPer(brReps, [&] { r = gateNAND(c, y, ek); });

    double fwdGate = steps * 2.0 * p.ell * fwd;
    double bwdGate = steps * 2.0 * bwd;
    double macGate = steps * 4.0 * p.ell * macKey;
    double decGate = steps * 2.0 * dcp;
    double rotGate = steps * rot;
    double other = gate - (fwdGate + bwdGate + macGate + decGate + rotGate + ks);

    std::fprintf(stderr, "primitives (us): forward %.3f backward %.3f mulAcc %.3f "
                 "(key from memory %.3f) decompose %.3f rotate %.3f\n",
                 fwd, bwd, mac, macKey, dcp, rot);
    std::fprintf(stderr, "externalProduct (us): same key sample %.2f, cycling over the key %.2f\n",
                 ext, extKey);
    std::fprintf(stderr, "steps %d blindRotate %.1f us keySwitch %.1f us NAND %.1f us\n",
                 steps, br, ks, gate);
    std::fprintf(stderr, "checksum %u\n", out[0] + prod.b[1] + r.b + accBR.b[2] +
                 static_cast<Torus>(acc[3].real()) + dec[0][5] + d.a[7]);

    std::printf("part us\n");
    std::printf("fft-forward %.1f\n", fwdGate);
    std::printf("fft-backward %.1f\n", bwdGate);
    std::printf("pointwise %.1f\n", macGate);
    std::printf("decompose %.1f\n", decGate);
    std::printf("rotate %.1f\n", rotGate);
    std::printf("keyswitch %.1f\n", ks);
    std::printf("other %.1f\n", other);
    std::printf("total %.1f\n", gate);
    return 0;
}
