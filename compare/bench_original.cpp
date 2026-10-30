// bench_original.cpp -- The measurements of apps/bench.cpp, repeated with the
// original TFHE library (https://github.com/tfhe/tfhe).
//
// Usage: ./build/bench_original [number of gates per measurement]
// Build: make compare TFHE_PREFIX=<install dir> [TFHE_FFT=spqlios-fma|nayuki-portable]
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <tfhe/tfhe.h>

using clk = std::chrono::steady_clock;

// Milliseconds between two instants.
static double ms(clk::time_point a, clk::time_point b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
}

int main(int argc, char** argv) {
    int reps = argc > 1 ? std::atoi(argv[1]) : 200;

    // The library starts from a fixed seed; without reseeding, every run would
    // generate the same keys and the same noise.
    std::random_device rd;
    uint32_t seed[] = {rd(), rd(), rd(), rd()};
    tfhe_random_generator_setSeed(seed, 4);
    std::srand(rd());

    // 128 selects n = 630, N = 1024, ell = 3, Bgbit = 7: the same as Params.
    TFheGateBootstrappingParameterSet* params =
        new_default_gate_bootstrapping_parameters(128);

    // Generates the secret key together with the evaluation ("cloud") key.
    auto t0 = clk::now();
    TFheGateBootstrappingSecretKeySet* key =
        new_random_gate_bootstrapping_secret_keyset(params);
    auto t1 = clk::now();
    const TFheGateBootstrappingCloudKeySet* cloud = &key->cloud;
    std::printf("key generation: %.2f ms\n", ms(t0, t1));

    LweSample* x = new_gate_bootstrapping_ciphertext(params);
    LweSample* y = new_gate_bootstrapping_ciphertext(params);
    LweSample* r = new_gate_bootstrapping_ciphertext(params);

    // time per gate
    bootsSymEncrypt(x, 1, key);
    bootsSymEncrypt(y, 0, key);
    auto a0 = clk::now();
    for (int i = 0; i < reps; ++i) bootsNAND(r, x, y, cloud);
    auto a1 = clk::now();
    std::printf("NAND: %.2f ms per gate (average of %d)\n", ms(a0, a1) / reps, reps);

    // output noise (phase minus +-1/8) and decryption error rate; Torus32 is
    // signed in the library, so the difference already carries its sign
    const Torus32 MU = modSwitchToTorus32(1, 8);
    double sum = 0, sum2 = 0;
    int errors = 0;
    for (int i = 0; i < reps; ++i) {
        int b1 = std::rand() % 2, b2 = std::rand() % 2, expected = !(b1 && b2);
        bootsSymEncrypt(x, b1, key);
        bootsSymEncrypt(y, b2, key);
        bootsNAND(r, x, y, cloud);
        Torus32 phase = lwePhase(r, key->lwe_key);
        double e = (phase - (expected ? MU : -MU)) / 4294967296.0;
        sum += e;
        sum2 += e * e;
        if (bootsSymDecrypt(r, key) != expected) ++errors;
    }
    double mean = sum / reps;
    std::printf("output error: mean %+.2e variance %.2e\n",
                mean, sum2 / reps - mean * mean);
    std::printf("decryption errors: %d out of %d gates\n", errors, reps);

    delete_gate_bootstrapping_ciphertext(r);
    delete_gate_bootstrapping_ciphertext(y);
    delete_gate_bootstrapping_ciphertext(x);
    delete_gate_bootstrapping_secret_keyset(key);
    delete_gate_bootstrapping_parameters(params);

    return errors;
}
