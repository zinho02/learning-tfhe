// noise.cpp -- Output noise of the gate bootstrapping over several keys
// (Section 5 of the report).
//
// Usage: ./noise [keys] [gates per key] [prefix of the output files]
//
// For each key, evaluates the bootstrapping of random NAND inputs and records
// the error (phase minus expected message) after blind rotation + extraction
// and after key switching. Writes
//   <prefix>_keys.dat  per-key mean and variance of both errors
//   <prefix>_hist.dat  histogram (density) of the error after key switching,
//                      centered on the mean of its own key, all keys pooled
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "../../src/bootstrap.hpp"

using namespace tfhe;

int main(int argc, char** argv) {
    int keys = argc > 1 ? std::atoi(argv[1]) : 10;
    int gates = argc > 2 ? std::atoi(argv[2]) : 1000;
    std::string prefix = argc > 3 ? argv[3] : "noise";
    setSeed(2026);
    Params p;

    std::FILE* fk = std::fopen((prefix + "_keys.dat").c_str(), "w");
    std::fprintf(fk, "key meanBR varBR meanKS varKS errors\n");
    std::vector<double> centered;
    for (int k = 0; k < keys; ++k) {
        SecretKey sk = secretKeyGen(p);
        EvalKey ek = evalKeyGen(sk);
        LweKey extracted{sk.trlwe.s};
        std::vector<double> eBR(gates), eKS(gates);
        int errors = 0;
        for (int i = 0; i < gates; ++i) {
            bool b1 = uniformBit(), b2 = uniformBit();
            LweSample t = lweTrivial(MU, p.n);
            lweSubTo(t, encryptBit(b1, sk));
            lweSubTo(t, encryptBit(b2, sk));
            Torus expected = !(b1 && b2) ? MU : -MU;
            LweSample v = bootstrapNoKS(MU, t, ek);
            eBR[i] = ttod(lwePhase(v, extracted) - expected);
            LweSample w = keySwitch(v, ek.ksk);
            eKS[i] = ttod(lwePhase(w, sk.lwe) - expected);
            if (decryptBit(w, sk) != !(b1 && b2)) ++errors;
        }
        double mBR = 0, mKS = 0, vBR = 0, vKS = 0;
        for (int i = 0; i < gates; ++i) mBR += eBR[i] / gates, mKS += eKS[i] / gates;
        for (int i = 0; i < gates; ++i) {
            vBR += (eBR[i] - mBR) * (eBR[i] - mBR) / gates;
            vKS += (eKS[i] - mKS) * (eKS[i] - mKS) / gates;
            centered.push_back(eKS[i] - mKS);
        }
        std::fprintf(fk, "%d %.4e %.4e %.4e %.4e %d\n", k, mBR, vBR, mKS, vKS, errors);
        std::fflush(fk);
        std::fprintf(stderr, "key %d: BR mean %+.2e var %.2e | KS mean %+.2e var %.2e | errors %d\n",
                     k, mBR, vBR, mKS, vKS, errors);
    }
    std::fclose(fk);

    const int bins = 48;
    const double width = 0.024 / bins;
    std::vector<int> count(bins, 0);
    for (double e : centered) {
        int b = static_cast<int>(std::floor((e + 0.012) / width));
        if (b >= 0 && b < bins) ++count[b];
    }
    std::FILE* fh = std::fopen((prefix + "_hist.dat").c_str(), "w");
    std::fprintf(fh, "center density\n");
    for (int b = 0; b < bins; ++b)
        std::fprintf(fh, "%.6e %.6e\n", -0.012 + (b + 0.5) * width,
                     count[b] / (centered.size() * width));
    std::fclose(fh);
    return 0;
}
