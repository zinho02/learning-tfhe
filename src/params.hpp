// params.hpp -- Parameter set of the scheme.
//
// The default values reproduce the ~128-bit security parameter set of the
// original TFHE library (default_128bit_gate_bootstrapping_parameters), so
// that the comparison in compare/bench_original.cpp uses identical settings.
#pragma once
#include <cmath>

namespace tfhe {

struct Params {
    int n = 630;                             // dimension of the LWE key s
    int N = 1024;                            // ring degree of Z[X]/(X^N + 1); power of two
    int ell = 3;                             // number of levels of the gadget decomposition
    int Bgbit = 7;                           // log2 of the gadget base Bg
    int ksT = 8;                             // number of digits used by key switching
    int ksBasebit = 2;                       // log2 of the key-switching base
    double lweStdev = std::ldexp(1.0, -15);  // noise of LWE samples and of the key-switching key
    double bkStdev = std::ldexp(1.0, -25);   // noise of the TRGSW samples of the bootstrapping key
};

}  // namespace tfhe
