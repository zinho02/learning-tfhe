# learning-tfhe

A from-scratch, header-only C++17 implementation of the TFHE fully homomorphic
encryption scheme (Chillotti, Gama, Georgieva and Izabachène, *Journal of
Cryptology* 33(1), 2020), written for the course MO422 (Cryptographic
Algorithms) at Unicamp. It implements gate bootstrapping over the 32-bit
discretized torus with the 128-bit parameter set of the original
[TFHE library](https://github.com/tfhe/tfhe), and compares the result against
that library.

The library has no external dependencies: every module is a header in `src/`
with `inline` functions, so each test and application is a single `.cpp` file.

## Layout

```
Makefile
src/            the library, one header per layer (each includes the previous ones)
  params.hpp        parameter set (n = 630, N = 1024, ell = 3, Bg = 2^7, ...)
  torus.hpp         torus T = R/Z as uint32_t, message encoding, randomness
  lwe.hpp           LWE samples: key generation, encryption, phase, linear operations
  poly.hpp          polynomials mod X^N + 1: schoolbook product, rotation by X^e
  fft.hpp           negacyclic FFT (fold + twist + size-N/2 complex FFT)
  trlwe.hpp         TRLWE samples (ring LWE over the torus, k = 1)
  trgsw.hpp         gadget decomposition, TRGSW samples, external product, CMux
  keys.hpp          secret key and evaluation key (bootstrapping + key-switching keys)
  blindrotate.hpp   modulus switching, blind rotation, sample extraction
  bootstrap.hpp     key switching, gate bootstrapping, bit encoding (+-1/8)
  gates.hpp         NAND, AND, OR, NOR, XOR, XNOR, NOT, MUX
tests/          one unit test per module (test01_torus ... test10_gates) and the CHECK macro
apps/
  adder.cpp         homomorphic n-bit ripple-carry adder
  bench.cpp         time, memory, noise and error-rate measurements
compare/
  bench_original.cpp  the same measurements using the original TFHE library
report/         final report (LaTeX, SBC template)
  main.tex, figures/  text and TikZ figures; build with `make` inside report/
  experiments/        profile.cpp, noise.cpp, run_all.sh and summarize.py
  data/               raw measurements (data/raw/) and the numbers derived from them
slides/         presentation (Beamer + TikZ, ~30 min); build with `make` inside slides/
```

The data flow of one gate is

```
inputs (LWE, key s) -> linear combination -> blind rotation -> sample extraction
                    -> key switching -> output (LWE, key s)
```

where blind rotation uses the bootstrapping key (n TRGSW samples) and key
switching uses the key-switching key (a table of LWE samples).

## Building and running

Requires a C++17 compiler (g++ 9 or newer, or clang++) and `make`. Always
build with optimizations (the default `CXXFLAGS` use `-O3 -march=native`).

```sh
make test                   # build and run the ten unit tests, in order
make                        # also build apps/adder and apps/bench
./build/adder 8 173 94      # homomorphic 8-bit addition of 173 and 94
./build/bench 500           # timing, key sizes, noise and error rate over 500 gates
```

### Comparison with the original library

Install the original library (both FFT back ends) under a prefix, e.g.:

```sh
git clone https://github.com/tfhe/tfhe.git && cd tfhe && mkdir build && cd build
cmake ../src -DCMAKE_BUILD_TYPE=optim -DENABLE_TESTS=off \
      -DENABLE_SPQLIOS_FMA=on -DENABLE_NAYUKI_PORTABLE=on \
      -DENABLE_SPQLIOS_AVX=off -DENABLE_NAYUKI_AVX=off \
      -DCMAKE_INSTALL_PREFIX=$HOME/tfhe-install
make -j4 && make install
```

and then:

```sh
make compare TFHE_PREFIX=$HOME/tfhe-install                      # spqlios-fma FFT
./build/bench_original 500
rm build/bench_original
make compare TFHE_PREFIX=$HOME/tfhe-install TFHE_FFT=nayuki-portable
./build/bench_original 500
```

### Report and slides

Every number in the report and in the slides comes from `report/data/`. To
measure everything again on another machine (about 15 minutes, single thread
pinned to one core) and rebuild both documents:

```sh
TFHE_PREFIX=$HOME/tfhe-install report/experiments/run_all.sh [cpu]
make -C report && make -C slides
```

## Parameters

| Parameter | Value | Role |
|-----------|-------|------|
| `n` | 630 | LWE dimension; number of external products per bootstrapping |
| `N` | 1024 | ring degree; the phase is rounded to 2N positions |
| `ell`, `Bgbit` | 3, 7 | gadget decomposition of the external product |
| `ksT`, `ksBasebit` | 8, 2 | decomposition used by key switching |
| `lweStdev` | 2^-15 | noise of LWE samples and of the key-switching key |
| `bkStdev` | 2^-25 | noise of the bootstrapping key |

## Limitations

- The random number generator (`std::mt19937_64`) is **not** cryptographically
  secure; it is adequate for studying and measuring the scheme only.
- The code makes no attempt at constant-time execution or side-channel
  resistance.
- Keys and ciphertexts are not serialized, so client and server run in the
  same process.
