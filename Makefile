# Build the unit tests (tests/), the applications (apps/) and, optionally, the
# benchmark linked against the original TFHE library (compare/).
#
#   make                 build every test and application into build/
#   make test            build and run the tests in order, stopping at the first failure
#   make compare         build build/bench_original (see TFHE_PREFIX and TFHE_FFT below)
#   make clean           remove build/

CXX ?= g++
CXXFLAGS ?= -std=c++17 -O3 -march=native -Wall -Wextra

HEADERS := $(wildcard src/*.hpp) tests/check.hpp
TESTS := test01_torus test02_lwe test03_poly test04_fft test05_trlwe \
			test06_trgsw test07_keys test08_blindrotate test09_bootstrap \
			test10_gates
APPS := adder bench

all: $(addprefix build/,$(TESTS) $(APPS))

# tests/NAME.cpp and apps/NAME.cpp become the executable build/NAME; every
# header is a dependency, since the whole library is header-only.
build/%: tests/%.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $< -o $@

build/%: apps/%.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $< -o $@

test: $(addprefix build/,$(TESTS))
	@for t in $(TESTS); do ./build/$$t || exit 1; done

clean:
	rm -rf build

# Original TFHE library, installed under TFHE_PREFIX; TFHE_FFT selects which of
# its FFT back ends to link (spqlios-fma or nayuki-portable).
TFHE_PREFIX ?= /usr/local
TFHE_FFT ?= spqlios-fma

compare: build/bench_original

build/bench_original: compare/bench_original.cpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $< -I$(TFHE_PREFIX)/include -L$(TFHE_PREFIX)/lib \
		-ltfhe-$(TFHE_FFT) -Wl,-rpath,$(TFHE_PREFIX)/lib -o $@

.PHONY: all test clean compare
