// fft.hpp -- Fast multiplication modulo X^N + 1 through a negacyclic FFT.
//
// A product modulo X^N + 1 is computed by evaluating both factors at the
// roots of X^N + 1, multiplying pointwise and interpolating back. Since the
// coefficients are real, only N/2 roots (one per conjugate pair) are needed:
// omega_k = exp(i*pi*(4k+1)/N), k = 0..N/2-1. Using omega_k^(N/2) = i, the
// evaluation reduces to a standard complex FFT of size M = N/2:
//
//     forward:  fold z_j = p_j + i*p_(j+N/2),  twist z_j *= exp(i*pi*j/N),  FFT
//     backward: inverse FFT,  untwist,  unfold (real -> p_j, imag -> p_(j+N/2))
//
// Precision: one factor always has small coefficients (|digit| <= 64 or key
// bits) and the other has torus coefficients read as signed 32-bit integers,
// so every output coefficient is below 2^6 * 2^31 * 2^10 = 2^47 in absolute
// value, which fits in the 53-bit mantissa of a double.
#pragma once
#include <cmath>
#include <complex>
#include <utility>
#include <vector>
#include "poly.hpp"

namespace tfhe {

using cplx = std::complex<double>;

// A polynomial in the frequency domain: its evaluations at N/2 roots of
// X^N + 1. Addition and multiplication are pointwise in this representation.
using FreqPoly = std::vector<cplx>;

// Complex product written by hand: operator* of std::complex performs
// NaN/infinity checks that make the FFT loops much slower.
inline cplx cmul(cplx x, cplx y) {
    return {x.real() * y.real() - x.imag() * y.imag(),
            x.real() * y.imag() + x.imag() * y.real()};
}

// Negacyclic transform for a fixed power-of-two N. The constructor
// precomputes every table that depends only on N.
class NegFFT {
    public:
        explicit NegFFT(int N) : N_(N), M_(N/2), twist_(M_), root_(M_/2), rootInv_(M_/2), rev_(M_) {
            const double pi = std::acos(-1.0);
            for (int j = 0; j < M_; ++j) twist_[j] = std::polar(1.0, pi * j / N_);
            for (int k = 0; k < M_/2; ++k) {
                root_[k] = std::polar(1.0, 2 * pi * k / M_);
                rootInv_[k] = std::conj(root_[k]);
            }

            int logM = 0;
            while ((1 << logM) < M_) ++logM;
            for (int i = 0; i < M_; ++i) {
                int r = 0;
                for (int b = 0; b < logM; ++b) if (i & (1 << b)) r |= 1 << (logM - 1 - b);
                rev_[i] = r;
            }
        }
        int N() const { return N_; }

        // Coefficients -> evaluations. Works for IntPoly and TorusPoly; torus
        // coefficients are read as signed integers.
        template <class T>
        void forward(FreqPoly& out, const std::vector<T>& p) const {
            out.resize(M_);
            for (int j = 0; j < M_; ++j) {
                cplx z(static_cast<int32_t>(p[j]), static_cast<int32_t>(p[j + M_]));
                out[j] = cmul(z, twist_[j]);
            }
            transform(out, root_);
        }

        // Evaluations -> torus coefficients. The input is used as scratch
        // space and is destroyed, which saves a copy on the hot path.
        void backward(TorusPoly& p, FreqPoly& in) const {
            transform(in, rootInv_);
            p.resize(N_);
            const double scale = 1.0 / M_;
            for (int j = 0; j < M_; ++j) {
                cplx z = cmul(in[j], std::conj(twist_[j]));
                p[j] = toTorus(z.real() * scale);
                p[j + M_] = toTorus(z.imag() * scale);
            }
        }
    private:
        // Rounds to the nearest integer and reduces modulo 2^32. The value goes
        // through int64_t because it may be as large as ~2^47.
        static Torus toTorus(double x) {
            return static_cast<Torus>(static_cast<int64_t>(x < 0 ? x - 0.5 : x + 0.5));
        }

        // In-place iterative radix-2 Cooley-Tukey FFT of size M: bit-reversal
        // permutation followed by log2(M) rounds of butterflies
        // (u, v) -> (u + w*v, u - w*v). The table w is root_ for the forward
        // transform and rootInv_ for the inverse one.
        void transform(FreqPoly& a, const std::vector<cplx>& w) const {
            for (int i = 0; i < M_; ++i)
                if (i < rev_[i]) std::swap(a[i], a[rev_[i]]);
            for (int len = 2; len <= M_; len <<= 1) {
                int half = len / 2, step = M_ / len;
                for (int i = 0; i < M_; i += len) {
                    for (int k = 0; k < half; ++k) {
                        cplx u = a[i + k];
                        cplx v = cmul(a[i + k + half], w[k * step]);
                        a[i + k] = u + v;
                        a[i + k + half] = u - v;
                    }
                }
            }
        }

        int N_, M_;                     // ring degree and FFT size (N/2)
        std::vector<cplx> twist_;       // exp(i*pi*j/N),    j = 0..M-1
        std::vector<cplx> root_;        // exp(+2*pi*i*k/M), k = 0..M/2-1
        std::vector<cplx> rootInv_;     // exp(-2*pi*i*k/M)
        std::vector<int> rev_;          // bit-reversal permutation
};

// acc <- acc + x * y, pointwise in the frequency domain. Used by the external
// product to accumulate several products before a single inverse transform.
inline void mulAcc(FreqPoly& acc, const FreqPoly& x, const FreqPoly& y) {
    for (size_t k = 0; k < acc.size(); ++k) acc[k] += cmul(x[k], y[k]);
}

// r = a * b mod X^N + 1 in O(N log N): two forward transforms, a pointwise
// product and one backward transform.
inline TorusPoly mulFFT(const NegFFT& fft, const IntPoly& a, const TorusPoly& b) {
    FreqPoly fa, fb;
    fft.forward(fa, a);
    fft.forward(fb, b);
    for (size_t k = 0; k < fa.size(); ++k) fa[k] = cmul(fa[k], fb[k]);
    TorusPoly r;
    fft.backward(r, fa);
    return r;
}

}  // namespace tfhe
