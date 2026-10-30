// trgsw.hpp -- Gadget decomposition, TRGSW samples, external product and CMux.
//
// Multiplying a TRLWE sample by a large integer amplifies its noise, and two
// torus elements cannot be multiplied at all. The external product solves
// both problems: the TRLWE operand is decomposed into small signed digits
// (gadget decomposition) and only those digits multiply the rows of a TRGSW
// sample, which hold encryptions of m and -m*s' scaled by every gadget power.
#pragma once
#include <vector>
#include "trlwe.hpp"

namespace tfhe {

// Signed decomposition in base Bg = 2^Bgbit with ell levels:
//     t ~ d_1 / Bg + d_2 / Bg^2 + ... + d_ell / Bg^ell,   d_i in [-Bg/2, Bg/2),
// with error at most 1 / (2 Bg^ell). With Bg = 2^7 and ell = 3 it reads the 21
// most significant bits of t in three 7-bit blocks.
//
// Offset trick: balanced digits are obtained with shifts and masks only. Add
// to t the constant with Bg/2 in every digit position, extract the unsigned
// digits of the sum and subtract Bg/2 from each of them. Adding and then
// removing Bg/2 per position leaves the represented value unchanged, but the
// carries of the addition propagate between positions and perform exactly the
// "borrow one from the next digit" adjustment of a balanced representation.
// Half of the least significant digit is also added, so that the discarded
// low-order bits are rounded instead of truncated.
struct Gadget {
    int ell, Bgbit;
    uint32_t halfBg;   // Bg / 2
    uint32_t mask;     // Bg - 1: selects the Bgbit bits of one digit
    uint32_t offset;   // constant added before extracting the digits

    Gadget(int ell_, int Bgbit_) : ell(ell_), Bgbit(Bgbit_) {
        halfBg = 1u << (Bgbit - 1);
        mask = (1u << Bgbit) - 1;
        offset = 0;

        for (int i = 0; i < ell; ++i) offset += halfBg << (32 - (i + 1) * Bgbit);

        offset += 1u << (32 - ell * Bgbit - 1);
    }

    // The torus element 1 / Bg^(i+1) (levels are counted from 0 in the code).
    Torus h(int i) const { return 1u << (32 - (i + 1) * Bgbit); }
};

// Decomposes every coefficient of p: out[i][j] is the level-i digit of
// coefficient j. The result is ell polynomials with small coefficients.
inline void decompose(std::vector<IntPoly>& out, const TorusPoly& p, const Gadget& g) {
    int N = static_cast<int>(p.size());
    out.assign(g.ell, IntPoly(N));
    for (int j = 0; j < N; ++j) {
        uint32_t t = p[j] + g.offset;
        for (int i = 0; i < g.ell; ++i) {
            uint32_t digit = (t >> (32 - (i + 1) * g.Bgbit)) & g.mask;
            out[i][j] = static_cast<int32_t>(digit) - static_cast<int32_t>(g.halfBg);
        }
    }
}

// A TRGSW encryption of an integer m is 2*ell TRLWE encryptions of zero plus
// m times the gadget matrix:
//     row i         (0 <= i < ell): (a_i + m / Bg^(i+1), b_i)
//     row ell + i   (0 <= i < ell): (a_i',               b_i' + m / Bg^(i+1))
// The first ell rows encrypt -m*s' / Bg^(i+1) and the last ell rows encrypt
// m / Bg^(i+1).
struct TrgswSample {
    std::vector<TrlweSample> rows;
};

// The same ciphertext with every polynomial already in the frequency domain.
// The bootstrapping key is stored in this form, so its forward transforms are
// paid once at key generation and never during gate evaluation.
struct TrgswFFT {
    std::vector<FreqPoly> a, b;
};

// Encrypts the integer m (in this project, always a bit of the LWE key).
// m / Bg^(i+1) is a constant polynomial, so only coefficient 0 changes.
inline TrgswSample trgswEncrypt(int32_t m, double stdev, const TrlweKey& key,
                                const Gadget& g, const NegFFT& fft) {
    TorusPoly zero(fft.N(), 0);
    TrgswSample C;
    for (int r = 0; r < 2 * g.ell; ++r)
        C.rows.push_back(trlweEncrypt(zero, stdev, key, fft));
    for (int i = 0; i < g.ell; ++i) {
        C.rows[i].a[0] += m * g.h(i);
        C.rows[g.ell + i].b[0] += m * g.h(i);
    }
    return C;
}

// Transforms every row of a TRGSW sample to the frequency domain.
inline TrgswFFT toFFT(const TrgswSample& C, const NegFFT& fft) {
    TrgswFFT F;
    F.a.resize(C.rows.size());
    F.b.resize(C.rows.size());
    for (size_t r = 0; r < C.rows.size(); ++r) {
        fft.forward(F.a[r], C.rows[r].a);
        fft.forward(F.b[r], C.rows[r].b);
    }
    return F;
}

// External product C [.] d: if C encrypts the integer m and d = (a, b) encrypts
// the polynomial mu, the result is a TRLWE encryption of m * mu, computed as
//     sum_i dec_i(a) * row_i  +  sum_i dec_i(b) * row_(ell+i).
// Products are accumulated in the frequency domain, so the cost is 2*ell
// forward transforms plus only 2 backward transforms.
inline TrlweSample externalProduct(const TrgswFFT& C, const TrlweSample& d,
                                    const Gadget& g, const NegFFT& fft) {
    int N = fft.N();
    FreqPoly accA(N / 2, 0), accB(N / 2, 0);
    FreqPoly tmp;
    std::vector<IntPoly> dec;

    decompose(dec, d.a, g);
    for (int i = 0; i < g.ell; ++i) {
        fft.forward(tmp, dec[i]);
        mulAcc(accA, tmp, C.a[i]);
        mulAcc(accB, tmp, C.b[i]);
    }

    decompose(dec, d.b, g);
    for (int i = 0; i < g.ell; ++i) {
        fft.forward(tmp, dec[i]);
        mulAcc(accA, tmp, C.a[g.ell + i]);
        mulAcc(accB, tmp, C.b[g.ell + i]);
    }

    TrlweSample r(N);
    fft.backward(r.a, accA);
    fft.backward(r.b, accB);
    return r;
}

// Encrypted selector: returns an encryption of d1 if C encrypts 1 and of d0 if
// C encrypts 0, computed as d0 + C [.] (d1 - d0).
inline TrlweSample cmux(const TrgswFFT& C, const TrlweSample& d1,
                        const TrlweSample& d0, const Gadget& g,
                        const NegFFT& fft) {
    TrlweSample diff = d1;
    polySubTo(diff.a, d0.a);
    polySubTo(diff.b, d0.b);
    TrlweSample r = externalProduct(C, diff, g, fft);
    polyAddTo(r.a, d0.a);
    polyAddTo(r.b, d0.b);
    return r;
}

}  // namespace tfhe
