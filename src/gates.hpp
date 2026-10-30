// gates.hpp -- Homomorphic logic gates.
//
// Every binary gate is a linear combination of its two inputs, chosen so that
// the phase is positive exactly when the gate output is true, followed by a
// bootstrapping that maps the result back to +-1/8 with fresh noise:
//
//     NAND  (0,  1/8) - x - y          AND   (0, -1/8) + x + y
//     OR    (0,  1/8) + x + y          NOR   (0, -1/8) - x - y
//     XOR   (0,  1/4) + 2(x + y)       XNOR  (0, -1/4) - 2(x + y)
//
// For XOR, x + y is -1/4, 0 or +1/4 and the middle case would sit on the
// decision boundary; doubling maps both equal-input cases to 1/2 = -1/2,
// and adding 1/4 separates them from the unequal-input case.
#pragma once
#include "bootstrap.hpp"

namespace tfhe {

inline LweSample gateNAND(const LweSample& x, const LweSample& y, const EvalKey& ek) {
    LweSample t = lweTrivial(MU, ek.p.n);
    lweSubTo(t, x);
    lweSubTo(t, y);
    return bootstrap(MU, t, ek);
}

inline LweSample gateAND(const LweSample& x, const LweSample& y, const EvalKey& ek) {
    LweSample t = lweTrivial(-MU, ek.p.n);
    lweAddTo(t, x);
    lweAddTo(t, y);
    return bootstrap(MU, t, ek);
}

inline LweSample gateOR(const LweSample& x, const LweSample& y, const EvalKey& ek) {
    LweSample t = lweTrivial(MU, ek.p.n);
    lweAddTo(t, x);
    lweAddTo(t, y);
    return bootstrap(MU, t, ek);
}

inline LweSample gateNOR(const LweSample& x, const LweSample& y, const EvalKey& ek) {
    LweSample t = lweTrivial(-MU, ek.p.n);
    lweSubTo(t, x);
    lweSubTo(t, y);
    return bootstrap(MU, t, ek);
}

inline LweSample gateXOR(const LweSample& x, const LweSample& y, const EvalKey& ek) {
    LweSample t = lweTrivial(2 * MU, ek.p.n);
    lweAddMulTo(t, 2, x);
    lweAddMulTo(t, 2, y);
    return bootstrap(MU, t, ek);
}

inline LweSample gateXNOR(const LweSample& x, const LweSample& y, const EvalKey& ek) {
    LweSample t = lweTrivial(-2 * MU, ek.p.n);
    lweAddMulTo(t, -2, x);
    lweAddMulTo(t, -2, y);
    return bootstrap(MU, t, ek);
}

// NOT needs no bootstrapping: negating the sample swaps +1/8 and -1/8, and the
// noise only changes sign.
inline LweSample gateNOT(const LweSample& x) {
    LweSample r(static_cast<int>(x.a.size()));
    lweSubTo(r, x);
    return r;
}

// MUX(s, x, y) = s ? x : y = (s AND x) OR ((NOT s) AND y). The two conjunctions
// are never true together, so the OR becomes a plain sum: two bootstrappings
// without key switching (outputs under s'), sum (-1/4 or 0), plus 1/8, and a
// single key switching. Two blind rotations instead of three.
inline LweSample gateMUX(const LweSample& s, const LweSample& x,
                        const LweSample& y, const EvalKey& ek) {
    LweSample t1 = lweTrivial(-MU, ek.p.n);
    lweAddTo(t1, s);
    lweAddTo(t1, x);
    LweSample t2 = lweTrivial(-MU, ek.p.n);
    lweSubTo(t2, s);
    lweAddTo(t2, y);

    LweSample u = bootstrapNoKS(MU, t1, ek);
    lweAddTo(u, bootstrapNoKS(MU, t2, ek));
    u.b += MU;
    return keySwitch(u, ek.ksk);
}

}  // namespace tfhe
