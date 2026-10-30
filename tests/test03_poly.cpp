// test03_poly.cpp -- Schoolbook product and negacyclic rotation modulo X^N + 1.
#include <cstdio>
#include "../src/poly.hpp"
#include "check.hpp"
using namespace tfhe;

int main() {
    setSeed(3);

    // 1) hand-checkable example with N = 4:
    //    (1 + 2X + 3X^3) * X^2 = X^2 + 2X^3 + 3X^5 = -3X + X^2 + 2X^3
    TorusPoly b = {1, 2, 0, 3};
    IntPoly x2 = {0, 0, 1, 0};
    TorusPoly expected = {0, Torus(-3), 1, 2};
    CHECK(mulNaive(x2, b) == expected);
    CHECK(mulByXe(b, 2) == expected);

    // 2) powers of X with N = 4: X^0 = 1, X^1, X^4 = -1, X^5 = -X
    CHECK(mulByXe(b, 0) == b);
    CHECK(mulByXe(b, 1) == (TorusPoly{Torus(-3), 1, 2, 0}));
    CHECK(mulByXe(b, 4) == (TorusPoly{Torus(-1), Torus(-2), 0, Torus(-3)}));
    CHECK(mulByXe(b, 5) == (TorusPoly{3, Torus(-1), Torus(-2), 0}));

    // 3) the same identities at full size (N = 1024):
    //    X^N = -1, X^(2N) = 1 and X^a * X^b = X^(a+b mod 2N)
    const int N = 1024;
    TorusPoly p(N);
    for (auto& c : p) c = uniformTorus();

    TorusPoly q = mulByXe(p, N);
    for (int j = 0; j < N; ++j) CHECK(q[j] == -p[j]);

    CHECK(mulByXe(mulByXe(p, 700), 2 * N - 700) == p);
    CHECK(mulByXe(mulByXe(p, 300), 500) == mulByXe(p, 800));
    CHECK(mulByXe(mulByXe(p, 900), 1500) == mulByXe(p, 2400 % (2 * N)));

    // 4) rotating by X^5 equals multiplying by the monomial X^5
    IntPoly mono(N, 0);
    mono[5] = 1;
    CHECK(mulNaive(mono, p) == mulByXe(p, 5));

    // 5) distributivity: a * (p + q) = a * p + a * q
    IntPoly a(N);
    for (auto& c : a) c = static_cast<int32_t>(rng()() % 128) - 64;
    TorusPoly pq = p;
    polyAddTo(pq, q);
    TorusPoly lhs = mulNaive(a, pq);
    TorusPoly rhs = mulNaive(a, p);
    polyAddTo(rhs, mulNaive(a, q));
    CHECK(lhs == rhs);

    std::printf("test03_poly: OK\n");
    return 0;
}
