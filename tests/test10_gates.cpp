// test10_gates.cpp -- Truth tables, deep composition and a small circuit.
#include <cstdio>
#include "../src/gates.hpp"
#include "check.hpp"
using namespace tfhe;

int main() {
    setSeed(10);
    Params p;
    SecretKey sk = secretKeyGen(p);
    EvalKey ek = evalKeyGen(sk);

    // 1) every gate on every input, 10 times with fresh ciphertexts
    for (int rep = 0; rep < 10; ++rep) {
        for (int x = 0; x <= 1; ++x) {
            for (int y = 0; y <= 1; ++y) {
                LweSample cx = encryptBit(x, sk), cy = encryptBit(y, sk);
                CHECK(decryptBit(gateNAND(cx, cy, ek), sk) == !(x && y));
                CHECK(decryptBit(gateAND(cx, cy, ek), sk) == (x && y));
                CHECK(decryptBit(gateOR(cx, cy, ek), sk) == (x || y));
                CHECK(decryptBit(gateNOR(cx, cy, ek), sk) == !(x || y));
                CHECK(decryptBit(gateXOR(cx, cy, ek), sk) == (x != y));
                CHECK(decryptBit(gateXNOR(cx, cy, ek), sk) == (x == y));
                CHECK(decryptBit(gateNOT(cx), sk) == !x);
                for (int s = 0; s <= 1; ++s) {
                    LweSample cs = encryptBit(s, sk);
                    CHECK(decryptBit(gateMUX(cs, cx, cy, ek), sk) == (s ? x : y));
                }
            }
        }
    }

    // 2) unbounded depth: NAND(x, 1) = NOT x applied 50 times in sequence
    bool m = true;
    LweSample c = encryptBit(m, sk), one = encryptBit(true, sk);
    for (int k = 0; k < 50; ++k) {
        c = gateNAND(c, one, ek);
        m = !m;
        CHECK(decryptBit(c, sk) == m);
    }

    // 3) majority of three bits: maj(x, y, z) = (x AND y) OR (z AND (x XOR y))
    for (int bits = 0; bits < 8; ++bits) {
        bool x = bits & 1, y = bits & 2, z = bits & 4;
        LweSample cx = encryptBit(x, sk), cy = encryptBit(y, sk), cz = encryptBit(z, sk);
        LweSample maj = gateOR(gateAND(cx, cy, ek),
                                gateAND(cz, gateXOR(cx, cy, ek), ek), ek);
        CHECK(decryptBit(maj, sk) == (x + y + z >= 2));
    }

    std::printf("test10_gates: OK\n");
    return 0;
}
