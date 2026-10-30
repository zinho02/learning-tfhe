// adder.cpp -- Homomorphic n-bit ripple-carry adder.
//
// Usage: ./build/adder [nbits] [a] [b]
//
// main() is split into the three moments of the protocol: the client generates
// the keys and encrypts the inputs, the server adds them using only the
// evaluation key and the ciphertexts, and the client decrypts the result.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "../src/gates.hpp"

using namespace tfhe;

// An encrypted integer: one encrypted bit per position, least significant first.
using Cipher = std::vector<LweSample>;

Cipher encryptInt(uint32_t x, int nbits, const SecretKey& sk) {
    Cipher c;
    for (int i = 0; i < nbits; ++i) c.push_back(encryptBit((x >> i) & 1, sk));
    return c;
}

uint32_t decryptInt(const Cipher& c, const SecretKey& sk) {
    uint32_t x = 0;
    for (size_t i = 0; i < c.size(); ++i) x |= uint32_t{decryptBit(c[i], sk)} << i;
    return x;
}

// For each position i:
//     sum_i = x_i XOR y_i XOR carry
//     carry = (x_i XOR y_i) ? carry : x_i
// If x_i = y_i the new carry equals x_i whatever the old carry; otherwise the
// old carry propagates. Two XORs and one MUX per bit (four blind rotations).
// The initial carry is a noiseless trivial encryption of false, which the
// server can build without any key. The result has nbits + 1 bits.
Cipher add(const Cipher& x, const Cipher& y, const EvalKey& ek) {
    Cipher sum;
    LweSample carry = lweTrivial(-MU, ek.p.n);
    for (size_t i = 0; i < x.size(); ++i) {
        LweSample xy = gateXOR(x[i], y[i], ek);
        sum.push_back(gateXOR(xy, carry, ek));
        carry = gateMUX(xy, carry, x[i], ek);
    }
    sum.push_back(carry);
    return sum;
}

int main(int argc, char** argv) {
    int nbits = argc > 1 ? std::atoi(argv[1]) : 8;
    uint32_t a = argc > 2 ? std::atoi(argv[2]) : 173;
    uint32_t b = argc > 3 ? std::atoi(argv[3]) : 94;

    // client: key generation and encryption
    Params p;
    SecretKey sk = secretKeyGen(p);
    EvalKey ek = evalKeyGen(sk);
    Cipher ca = encryptInt(a, nbits, sk), cb = encryptInt(b, nbits, sk);

    // server: homomorphic addition, without access to sk
    auto t0 = std::chrono::steady_clock::now();
    Cipher cs = add(ca, cb, ek);
    auto t1 = std::chrono::steady_clock::now();

    // client: decryption
    uint32_t s = decryptInt(cs, sk);
    double sec = std::chrono::duration<double>(t1 - t0).count();
    std::printf("%u + %u = %u (%d bits, %.2f s)\n", a, b, s, nbits, sec);
    return s == a + b ? 0 : 1;
}
