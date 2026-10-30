// test01_torus.cpp -- Torus conversions, arithmetic, encoding and noise sampling.
#include <cmath>
#include <cstdio>
#include "../src/torus.hpp"
#include "check.hpp"
using namespace tfhe;

int main() {
    setSeed(1);

    // 1) real -> torus, including reduction mod 1
    CHECK(dtot(0.25) == 0x40000000u);
    CHECK(dtot(-0.25) == 0xC0000000u);
    CHECK(dtot(1.25) == 0x40000000u);

    // 2) torus arithmetic: addition wraps around; integer times torus
    CHECK(dtot(0.75) + dtot(0.5) == dtot(0.25));
    CHECK(3 * dtot(0.375) == dtot(0.125));

    // 3) torus -> real in [-1/2, 1/2)
    CHECK(ttod(dtot(0.75)) == -0.25);
    CHECK(std::fabs(ttod(dtot(0.123)) - 0.123) < 1e-9);

    // 4) encoding in Z_8: decoding tolerates errors below 1/16
    for (int mu = 0; mu < 8; ++mu) {
        Torus t = encode(mu, 8);
        CHECK(decode(t, 8) == mu);
        CHECK(decode(t + dtot(0.06), 8) == mu);
        CHECK(decode(t - dtot(0.06), 8) == mu);
        CHECK(decode(t + dtot(0.07), 8) == (mu + 1) % 8);
    }

    // 5) Gaussian noise: mean close to 0 and the requested standard deviation
    const double stdev = std::ldexp(1.0, -15);
    const int reps = 100000;
    double sum = 0, sum2 = 0;
    for (int i = 0; i < reps; ++i) {
        double e = ttod(gaussianTorus(stdev));
        sum += e;
        sum2 += e * e;
    }
    double mean = sum / reps, sd = std::sqrt(sum2 / reps - mean * mean);
    std::printf("noise: requested standard deviation %.3e, measured %.3e\n", stdev, sd);
    CHECK(std::fabs(mean) < 0.02 * stdev);
    CHECK(std::fabs(sd / stdev - 1.0) < 0.02);

    // 6) uniform bits: about half of them are 1
    int ones = 0;
    for (int i = 0; i < reps; ++i) ones += uniformBit();
    CHECK(ones > 0.49 * reps && ones < 0.51 * reps);

    std::printf("test01_torus: OK\n");
    return 0;
}
