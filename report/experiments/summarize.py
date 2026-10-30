#!/usr/bin/env python3
"""summarize.py -- Turns the raw measurements into the data used by the report.

Reads report/data/raw/ (written by run_all.sh) and writes
  report/data/results.tex    LaTeX macros with every number quoted in the text
  report/data/time.dat       time per gate and key generation (Figure 'time')
  report/data/profile.dat    cost breakdown of one gate (Figure 'profile')
  report/data/noise_hist.dat histogram of the output error (Figure 'noise')
  report/data/noise_keys.dat per-key mean and variance (Figure 'noise')
"""
import glob
import math
import os
import re
import shutil
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
RAW = os.path.join(DATA, "raw")
SRC = os.path.join(HERE, "..", "..", "src")

# Parameters (src/params.hpp).
n, N, ell, Bgbit, t, basebit = 630, 1024, 3, 7, 8, 2
var_lwe, var_bk = 2.0 ** -30, 2.0 ** -50
Bks = 2 ** basebit


def num(pattern, text, k=1):
    return float(re.search(pattern, text).group(k))


def read(name):
    with open(os.path.join(RAW, name)) as f:
        return f.read()


def rss_mib(name):
    return num(r"Maximum resident set size \(kbytes\): (\d+)", read(name)) / 1024


def runs(prefix):
    return sorted(glob.glob(os.path.join(RAW, prefix + "_*.txt")))


# ---------------------------------------------------------------- this work
ours = []
for path in runs("bench"):
    txt = open(path).read()
    ours.append(dict(
        keygen=num(r"evaluation key: ([\d.]+) ms", txt) / 1000,
        nand=num(r"NAND: ([\d.]+) ms", txt),
        br=num(r"extraction: ([\d.]+) ms", txt),
        ks=num(r"key switching: ([\d.]+) ms", txt),
        meanKS=num(r"after key switching: mean ([-+\d.e]+)", txt),
        varKS=num(r"after key switching: mean [-+\d.e]+ variance ([\d.e-]+)", txt),
        errors=int(num(r"decryption errors: (\d+) out of (\d+)", txt)),
        gates=int(num(r"decryption errors: (\d+) out of (\d+)", txt, 2)),
        rss=rss_mib(os.path.basename(path).replace(".txt", ".time")),
    ))

# ---------------------------------------------------------- reference library
def orig(prefix):
    res = []
    for path in runs(prefix):
        txt = open(path).read()
        res.append(dict(
            keygen=num(r"key generation: ([\d.]+) ms", txt) / 1000,
            nand=num(r"NAND: ([\d.]+) ms", txt),
            mean=num(r"output error: mean ([-+\d.e]+)", txt),
            var=num(r"output error: mean [-+\d.e]+ variance ([\d.e-]+)", txt),
            errors=int(num(r"decryption errors: (\d+) out of (\d+)", txt)),
            gates=int(num(r"decryption errors: (\d+) out of (\d+)", txt, 2)),
            rss=rss_mib(os.path.basename(path).replace(".txt", ".time")),
        ))
    return res


spq, nay = orig("orig_spq"), orig("orig_nay")

# ------------------------------------------------------------ other experiments
adder = re.findall(r"(\d+) \+ (\d+) = (\d+) \((\d+) bits, ([\d.]+) s\)", read("adder.txt"))
fft = read("test04_fft.txt")
fft_naive = num(r"schoolbook product: ([\d.]+) us", fft)
fft_fast = num(r"FFT product: ([\d.]+) us", fft)
t09 = read("test09_bootstrap.txt")
boot_in = num(r"noise before bootstrapping: ([\d.e+-]+)", t09)
boot_out = num(r"after: ([\d.e+-]+)", t09)
profile = dict(line.split() for line in read("profile.txt").splitlines()[1:])
profile = {k: float(v) for k, v in profile.items()}
plog = read("profile.log")
mac_hot = num(r"mulAcc ([\d.]+)", plog)
mac_key = num(r"key from memory ([\d.]+)", plog)
ext_hot = num(r"same key sample ([\d.]+)", plog)
ext_key = num(r"cycling over the key ([\d.]+)", plog)
prim = re.search(r"forward ([\d.]+) backward ([\d.]+) mulAcc [\d.]+ \(key from memory [\d.]+\) "
                 r"decompose ([\d.]+)", plog)
fwd_us, bwd_us, dec_us = (float(prim.group(i)) for i in (1, 2, 3))
ext_parts = 2 * ell * fwd_us + 2 * bwd_us + 4 * ell * mac_hot + 2 * dec_us

keys = [line.split() for line in read("noise_keys.dat").splitlines()[1:]]
noise_var_br = st.mean(float(r[2]) for r in keys)
noise_var_ks = st.mean(float(r[4]) for r in keys)
noise_means = [float(r[3]) for r in keys]
noise_bias_std = st.stdev(noise_means)
noise_errors = sum(int(r[5]) for r in keys)
noise_gates = 1000

cpu = read("cpu.txt").strip().replace("(R)", "").replace("(TM)", "")
compiler = read("compiler.txt").strip()
m = re.match(r"g\+\+ \(GCC\) ([\d.]+)", compiler)
compiler = "GCC " + m.group(1) if m else compiler

# ------------------------------------------------------------------- theory
beta = 2 ** (Bgbit - 1)
eps = 2.0 ** -(ell * Bgbit + 1)
bound_br = n * 2 * ell * N * beta ** 2 * var_bk + n * (1 + N) * eps ** 2
bound_ks = t * N * var_lwe + N * 2.0 ** (-2 * (t * basebit + 1))
Ed2 = sum(d * d for d in range(-beta, beta)) / (2 * beta)
avg_br = n * 2 * ell * N * Ed2 * var_bk
avg_ks = t * N * ((Bks - 1) / Bks) ** 2 * var_lwe
bias_std = math.sqrt(t * N * (Bks - 1) / Bks ** 2 * var_lwe)

# Lines of library code, without comments and blank lines.
lib_lines = 0
for path in glob.glob(os.path.join(SRC, "*.hpp")):
    for line in open(path):
        s = line.strip()
        if s and not s.startswith("//"):
            lib_lines += 1


# ----------------------------------------------------------------- formatting
def sci(x, digits=2):
    """x as a LaTeX number in scientific notation."""
    e = int(math.floor(math.log10(abs(x))))
    return r"\ensuremath{%.*f\cdot10^{%d}}" % (digits - 1, x / 10 ** e, e)


def rng(vals, fmt):
    lo, hi = min(vals), max(vals)
    a, b = fmt % lo, fmt % hi
    return a if a == b else a + "--" + b


def thousands(k):
    return "{:,}".format(k).replace(",", "{,}")


mean = lambda rs, key: st.mean(r[key] for r in rs)
our_nand, spq_nand, nay_nand = mean(ours, "nand"), mean(spq, "nand"), mean(nay, "nand")
total_gates = sum(r["gates"] for r in ours) + noise_gates * len(keys)
assert sum(r["errors"] for r in ours) == 0 and noise_errors == 0

total = profile["total"]
profile_steps = num(r"steps (\d+)", plog)
share = lambda *parts: 100 * sum(profile[p] for p in parts) / total

macros = {
    "OurNandMean": "%.1f" % our_nand,
    "SlowdownSpq": "%.1f" % (our_nand / spq_nand),
    "SpeedupNay": "%.1f" % (nay_nand / our_nand),
    "TotalGates": thousands(total_gates),
    "OurBenchGates": thousands(sum(r["gates"] for r in ours)),
    "OrigGates": thousands(sum(r["gates"] for r in spq + nay)),
    "OrigGatesEach": thousands(sum(r["gates"] for r in spq)),
    "OurRuns": str(len(ours)),
    "OrigRuns": str(len(spq)),
    "OurKS": "%.1f" % max(r["ks"] for r in ours),
    "OurNandRange": rng([r["nand"] for r in ours], "%.1f"),
    "SpqNandRange": rng([r["nand"] for r in spq], "%.1f"),
    "NayNandRange": rng([r["nand"] for r in nay], "%.1f"),
    "OurKeygenRange": rng([r["keygen"] for r in ours], "%.2f"),
    "SpqKeygenRange": rng([r["keygen"] for r in spq], "%.2f"),
    "NayKeygenRange": rng([r["keygen"] for r in nay], "%.2f"),
    "OurVarRange": rng([r["varKS"] * 1e5 for r in ours], "%.2f"),
    "SpqVarRange": rng([r["var"] * 1e5 for r in spq], "%.2f"),
    "NayVarRange": rng([r["var"] * 1e5 for r in nay], "%.2f"),
    "OurMeanMax": sci(max(abs(r["meanKS"]) for r in ours), 2),
    "SpqMeanMax": sci(max(abs(r["mean"]) for r in spq), 2),
    "NayMeanMax": sci(max(abs(r["mean"]) for r in nay), 2),
    "OurRSS": "%.0f" % mean(ours, "rss"),
    "SpqRSS": "%.0f" % mean(spq, "rss"),
    "NayRSS": "%.0f" % mean(nay, "rss"),
    "BkMB": "%.1f" % (n * 2 * ell * 2 * (N // 2) * 16 / 2 ** 20),
    "KskMB": "%.1f" % (N * t * Bks * (n + 1) * 4 / 2 ** 20),
    "AdderEight": adder[0][4], "AdderSixteen": adder[1][4], "AdderTwentyFour": adder[2][4],
    "AdderEightText": r"$%s+%s=%s$ (%s bits)" % adder[0][:4],
    "AdderSixteenText": r"$%s+%s=%s$ (%s bits)" % adder[1][:4],
    "AdderTwentyFourText": r"$%s+%s=%s$ (%s bits)" % adder[2][:4],
    "FftNaive": "%.0f" % fft_naive,
    "FftFast": "%.1f" % fft_fast,
    "FftSpeedup": "%.0f" % (fft_naive / fft_fast),
    "BootSdIn": sci(boot_in), "BootSdOut": sci(boot_out),
    "BootSdInRaw": "%.6e" % boot_in, "BootSdOutRaw": "%.6e" % boot_out,
    "BRShare": "%.1f" % (100 - share("keyswitch")),
    "FftShare": "%.0f" % share("fft-forward", "fft-backward"),
    "MacShare": "%.0f" % share("pointwise"),
    "FreqShare": "%.0f" % share("fft-forward", "fft-backward", "pointwise"),
    "DecShare": "%.0f" % share("decompose"),
    "OtherShare": "%.0f" % share("other"),
    "ExtHot": "%.1f" % ext_hot, "ExtKey": "%.1f" % ext_key,
    "ExtParts": "%.1f" % ext_parts, "ExtAlloc": "%.0f" % (ext_hot - ext_parts),
    "MacHot": "%.2f" % mac_hot, "MacKey": "%.2f" % mac_key,
    "MacRatio": "%.1f" % (mac_key / mac_hot),
    "MemShare": "%.0f" % (100 * (ext_key - ext_hot) * profile_steps / total),
    "BoundBR": sci(bound_br), "BoundTotal": sci(bound_br + bound_ks),
    "AvgBR": sci(avg_br), "AvgKS": sci(avg_ks), "AvgTotal": sci(avg_br + avg_ks),
    "BiasStd": sci(bias_std),
    "NoiseKeys": str(len(keys)), "NoiseGates": thousands(noise_gates),
    "NoiseVarBR": sci(noise_var_br), "NoiseVarKS": sci(noise_var_ks),
    "NoiseBiasStd": sci(noise_bias_std),
    "CPU": cpu, "Compiler": compiler,
    "LibLines": str(lib_lines),
    # raw values used by the plots
    "AvgTotalRaw": "%.6e" % (avg_br + avg_ks), "NoiseVarKSRaw": "%.6e" % noise_var_ks,
    "BiasStdRaw": "%.6e" % bias_std,
}

with open(os.path.join(DATA, "results.tex"), "w") as f:
    f.write("% Generated by experiments/summarize.py -- do not edit.\n")
    for k, v in macros.items():
        f.write("\\newcommand{\\%s}{%s}\n" % (k, v))

with open(os.path.join(DATA, "time.dat"), "w") as f:
    f.write("impl x nand nandlo nandhi keygen keygenlo keygenhi\n")
    for i, (name, rs, kg) in enumerate([("ours", ours, "keygen"), ("spq", spq, "keygen"),
                                         ("nay", nay, "keygen")]):
        v = [r["nand"] for r in rs]
        g = [r[kg] for r in rs]
        f.write("%s %d %.3f %.3f %.3f %.4f %.4f %.4f\n" % (
            name, i, st.mean(v), st.mean(v) - min(v), max(v) - st.mean(v),
            st.mean(g), st.mean(g) - min(g), max(g) - st.mean(g)))

order = ["fft-forward", "pointwise", "fft-backward", "other", "keyswitch", "rotate", "decompose"]
with open(os.path.join(DATA, "profile.dat"), "w") as f:
    f.write("part y ms pct\n")
    for i, p in enumerate(reversed(order)):
        f.write("%s %d %.3f %.1f\n" % (p, i, profile[p] / 1000, 100 * profile[p] / total))

shutil.copy(os.path.join(RAW, "noise_hist.dat"), os.path.join(DATA, "noise_hist.dat"))
with open(os.path.join(DATA, "noise_keys.dat"), "w") as f:
    f.write("key meanKS\n")
    for r in keys:
        f.write("%d %.6e\n" % (int(r[0]) + 1, float(r[3])))

for k in ["OurNandMean", "SlowdownSpq", "SpeedupNay", "TotalGates", "NoiseVarBR",
          "NoiseVarKS", "NoiseBiasStd", "AvgBR", "AvgTotal", "BiasStd", "FreqShare",
          "OtherShare", "LibLines"]:
    print("%-14s %s" % (k, macros[k]))
