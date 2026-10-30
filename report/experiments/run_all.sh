#!/bin/bash
# run_all.sh -- Reproduces every measurement used in the report.
#
# Usage (from the repository root):
#   TFHE_PREFIX=$HOME/tfhe-install report/experiments/run_all.sh [cpu]
#
# Runs single-threaded and pinned to one core (default: CPU 2; pick a
# performance core on hybrid processors). Raw outputs go to report/data/raw/;
# experiments/summarize.py then turns them into report/data/*.dat and
# report/data/results.tex. Takes about 15 minutes.
set -e
CPU=${1:-2}
TFHE_PREFIX=${TFHE_PREFIX:-$HOME/tfhe-install}
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
RAW=$ROOT/report/data/raw
BIN=$ROOT/build
cd "$ROOT"
mkdir -p "$RAW"
PIN="taskset -c $CPU"
CXXFLAGS="-std=c++17 -O3 -march=native -Wall -Wextra"

make -s all
g++ $CXXFLAGS compare/bench_original.cpp -I"$TFHE_PREFIX/include" -L"$TFHE_PREFIX/lib" \
    -ltfhe-spqlios-fma -Wl,-rpath,"$TFHE_PREFIX/lib" -o "$BIN/bench_original_spq"
g++ $CXXFLAGS compare/bench_original.cpp -I"$TFHE_PREFIX/include" -L"$TFHE_PREFIX/lib" \
    -ltfhe-nayuki-portable -Wl,-rpath,"$TFHE_PREFIX/lib" -o "$BIN/bench_original_nay"
g++ $CXXFLAGS report/experiments/profile.cpp -o "$BIN/profile"
g++ $CXXFLAGS report/experiments/noise.cpp -o "$BIN/noise"

for k in 1 2 3 4 5; do
    /usr/bin/time -v $PIN "$BIN/bench" 500 > "$RAW/bench_$k.txt" 2> "$RAW/bench_$k.time"
done
for k in 1 2 3; do
    /usr/bin/time -v $PIN "$BIN/bench_original_spq" 500 > "$RAW/orig_spq_$k.txt" 2> "$RAW/orig_spq_$k.time"
done
for k in 1 2 3; do
    /usr/bin/time -v $PIN "$BIN/bench_original_nay" 500 > "$RAW/orig_nay_$k.txt" 2> "$RAW/orig_nay_$k.time"
done
$PIN "$BIN/adder" 8 173 94 > "$RAW/adder.txt"
$PIN "$BIN/adder" 16 40000 30000 >> "$RAW/adder.txt"
$PIN "$BIN/adder" 24 9999999 7777777 >> "$RAW/adder.txt"
$PIN "$BIN/test04_fft" > "$RAW/test04_fft.txt"
$PIN "$BIN/test09_bootstrap" > "$RAW/test09_bootstrap.txt"
$PIN "$BIN/profile" 20000 > "$RAW/profile.txt" 2> "$RAW/profile.log"
$PIN "$BIN/noise" 10 1000 "$RAW/noise" 2> "$RAW/noise.log"
lscpu | grep "Model name" | sed 's/.*: *//' > "$RAW/cpu.txt"
g++ --version | head -1 > "$RAW/compiler.txt"

python3 report/experiments/summarize.py
