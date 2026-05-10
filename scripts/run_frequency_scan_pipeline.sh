#!/usr/bin/env bash
set -euo pipefail

CONFIG=${1:-configs/base.yaml}
F_MIN_HZ=${2:-20000}
F_MAX_HZ=${3:-200000}
STEPS=${4:-25}
MAX_CANDIDATES=${5:-10}
BUILD_DIR=${BUILD_DIR:-build}
OUTPUT_DIR=${OUTPUT_DIR:-results/frequency_scan_${F_MIN_HZ}_${F_MAX_HZ}_${STEPS}}

cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_DIR" -j --target build_geometry scan_frequency solve_bubble_candidates

mkdir -p "$OUTPUT_DIR"
"$BUILD_DIR/build_geometry" "$CONFIG" "$OUTPUT_DIR/mesh.msh"
"$BUILD_DIR/scan_frequency" "$CONFIG" "$OUTPUT_DIR/mesh.msh" "$OUTPUT_DIR" "$F_MIN_HZ" "$F_MAX_HZ" "$STEPS" "$MAX_CANDIDATES" 0
"$BUILD_DIR/solve_bubble_candidates" "$CONFIG" "$OUTPUT_DIR/bubble_excitations_all.csv" "$OUTPUT_DIR/bubble_results_all.csv" "$OUTPUT_DIR/timeseries"
python3 scripts/plot_frequency_scan.py "$OUTPUT_DIR" "$OUTPUT_DIR/bubble_results_all.csv"

printf '\nDone. Results: %s\n' "$OUTPUT_DIR"
