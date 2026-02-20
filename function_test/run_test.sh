#!/bin/bash
# Run both C++ and R score tests and compare results
set -e

DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"

echo "==================================================="
echo "  Score Test Comparison: C++ vs R"
echo "==================================================="
echo ""

# --- Step 1: Build C++ test if needed ---
if [ ! -f "$DIR/test_scoretest" ] || [ "$DIR/test_scoretest.cpp" -nt "$DIR/test_scoretest" ]; then
    echo "Building C++ test..."
    make -C "$DIR" clean && make -C "$DIR"
    echo ""
fi

# --- Step 2: Run C++ test ---
echo "Running C++ score test..."
"$DIR/test_scoretest"
echo ""

# --- Step 3: Run R test ---
echo "Running R score test..."
~/.pixi/bin/pixi run --manifest-path=/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/code_copy/SAIGE_isolated/pixi.toml Rscript "$DIR/test_scoretest.R"
echo ""

# --- Step 4: Run comparison ---
echo "Running comparison..."
~/.pixi/bin/pixi run --manifest-path=/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/code_copy/SAIGE_isolated/pixi.toml Rscript "$DIR/compare_results.R"
echo ""

echo "==================================================="
echo "  Done. Results in $DIR/RESULTS_scoretest.txt"
echo "==================================================="
