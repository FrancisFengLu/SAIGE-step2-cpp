#!/bin/bash
# ==============================================================================
# SAIGE Step 2: End-to-End R vs C++ Comparison
# ==============================================================================
# Steps:
#   1. Convert .rda null model to .arma format
#   2. Run R SAIGE Step 2 (via pixi)
#   3. Run C++ saige-step2 (with converted model)
#   4. Compare results
#
# Usage: cd Step_2_Feb_11 && bash test/run_comparison.sh
# ==============================================================================

set -e

# Resolve paths
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CPP_DIR="$BASE_DIR/code_copy/cpp_standalone"

PIXI="$HOME/.pixi/bin/pixi run --manifest-path=/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/code_copy/SAIGE_isolated/pixi.toml"

echo "============================================"
echo "  SAIGE Step 2: R vs C++ Comparison"
echo "============================================"
echo ""
echo "Base dir: $BASE_DIR"
echo ""

# --- Step 1: Convert .rda to .arma -------------------------------------------
echo "--- Step 1: Converting .rda null model to .arma format ---"
mkdir -p "$BASE_DIR/test/data/nullmodel_from_rda"
$PIXI Rscript "$BASE_DIR/test/R/convert_rda_to_arma.R" 2>&1
echo ""

# Verify conversion produced files
if [ ! -f "$BASE_DIR/test/data/nullmodel_from_rda/nullmodel.json" ]; then
    echo "ERROR: Conversion failed — nullmodel.json not found"
    exit 1
fi
echo "Conversion successful."
echo ""

# --- Step 2: Run R SAIGE Step 2 ----------------------------------------------
echo "--- Step 2: Running R SAIGE Step 2 ---"
$PIXI Rscript "$BASE_DIR/test/R/run_r_step2.R" 2>&1
echo ""

# Verify R output
R_RESULTS="$BASE_DIR/test/R/r_step2_results.txt"
if [ ! -f "$R_RESULTS" ]; then
    echo "ERROR: R Step 2 failed — output file not found"
    exit 1
fi
R_LINES=$(wc -l < "$R_RESULTS")
echo "R results: $R_LINES lines"
echo ""

# --- Step 3: Run C++ saige-step2 ---------------------------------------------
echo "--- Step 3: Running C++ saige-step2 ---"
mkdir -p "$BASE_DIR/test/output"
mkdir -p "$BASE_DIR/test/Cpp_compare"

# Build if needed
if [ ! -f "$CPP_DIR/saige-step2" ]; then
    echo "Building saige-step2..."
    make -C "$CPP_DIR" -j4
fi

"$CPP_DIR/saige-step2" "$CPP_DIR/config_compare.yaml" 2>&1
echo ""

# Verify C++ output
CPP_RESULTS="$BASE_DIR/test/output/cpp_compare_results.txt"
if [ ! -f "$CPP_RESULTS" ]; then
    echo "ERROR: C++ Step 2 failed — output file not found"
    exit 1
fi
CPP_LINES=$(wc -l < "$CPP_RESULTS")
echo "C++ results: $CPP_LINES lines"
echo ""

# --- Step 4: Compare ---------------------------------------------------------
echo "--- Step 4: Comparing R vs C++ results ---"
$PIXI Rscript "$BASE_DIR/test/compare_step2.R" \
    "$R_RESULTS" \
    "$CPP_RESULTS" \
    "$BASE_DIR/test/RESULTS_step2.txt" 2>&1

echo ""
echo "============================================"
echo "  Done. See test/RESULTS_step2.txt"
echo "============================================"
