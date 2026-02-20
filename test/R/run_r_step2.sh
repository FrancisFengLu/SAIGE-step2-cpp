#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PIXI="$HOME/.pixi/bin/pixi run --manifest-path=/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/code_copy/SAIGE_isolated/pixi.toml"

echo "=== Running R SAIGE Step 2 ==="
echo "Script directory: $SCRIPT_DIR"
echo "Pixi command: $PIXI"
echo ""

$PIXI Rscript "$SCRIPT_DIR/run_r_step2.R" 2>&1 | tee "$SCRIPT_DIR/r_step2_log.txt"

echo ""
echo "=== Done ==="
echo "Log saved to: $SCRIPT_DIR/r_step2_log.txt"
echo "Results:      $SCRIPT_DIR/r_step2_results.txt"
echo "Checkpoints:  $SCRIPT_DIR/ckpt_*.txt"
