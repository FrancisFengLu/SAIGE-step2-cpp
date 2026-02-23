#!/usr/bin/env Rscript
# ==============================================================================
# SAIGE Step 2 -- Run R Reference with Conditional Analysis
# ==============================================================================
#
# Purpose:
#   Run SAIGE Step 2 single-variant association with conditional analysis
#   (conditioning on specified markers). Produces reference output for
#   comparison with the standalone C++ port.
#
# Usage:
#   Rscript run_r_step2_conditional.R
#
# Output:
#   r_step2_conditional_results.txt -- SPAGMMATtest results with condition columns
#
# ==============================================================================

options(digits = 15)

# --- Timing -------------------------------------------------------------------
time_start <- proc.time()

# --- Paths --------------------------------------------------------------------

# Resolve this script's directory
script_dir <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) {
  getwd()
})

checkpoint_dir <- normalizePath(script_dir, mustWork = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

# Input paths (absolute)
base_dir    <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi"
model_file  <- file.path(base_dir, "SAIGE/extdata/output/example_quantitative.rda")
vr_file     <- file.path(base_dir, "SAIGE/extdata/output/example_quantitative.varianceRatio.txt")
plink_stem  <- file.path(base_dir, "SAIGE/extdata/input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly")
bed_file    <- paste0(plink_stem, ".bed")
bim_file    <- paste0(plink_stem, ".bim")
fam_file    <- paste0(plink_stem, ".fam")

# Output paths
output_file <- file.path(checkpoint_dir, "r_step2_conditional_results.txt")

# Conditioning markers: use rsIDs from the .bim file
# These are common markers with good allele frequencies:
#   rs1e+05 (AF=0.239), rs8 (AF=0.3795), rs9 (AF=0.186)
condition_string <- "rs1e+05,rs8,rs9"

# --- Preflight checks ---------------------------------------------------------

cat("=== SAIGE Step 2: R Reference Run (Conditional Analysis) ===\n\n")
cat("Script directory :", checkpoint_dir, "\n")
cat("Model file       :", model_file, "\n")
cat("VR file          :", vr_file, "\n")
cat("PLINK stem       :", plink_stem, "\n")
cat("Output file      :", output_file, "\n")
cat("Condition        :", condition_string, "\n\n")

# Verify input files exist
for (f in c(model_file, vr_file, bed_file, bim_file, fam_file)) {
  if (!file.exists(f)) {
    stop("Required file not found: ", f)
  }
}
cat("All input files verified.\n\n")

# --- Load SAIGE ---------------------------------------------------------------

if (!requireNamespace("SAIGE", quietly = TRUE)) {
  stop("SAIGE R package is not installed. Please install it first.")
}
library(SAIGE)
cat("SAIGE version:", as.character(packageVersion("SAIGE")), "\n\n")

# ==============================================================================
# Run SAIGE Step 2 with conditional analysis
# ==============================================================================

cat("--- Running SPAGMMATtest with condition ---\n\n")

time_step2_start <- proc.time()

tryCatch({
  SPAGMMATtest(
    bedFile              = bed_file,
    bimFile              = bim_file,
    famFile              = fam_file,
    AlleleOrder          = "alt-first",
    chrom                = "1",
    GMMATmodelFile       = model_file,
    varianceRatioFile    = vr_file,
    SAIGEOutputFile      = output_file,
    min_MAF              = 0,
    min_MAC              = 0.5,
    LOCO                 = TRUE,
    is_Firth_beta        = TRUE,
    pCutoffforFirth      = 0.01,
    is_output_moreDetails = TRUE,
    is_overwrite_output  = TRUE,
    condition            = condition_string
  )
}, error = function(e) {
  cat("ERROR in SPAGMMATtest:", conditionMessage(e), "\n")
  cat("Traceback:\n")
  traceback()
})

time_step2_end <- proc.time()

# ==============================================================================
# Summary
# ==============================================================================

cat("\n--- Summary ---\n\n")

if (file.exists(output_file)) {
  results <- read.table(output_file, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
  cat("Total markers tested:", nrow(results), "\n")

  # Show columns
  cat("Output columns:", paste(names(results), collapse = ", "), "\n\n")

  # Show first 5 rows
  cat("First 5 results:\n")
  print(head(results, 5))

  # Check for conditional columns
  cond_cols <- grep("_c$", names(results), value = TRUE)
  cat("\nConditional columns found:", paste(cond_cols, collapse = ", "), "\n")

  # Count significant hits
  if ("p.value" %in% names(results)) {
    pvals <- as.numeric(results[["p.value"]])
    pvals <- pvals[!is.na(pvals)]
    cat("\nSignificant hits (unconditional):\n")
    cat("  p < 0.05  :", sum(pvals < 0.05), "\n")
    cat("  p < 0.01  :", sum(pvals < 0.01), "\n")
    cat("  p < 1e-3  :", sum(pvals < 1e-3), "\n")
  }

  if ("p.value_c" %in% names(results)) {
    pvals_c <- as.numeric(results[["p.value_c"]])
    pvals_c <- pvals_c[!is.na(pvals_c)]
    cat("\nSignificant hits (conditional):\n")
    cat("  p < 0.05  :", sum(pvals_c < 0.05), "\n")
    cat("  p < 0.01  :", sum(pvals_c < 0.01), "\n")
    cat("  p < 1e-3  :", sum(pvals_c < 1e-3), "\n")
  }
} else {
  cat("WARNING: Output file not found at", output_file, "\n")
}

# Timing
elapsed_total <- (proc.time() - time_start)[["elapsed"]]
elapsed_step2 <- (time_step2_end - time_step2_start)[["elapsed"]]
cat("\nTiming:\n")
cat("  SPAGMMATtest :", sprintf("%.2f", elapsed_step2), "seconds\n")
cat("  Total script :", sprintf("%.2f", elapsed_total), "seconds\n")

cat("\n=== Done ===\n")
