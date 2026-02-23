#!/usr/bin/env Rscript
# ==============================================================================
# SAIGE Step 2 -- Run R Reference for sparseGRM model with is_noadjCov = TRUE
# ==============================================================================
#
# Purpose:
#   Run SAIGE Step 2 single-variant association via SPAGMMATtest()
#   using the quantitative sparse GRM example model with is_noadjCov = TRUE.
#   This exercises the scoreTestFast_noadjCov code path, which skips
#   covariate adjustment in the fast-pass score test.
#
# Model: example_quantitative_sparseGRM.rda (LOCO=FALSE)
# Input: PLINK files (nfam_100_nindep_0_step1_includeMoreRareVariants_poly)
# sparseGRM: sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx
#
# Usage:
#   Rscript run_r_step2_sparse_noadjcov.R
#   OR via pixi:
#   ~/.pixi/bin/pixi run --manifest-path=.../pixi.toml Rscript run_r_step2_sparse_noadjcov.R
#
# Output:
#   r_step2_sparse_noadjcov_results.txt -- SPAGMMATtest single-variant results
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

# Input paths (absolute)
base_dir    <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi"
model_file  <- file.path(base_dir, "SAIGE/extdata/output/example_quantitative_sparseGRM.rda")
vr_file     <- file.path(base_dir, "Step_2_Feb_11/test/data/varianceRatio_sparse_noadjcov.txt")
plink_stem  <- file.path(base_dir, "SAIGE/extdata/input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly")
bed_file    <- paste0(plink_stem, ".bed")
bim_file    <- paste0(plink_stem, ".bim")
fam_file    <- paste0(plink_stem, ".fam")

# Sparse GRM files
sparseGRM_file <- file.path(base_dir, "SAIGE/extdata/output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
sparseGRM_sampleID_file <- file.path(base_dir, "SAIGE/extdata/output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")

# Output path
output_dir  <- file.path(base_dir, "Step_2_Feb_11/test/R")
output_file <- file.path(output_dir, "r_step2_sparse_noadjcov_results.txt")

# --- Preflight checks ---------------------------------------------------------

cat("=== SAIGE Step 2: R Reference Run (Sparse GRM, Quantitative, noadjCov=TRUE) ===\n\n")
cat("Model file       :", model_file, "\n")
cat("VR file          :", vr_file, "\n")
cat("PLINK stem       :", plink_stem, "\n")
cat("sparseGRM file   :", sparseGRM_file, "\n")
cat("sparseGRM IDs    :", sparseGRM_sampleID_file, "\n")
cat("Output file      :", output_file, "\n")
cat("is_noadjCov      : TRUE\n\n")

# Verify input files exist
for (f in c(model_file, vr_file, bed_file, bim_file, fam_file,
            sparseGRM_file, sparseGRM_sampleID_file)) {
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
# Run SAIGE Step 2 (single-variant association test) with is_noadjCov = TRUE
# ==============================================================================

cat("--- Running SPAGMMATtest (sparse GRM, LOCO=FALSE, is_noadjCov=TRUE) ---\n\n")

time_step2_start <- proc.time()

tryCatch({
  SPAGMMATtest(
    bedFile              = bed_file,
    bimFile              = bim_file,
    famFile              = fam_file,
    AlleleOrder          = "alt-first",
    GMMATmodelFile       = model_file,
    varianceRatioFile    = vr_file,
    SAIGEOutputFile      = output_file,
    sparseGRMFile        = sparseGRM_file,
    sparseGRMSampleIDFile = sparseGRM_sampleID_file,
    min_MAF              = 0,
    min_MAC              = 0.5,
    LOCO                 = FALSE,
    is_fastTest          = TRUE,
    pval_cutoff_for_fastTest = 0.05,
    is_output_moreDetails = TRUE,
    is_overwrite_output  = TRUE,
    is_noadjCov          = TRUE
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

  # Count significant hits at different thresholds
  if ("p.value" %in% names(results)) {
    pval_col <- "p.value"
  } else if ("pval" %in% names(results)) {
    pval_col <- "pval"
  } else {
    pval_col <- NULL
  }

  if (!is.null(pval_col)) {
    pvals <- as.numeric(results[[pval_col]])
    pvals <- pvals[!is.na(pvals)]
    cat("\nSignificant hits:\n")
    cat("  p < 0.05  :", sum(pvals < 0.05), "\n")
    cat("  p < 0.01  :", sum(pvals < 0.01), "\n")
    cat("  p < 1e-3  :", sum(pvals < 1e-3), "\n")
    cat("  p < 5e-8  :", sum(pvals < 5e-8), "\n")
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
