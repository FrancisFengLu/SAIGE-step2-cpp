#!/usr/bin/env Rscript
# ==============================================================================
# SAIGE Step 2 -- Run R Reference for Binary Trait
# ==============================================================================
#
# Purpose:
#   Run SAIGE Step 2 single-variant association via SPAGMMATtest() on the
#   binary example data (example_binary.rda) and save reference output.
#   This exercises the SPA and Firth correction code paths.
#
# Usage:
#   Rscript run_r_step2_binary.R
#
# Output:
#   test/R/r_step2_binary_results.txt  -- SPAGMMATtest single-variant results
#
# ==============================================================================

options(digits = 15)

# --- Timing -------------------------------------------------------------------
time_start <- proc.time()

# --- Paths --------------------------------------------------------------------

# Input paths (absolute)
base_dir    <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi"
model_file  <- file.path(base_dir, "SAIGE/extdata/output/example_binary.rda")
vr_file     <- file.path(base_dir, "SAIGE/extdata/output/example_binary.varianceRatio.txt")

# genotype_100markers PLINK files (all chr 1, 10000 samples, 100 markers)
plink_stem  <- file.path(base_dir, "SAIGE/extdata/input/genotype_100markers")
bed_file    <- paste0(plink_stem, ".bed")
bim_file    <- paste0(plink_stem, ".bim")
fam_file    <- paste0(plink_stem, ".fam")

# Output paths
output_dir  <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Step_2_Feb_11/test/R"
output_file <- file.path(output_dir, "r_step2_binary_results.txt")

# --- Preflight checks ---------------------------------------------------------

cat("=== SAIGE Step 2: R Reference Run (Binary Trait) ===\n\n")
cat("Model file       :", model_file, "\n")
cat("VR file          :", vr_file, "\n")
cat("PLINK stem       :", plink_stem, "\n")
cat("Output file      :", output_file, "\n\n")

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

# --- Quick model inspection ---------------------------------------------------

cat("--- Inspecting binary null model ---\n\n")

env <- new.env()
load(model_file, envir = env)
modglmm <- env$modglmm

cat("traitType:", modglmm$traitType, "\n")
cat("LOCO:", modglmm$LOCO, "\n")
cat("theta:", paste(modglmm$theta, collapse = ", "), "\n")
cat("N:", length(modglmm$y), "\n")
cat("X dims:", nrow(modglmm$X), "x", ncol(modglmm$X), "\n")
cat("sampleID count:", length(modglmm$sampleID), "\n")
cat("y summary: 0s =", sum(modglmm$y == 0), ", 1s =", sum(modglmm$y == 1), "\n\n")

# Show variance ratio file
vr_raw <- readLines(vr_file)
cat("Variance ratio file contents:\n")
for (line in vr_raw) cat("  ", line, "\n")
cat("\n")

rm(env, modglmm)
gc(verbose = FALSE)

# ==============================================================================
# Run SAIGE Step 2 (single-variant, binary trait, with SPA + Firth)
# ==============================================================================

cat("--- Running SPAGMMATtest (binary, LOCO=TRUE, SPA, Firth) ---\n\n")

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
    is_overwrite_output  = TRUE
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

  # Check for SPA-adjusted results
  if ("p.value.NA" %in% names(results)) {
    cat("\nSPA-adjusted markers (where p.value.NA is not NA):\n")
    spa_markers <- sum(!is.na(results$p.value.NA))
    cat("  Markers with SPA adjustment:", spa_markers, "\n")
  }

  # Check for Firth results
  if ("Is.SPA.converge" %in% names(results)) {
    cat("  SPA converged:", sum(results$Is.SPA.converge == 1, na.rm=TRUE), "\n")
    cat("  SPA did not converge:", sum(results$Is.SPA.converge == 0, na.rm=TRUE), "\n")
  }

  # Show top 10 by p-value
  cat("\nTop 10 by p-value:\n")
  top <- head(results[order(as.numeric(results[[pval_col]])), ], 10)
  cols_to_show <- intersect(names(results),
    c("CHR", "POS", "MarkerID", "Allele1", "Allele2", "AF_Allele2", "BETA", "SE", "p.value"))
  print(top[, cols_to_show])

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
