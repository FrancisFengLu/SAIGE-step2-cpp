#!/usr/bin/env Rscript
# ==============================================================================
# Test 8: Multi-Category Variance Ratio -- R Reference
# ==============================================================================
#
# Tests that R and C++ produce identical results when multiple VR categories
# are used (different VR values for different MAC ranges).
#
# VR file: multiVR_varianceRatio.txt (2 "null" rows)
#   VR[1] = 1.15000000000000 for MAC in (10.5, 20.5]  (low-MAC category)
#   VR[2] = 1.22420115246319 for MAC > 20.5             (high-MAC category)
#
# MAC category boundaries (R defaults):
#   cateVarRatioMinMACVecExclude = c(10.5, 20.5)
#   cateVarRatioMaxMACVecInclude = c(20.5)   -- R appends N=1000
#
# Usage:
#   Rscript run_r_step2_multiVR.R
#
# ==============================================================================

options(digits = 15)

time_start <- proc.time()

# --- Paths ---
base_dir    <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi"
test_dir    <- file.path(base_dir, "Step_2_Feb_11/test")
model_file  <- file.path(base_dir, "SAIGE/extdata/output/example_quantitative.rda")
vr_file     <- file.path(test_dir, "data/multiVR_varianceRatio.txt")
plink_stem  <- file.path(base_dir, "SAIGE/extdata/input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly")
bed_file    <- paste0(plink_stem, ".bed")
bim_file    <- paste0(plink_stem, ".bim")
fam_file    <- paste0(plink_stem, ".fam")
output_file <- file.path(test_dir, "output/r_multiVR_results.txt")

cat("=== Test 8: Multi-Category VR -- R Reference ===\n\n")
cat("Model file  :", model_file, "\n")
cat("VR file     :", vr_file, "\n")
cat("Output file :", output_file, "\n\n")

# Verify input files
for (f in c(model_file, vr_file, bed_file, bim_file, fam_file)) {
  if (!file.exists(f)) {
    stop("Required file not found: ", f)
  }
}
cat("All input files verified.\n\n")

# Show VR file contents
cat("VR file contents:\n")
vr_lines <- readLines(vr_file)
for (line in vr_lines) {
  cat("  ", line, "\n")
}
cat("\n")

# --- Load SAIGE ---
library(SAIGE)
cat("SAIGE version:", as.character(packageVersion("SAIGE")), "\n\n")

# --- Run SPAGMMATtest with multi-category VR ---
cat("--- Running SPAGMMATtest with multi-category VR ---\n")
cat("  cateVarRatioMinMACVecExclude = c(10.5, 20.5)\n")
cat("  cateVarRatioMaxMACVecInclude = c(20.5)  [R will append N=1000]\n\n")

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
    cateVarRatioMinMACVecExclude = c(10.5, 20.5),
    cateVarRatioMaxMACVecInclude = c(20.5)
  )
}, error = function(e) {
  cat("ERROR in SPAGMMATtest:", conditionMessage(e), "\n")
  traceback()
})

# --- Summary ---
cat("\n--- Summary ---\n\n")

if (file.exists(output_file)) {
  results <- read.table(output_file, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
  cat("Total markers tested:", nrow(results), "\n")
  cat("Output columns:", paste(names(results), collapse = ", "), "\n\n")

  # Show first 5 rows
  cat("First 5 results:\n")
  print(head(results, 5))

  # Analyze MAC distribution to verify multi-VR selection
  if ("AC_Allele2" %in% names(results)) {
    mac <- pmin(results$AC_Allele2, 2 * results$N - results$AC_Allele2)
    cat("\n\nMAC distribution across VR categories:\n")
    cat("  MAC <= 10.5           :", sum(mac <= 10.5), "markers (fallback to VR[1] = 1.15)\n")
    cat("  MAC in (10.5, 20.5]   :", sum(mac > 10.5 & mac <= 20.5), "markers (VR[1] = 1.15)\n")
    cat("  MAC > 20.5            :", sum(mac > 20.5), "markers (VR[2] = 1.22420)\n")
  }
} else {
  cat("WARNING: Output file not found at", output_file, "\n")
}

elapsed <- (proc.time() - time_start)[["elapsed"]]
cat("\nTotal time:", sprintf("%.2f", elapsed), "seconds\n")
cat("\n=== Done ===\n")
