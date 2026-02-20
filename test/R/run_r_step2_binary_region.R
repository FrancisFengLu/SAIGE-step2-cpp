#!/usr/bin/env Rscript
# ==============================================================================
# SAIGE Step 2 -- Region/Gene-based Testing (Binary Trait)
# ==============================================================================
#
# Purpose:
#   Run SAIGE Step 2 region-based association (BURDEN/SKAT/SKAT-O) via
#   SPAGMMATtest() on the binary example data and save reference output.
#   This exercises the SPA Phi adjustment code path for binary traits
#   in region testing.
#
# Usage:
#   Rscript run_r_step2_binary_region.R
#
# Output:
#   r_step2_binary_region_results.txt          -- Region-level p-values
#   r_step2_binary_region_results.txt.singleAssoc.txt -- Single-variant within groups
#
# ==============================================================================

options(digits = 15)

# --- Timing -------------------------------------------------------------------
time_start <- proc.time()

# --- Paths --------------------------------------------------------------------
script_dir <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) {
  getwd()
})

checkpoint_dir <- normalizePath(script_dir, mustWork = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

# Input paths (absolute)
base_dir    <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi"
model_file  <- file.path(base_dir, "SAIGE/extdata/output/example_binary.rda")
vr_file     <- file.path(base_dir, "SAIGE/extdata/output/example_binary.varianceRatio.txt")
plink_stem  <- file.path(base_dir, "SAIGE/extdata/input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly")
bed_file    <- paste0(plink_stem, ".bed")
bim_file    <- paste0(plink_stem, ".bim")
fam_file    <- paste0(plink_stem, ".fam")

# Group file for R (uses rsid format matching BIM column 2)
group_file  <- file.path(base_dir, "Step_2_Feb_11/test/data/group_region_test_R.txt")

# Output paths
output_file <- file.path(checkpoint_dir, "r_step2_binary_region_results.txt")

# --- Preflight checks ---------------------------------------------------------

cat("=== SAIGE Step 2: R Region-Based Test (Binary Trait) ===\n\n")
cat("Script directory :", checkpoint_dir, "\n")
cat("Model file       :", model_file, "\n")
cat("VR file          :", vr_file, "\n")
cat("PLINK stem       :", plink_stem, "\n")
cat("Group file       :", group_file, "\n")
cat("Output file      :", output_file, "\n\n")

# Verify input files exist
for (f in c(model_file, vr_file, bed_file, bim_file, fam_file, group_file)) {
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
# Run SAIGE Step 2: Region-based association test (binary trait)
# ==============================================================================

cat("--- Running SPAGMMATtest (binary region-based, LOCO=TRUE, SPA) ---\n\n")

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
    is_Firth_beta        = FALSE,
    pCutoffforFirth      = 0.01,
    is_output_moreDetails = FALSE,
    is_overwrite_output  = TRUE,
    # Region-specific parameters
    groupFile            = group_file,
    annotation_in_groupTest = c("lof", "missense;lof"),
    maxMAF_in_groupTest  = c(0.0001, 0.001, 0.01),
    r.corr               = 0,        # SKAT-O (includes BURDEN + SKAT)
    MACCutoff_to_CollapseUltraRare = 10,
    markers_per_chunk_in_groupTest = 100,
    groups_per_chunk     = 100,
    weights.beta         = c(1, 25),
    is_single_in_groupTest = TRUE,
    is_no_weight_in_groupTest = FALSE,
    is_output_markerList_in_groupTest = TRUE
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
  cat("Total region-annotation-MAF rows:", nrow(results), "\n")
  cat("Output columns:", paste(names(results), collapse = ", "), "\n\n")
  cat("Results:\n")
  print(results)
} else {
  cat("WARNING: Output file not found at", output_file, "\n")
}

# Check single-variant in group output
single_file <- paste0(output_file, ".singleAssoc.txt")
if (file.exists(single_file)) {
  single_results <- read.table(single_file, header = TRUE, sep = " ",
                               stringsAsFactors = FALSE)
  cat("\nSingle-variant within groups:", nrow(single_results), "markers\n")
  cat("Single columns:", paste(names(single_results), collapse = ", "), "\n")
  cat("First 5 single-variant results:\n")
  print(head(single_results, 5))
} else {
  cat("\nWARNING: Single-in-group output file not found at", single_file, "\n")
}

# Check marker list output
marker_file <- paste0(output_file, ".markerList.txt")
if (file.exists(marker_file)) {
  cat("\nMarker list file exists:", marker_file, "\n")
  marker_lines <- readLines(marker_file)
  cat("  Lines:", length(marker_lines), "\n")
  if (length(marker_lines) > 0) {
    cat("  First 3 lines:\n")
    for (l in head(marker_lines, 3)) cat("    ", l, "\n")
  }
} else {
  cat("\nWARNING: Marker list output file not found at", marker_file, "\n")
}

# Timing
elapsed_total <- (proc.time() - time_start)[["elapsed"]]
elapsed_step2 <- (time_step2_end - time_step2_start)[["elapsed"]]
cat("\nTiming:\n")
cat("  SPAGMMATtest :", sprintf("%.2f", elapsed_step2), "seconds\n")
cat("  Total script :", sprintf("%.2f", elapsed_total), "seconds\n")

cat("\n=== Done ===\n")
