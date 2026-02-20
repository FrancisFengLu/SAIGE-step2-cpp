#!/usr/bin/env Rscript
# ==============================================================================
# SAIGE Step 2 -- Run R Reference + Checkpoint Extraction
# ==============================================================================
#
# Purpose:
#   1. Run SAIGE Step 2 single-variant association via SPAGMMATtest()
#      on the quantitative example data and save reference output.
#   2. Extract checkpoint files (null model internals, variance ratios,
#      key intermediate values) for comparison with the standalone C++ port.
#
# Usage:
#   Rscript run_r_step2.R
#
# Output:
#   r_step2_results.txt          -- SPAGMMATtest single-variant results
#   ckpt_01_null_model.txt       -- null model scalars
#   ckpt_02_mu.txt               -- first 10 mu values (1-indexed)
#   ckpt_03_res.txt              -- first 10 residual values (1-indexed)
#   ckpt_04_XVX.txt              -- XVX matrix [p x p]
#   ckpt_05_variance_ratios.txt  -- variance ratio values
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
output_file <- file.path(checkpoint_dir, "r_step2_results.txt")

# --- Helper: write checkpoint with full precision -----------------------------

fmt <- function(x) {
  if (is.numeric(x)) {
    formatC(x, digits = 15, format = "g")
  } else {
    as.character(x)
  }
}

write_ckpt <- function(filename, lines) {
  filepath <- file.path(checkpoint_dir, filename)
  writeLines(lines, filepath)
  cat("  Wrote checkpoint:", filepath, "\n")
}

# --- Preflight checks ---------------------------------------------------------

cat("=== SAIGE Step 2: R Reference Run ===\n\n")
cat("Script directory :", checkpoint_dir, "\n")
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

# ==============================================================================
# PART 1: Extract checkpoints from the null model BEFORE running Step 2
# ==============================================================================

cat("--- Extracting null model checkpoints ---\n\n")

# Load the .rda file
env <- new.env()
load(model_file, envir = env)
model_names <- ls(env)
cat("Objects in .rda:", paste(model_names, collapse = ", "), "\n")

# The model object is typically named 'modglmm'
modglmm <- NULL
for (nm in model_names) {
  obj <- get(nm, envir = env)
  if (is.list(obj) && ("theta" %in% names(obj) || "tau" %in% names(obj))) {
    modglmm <- obj
    cat("Using model object: '", nm, "'\n")
    break
  }
}

if (is.null(modglmm)) {
  stop("Could not find null model object in .rda file")
}

# --- Replicate what ReadModel() does for LOCO=TRUE, quantitative trait --------
#
# From readInGLMM.R: When LOCO=TRUE and model has LOCO, SAIGE uses
# per-chromosome fitted.values and residuals from LOCOResult.
# For chr 1 markers, it uses LOCOResult[[1]].

# Check LOCO status
cat("modglmm$LOCO =", modglmm$LOCO, "\n")

# For LOCO=TRUE: use chr 1 LOCO values (matches what convert_rda_to_arma.R extracts)
if (!is.null(modglmm$LOCO) && modglmm$LOCO &&
    !is.null(modglmm$LOCOResult) && !is.null(modglmm$LOCOResult[[1]])) {
  cat("Using chr 1 LOCO model values\n")
  mu  <- as.vector(modglmm$LOCOResult[[1]]$fitted.values)
  res <- as.vector(modglmm$LOCOResult[[1]]$residuals)
} else {
  cat("No LOCO results, using base model values\n")
  mu  <- as.vector(modglmm$fitted.values)
  res <- as.vector(modglmm$residuals)
}
y   <- modglmm$y
X   <- modglmm$X
tau <- modglmm$theta
N   <- length(mu)
p   <- ncol(X)

traitType <- modglmm$traitType
cat("traitType =", traitType, "\n")
cat("N =", N, "\n")
cat("p =", p, "\n")
cat("tau =", paste(fmt(tau), collapse = ", "), "\n\n")

# Compute mu2 for quantitative trait (from readInGLMM.R line 128)
if (traitType == "quantitative") {
  mu2 <- (1 / tau[1]) * rep(1, N)
} else if (traitType == "binary") {
  mu2 <- mu * (1 - mu)
} else {
  stop("Unsupported traitType: ", traitType)
}

# Compute obj.noK matrices using ScoreTest_NULL_Model logic
# (from SAIGE_fitGLMM_fast.R lines 802-814)
V  <- as.vector(mu2)
XV <- t(X * V)
XVX <- t(X) %*% t(XV)
XVX_inv <- solve(XVX)
XXVX_inv <- X %*% XVX_inv
XVX_inv_XV <- XXVX_inv * V
S_a <- colSums(X * res)

cat("XVX dimensions:", nrow(XVX), "x", ncol(XVX), "\n")
cat("S_a:", paste(fmt(S_a), collapse = ", "), "\n\n")

# But also check if obj.noK is already stored in the model (from Step 1)
if (!is.null(modglmm$obj.noK)) {
  cat("obj.noK found in model - comparing with recomputed values\n")
  cat("  obj.noK$XVX[1,1] =", fmt(modglmm$obj.noK$XVX[1,1]), "\n")
  cat("  recomputed XVX[1,1] =", fmt(XVX[1,1]), "\n")
  # Use the stored obj.noK values (they include LOCO adjustments if any)
  XVX_stored <- modglmm$obj.noK$XVX
  if (!is.null(XVX_stored)) {
    cat("  Using stored obj.noK$XVX for checkpoint\n")
    XVX <- XVX_stored
  }
  if (!is.null(modglmm$obj.noK$S_a)) {
    S_a <- modglmm$obj.noK$S_a
    cat("  Using stored obj.noK$S_a:", paste(fmt(S_a), collapse = ", "), "\n")
  }
  cat("\n")
}

# --- Checkpoint 01: null model scalars ----------------------------------------

ckpt01_lines <- c(
  "parameter\tvalue",
  paste0("n\t", N),
  paste0("p\t", p),
  paste0("traitType\t", traitType),
  paste0("tau0\t", fmt(tau[1])),
  paste0("tau1\t", if (length(tau) > 1) fmt(tau[2]) else "NA"),
  paste0("SPA_Cutoff\t", 2)
)
write_ckpt("ckpt_01_null_model.txt", ckpt01_lines)

# --- Checkpoint 02: first 10 mu values (1-indexed) ---------------------------

n_show <- min(10, N)
ckpt02_lines <- c("index\tmu")
for (i in seq_len(n_show)) {
  ckpt02_lines <- c(ckpt02_lines, paste0(i, "\t", fmt(mu[i])))
}
write_ckpt("ckpt_02_mu.txt", ckpt02_lines)

# --- Checkpoint 03: first 10 residual values (1-indexed) ---------------------

ckpt03_lines <- c("index\tres")
for (i in seq_len(n_show)) {
  ckpt03_lines <- c(ckpt03_lines, paste0(i, "\t", fmt(res[i])))
}
write_ckpt("ckpt_03_res.txt", ckpt03_lines)

# --- Checkpoint 04: XVX matrix -----------------------------------------------

ckpt04_lines <- c(paste0("[", nrow(XVX), " x ", ncol(XVX), "]"))
for (i in seq_len(nrow(XVX))) {
  row_vals <- sapply(XVX[i, ], fmt)
  ckpt04_lines <- c(ckpt04_lines, paste(row_vals, collapse = "\t"))
}
write_ckpt("ckpt_04_XVX.txt", ckpt04_lines)

# --- Checkpoint 05: variance ratios -------------------------------------------

# Read the VR file. The format may vary by SAIGE version:
# Newer (>= 1.0.6): 3 columns (value, type, MAC_category)
# Older: 1 column (just the ratio values)
vr_raw <- readLines(vr_file)
cat("Raw VR file contents:\n")
for (line in vr_raw) {
  cat("  ", line, "\n")
}

# Parse variance ratio data
vr_data <- read.table(vr_file, header = FALSE, stringsAsFactors = FALSE)
cat("VR data: ", nrow(vr_data), "rows x", ncol(vr_data), "cols\n")

if (ncol(vr_data) == 3) {
  # Newer format: value, type ("null" or "sparse"), MAC_category
  ckpt05_lines <- c("index\tVR_value\tVR_type\tMAC_category")
  for (i in seq_len(nrow(vr_data))) {
    ckpt05_lines <- c(ckpt05_lines, paste0(
      i, "\t", fmt(as.numeric(vr_data[i, 1])), "\t",
      vr_data[i, 2], "\t", vr_data[i, 3]
    ))
  }
} else if (ncol(vr_data) == 1) {
  # Older format: just values
  ckpt05_lines <- c("index\tVR_value")
  for (i in seq_len(nrow(vr_data))) {
    ckpt05_lines <- c(ckpt05_lines, paste0(
      i, "\t", fmt(as.numeric(vr_data[i, 1]))
    ))
  }
} else {
  # Fallback: dump raw
  ckpt05_lines <- c(paste0("# columns: ", ncol(vr_data)))
  for (i in seq_len(nrow(vr_data))) {
    ckpt05_lines <- c(ckpt05_lines, paste(vr_data[i, ], collapse = "\t"))
  }
}
write_ckpt("ckpt_05_variance_ratios.txt", ckpt05_lines)

# Clean up model environment
rm(env, modglmm)
gc()

# ==============================================================================
# PART 2: Run SAIGE Step 2 (single-variant association test)
# ==============================================================================

cat("\n--- Running SPAGMMATtest ---\n\n")

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
# PART 3: Summary
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

# List all checkpoint files
cat("\nCheckpoint files:\n")
ckpt_files <- list.files(checkpoint_dir, pattern = "^ckpt_.*\\.txt$")
for (f in ckpt_files) {
  fpath <- file.path(checkpoint_dir, f)
  fsize <- file.info(fpath)$size
  cat(sprintf("  %-40s  %d bytes\n", f, fsize))
}

cat("\n=== Done ===\n")
