#!/usr/bin/env Rscript
# ==============================================================================
# convert_rda_to_arma.R
# ==============================================================================
#
# Purpose:
#   Convert a SAIGE .rda null model to Armadillo binary files (.arma) and a
#   nullmodel.json, so that the standalone C++ Step 2 binary can load the
#   exact same null model that R SAIGE Step 2 uses.
#
# Usage:
#   Rscript convert_rda_to_arma.R
#
#   Or via pixi:
#   ~/.pixi/bin/pixi run --manifest-path=.../pixi.toml Rscript convert_rda_to_arma.R
#
# Input:
#   SAIGE extdata example_quantitative.rda
#
# Output:
#   Directory: test/data/nullmodel_from_rda/
#     mu.arma, res.arma, y.arma, V.arma, S_a.arma,
#     X.arma, XVX.arma, XVX_inv.arma, XXVX_inv.arma, XV.arma, XVX_inv_XV.arma,
#     nullmodel.json
#
# ==============================================================================

options(digits = 15)

cat("=== convert_rda_to_arma.R ===\n\n")

# --- Paths -------------------------------------------------------------------

rda_path <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/SAIGE/extdata/output/example_quantitative.rda"
output_dir <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Step_2_Feb_11/test/data/nullmodel_from_rda"

cat("Input .rda:  ", rda_path, "\n")
cat("Output dir:  ", output_dir, "\n\n")

if (!file.exists(rda_path)) {
  stop("ERROR: .rda file not found: ", rda_path)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Helper: write Armadillo binary format -----------------------------------

# Write a vector as ARMA_COL_BIN_FN008 format
write_arma_vec <- function(filepath, vec) {
  vec <- as.double(vec)
  n <- length(vec)
  con <- file(filepath, "wb")
  on.exit(close(con))
  writeLines(c("ARMA_MAT_BIN_FN008", paste(n, 1)), con)
  writeBin(vec, con, size = 8, endian = "little")
}

# Write a matrix as ARMA_MAT_BIN_FN008 format (column-major)
write_arma_mat <- function(filepath, mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  nr <- nrow(mat)
  nc <- ncol(mat)
  con <- file(filepath, "wb")
  on.exit(close(con))
  writeLines(c("ARMA_MAT_BIN_FN008", paste(nr, nc)), con)
  # Armadillo stores column-major, same as R's internal storage
  writeBin(as.double(mat), con, size = 8, endian = "little")
}

# --- Helper: read back Armadillo binary for verification ---------------------

load_arma_bin <- function(filepath) {
  con <- file(filepath, "rb")
  on.exit(close(con))
  header <- readLines(con, n = 1)
  dims_line <- readLines(con, n = 1)
  dims <- as.integer(strsplit(trimws(dims_line), "\\s+")[[1]])
  nrow <- dims[1]; ncol <- dims[2]
  n_elements <- nrow * ncol
  values <- readBin(con, what = "double", n = n_elements, size = 8, endian = "little")
  if (ncol == 1) return(as.numeric(values))
  matrix(values, nrow = nrow, ncol = ncol, byrow = FALSE)
}

# --- Step 1: Load the .rda file ----------------------------------------------

cat("--- Loading .rda file ---\n")
env <- new.env()
load(rda_path, envir = env)
model_objects <- ls(env)
cat("Objects in .rda: ", paste(model_objects, collapse = ", "), "\n")

# Find the main model object (usually named 'modglmm')
modglmm <- NULL
for (obj_name in model_objects) {
  obj <- get(obj_name, envir = env)
  if (is.list(obj) && any(c("theta", "fitted.values", "residuals", "y", "X") %in% names(obj))) {
    modglmm <- obj
    cat("Found model object: '", obj_name, "'\n")
    break
  }
}

if (is.null(modglmm)) {
  stop("ERROR: Could not find a SAIGE model object in the .rda file")
}

cat("\nTop-level fields: ", paste(names(modglmm), collapse = ", "), "\n\n")

# --- Step 2: Determine LOCO vs non-LOCO -------------------------------------

is_loco <- !is.null(modglmm$LOCO) && modglmm$LOCO
cat("LOCO: ", is_loco, "\n")

if (is_loco) {
  cat("Model is LOCO. Checking available chromosomes...\n")
  avail_chrs <- c()
  for (chr in 1:22) {
    if (!is.null(modglmm$LOCOResult[[chr]]) &&
        !is.null(modglmm$LOCOResult[[chr]]$fitted.values)) {
      avail_chrs <- c(avail_chrs, chr)
    }
  }
  cat("Available LOCO chromosomes: ", paste(avail_chrs, collapse = ", "), "\n")

  # Use chromosome 1 for testing
  target_chr <- 1
  if (!(target_chr %in% avail_chrs)) {
    if (length(avail_chrs) > 0) {
      target_chr <- avail_chrs[1]
      cat("WARNING: Chr 1 not available, using chr", target_chr, "instead\n")
    } else {
      stop("ERROR: No LOCO results found in model")
    }
  }
  cat("Extracting data for chromosome: ", target_chr, "\n\n")

  # Extract per-chromosome fitted values and residuals
  loco_result <- modglmm$LOCOResult[[target_chr]]
  mu <- as.vector(loco_result$fitted.values)
  res <- as.vector(loco_result$residuals)

  # obj.noK comes from the LOCO result
  if (!is.null(loco_result$obj.noK)) {
    obj.noK <- loco_result$obj.noK
    cat("Using obj.noK from LOCOResult[[", target_chr, "]]\n")
  } else {
    cat("WARNING: obj.noK not in LOCOResult, will recompute\n")
    obj.noK <- NULL
  }
} else {
  cat("Model is non-LOCO.\n\n")
  mu <- as.vector(modglmm$fitted.values)
  res <- as.vector(modglmm$residuals)
  obj.noK <- modglmm$obj.noK
}

# --- Step 3: Extract common fields ------------------------------------------

tau <- modglmm$theta
traitType <- modglmm$traitType
y <- as.vector(modglmm$y)
X <- as.matrix(modglmm$X)
sampleIDs <- modglmm$sampleID

N <- length(y)
p <- ncol(X)

cat("traitType:  ", traitType, "\n")
cat("N (samples):", N, "\n")
cat("p (covars): ", p, "\n")
cat("tau:        ", paste(formatC(tau, digits = 15, format = "g"), collapse = ", "), "\n")

# Compute mu2 and V based on trait type
if (traitType == "binary") {
  mu2 <- mu * (1 - mu)
  cat("Binary trait: mu2 = mu*(1-mu)\n")
} else if (traitType == "quantitative") {
  # tau[1] in R is the first element (1-indexed), corresponding to tau_0
  mu2 <- rep(1.0 / tau[1], N)
  cat("Quantitative trait: mu2 = 1/tau[1] = ", 1.0/tau[1], "\n")
} else if (traitType == "survival") {
  mu2 <- mu
  cat("Survival trait: mu2 = mu\n")
} else {
  stop("ERROR: Unknown traitType: ", traitType)
}

V <- mu2

# --- Step 4: Compute or extract obj.noK matrices ----------------------------

# ScoreTest_NULL_Model computes:
#   V = mu2
#   XV = t(X * V)                    [p x N]
#   XVX = t(X) %*% t(XV)            [p x p]
#   XVX_inv = solve(XVX)             [p x p]
#   XXVX_inv = X %*% XVX_inv         [N x p]
#   XVX_inv_XV = XXVX_inv * V        [N x p] -- element-wise multiply by V (broadcast)
#   S_a = colSums(X * res)           [p]

if (!is.null(obj.noK)) {
  cat("\nUsing obj.noK from model (pre-computed matrices)\n")
  XV <- obj.noK$XV
  XVX <- obj.noK$XVX
  XVX_inv <- obj.noK$XVX_inv
  XXVX_inv <- obj.noK$XXVX_inv
  XVX_inv_XV <- obj.noK$XVX_inv_XV
  S_a <- obj.noK$S_a

  # Some older models may not have XVX_inv stored; compute if needed
  if (is.null(XVX_inv)) {
    cat("  XVX_inv not in obj.noK, computing as solve(XVX)\n")
    XVX_inv <- solve(XVX)
  }
} else {
  cat("\nRecomputing obj.noK matrices (ScoreTest_NULL_Model equivalent)\n")
  XV <- t(X * V)                      # [p x N]
  XVX <- t(X) %*% t(XV)               # [p x p]
  XVX_inv <- solve(XVX)               # [p x p]
  XXVX_inv <- X %*% XVX_inv           # [N x p]
  XVX_inv_XV <- XXVX_inv * V          # [N x p] element-wise multiply each row by V
  S_a <- as.vector(colSums(X * res))   # [p]
}

# --- Step 5: Diagnostics ----------------------------------------------------

cat("\n--- Dimensions of extracted fields ---\n")
cat("  mu:          length =", length(mu), "\n")
cat("  res:         length =", length(res), "\n")
cat("  y:           length =", length(y), "\n")
cat("  V:           length =", length(V), "\n")
cat("  S_a:         length =", length(S_a), "\n")
cat("  X:           ", nrow(X), "x", ncol(X), "\n")
cat("  XV:          ", nrow(XV), "x", ncol(XV), "\n")
cat("  XVX:         ", nrow(XVX), "x", ncol(XVX), "\n")
cat("  XVX_inv:     ", nrow(XVX_inv), "x", ncol(XVX_inv), "\n")
cat("  XXVX_inv:    ", nrow(XXVX_inv), "x", ncol(XXVX_inv), "\n")
cat("  XVX_inv_XV:  ", nrow(XVX_inv_XV), "x", ncol(XVX_inv_XV), "\n")

cat("\n--- First 5 values of each vector ---\n")
cat("  mu[1:5]:  ", formatC(mu[1:min(5, N)], digits = 15, format = "g"), "\n")
cat("  res[1:5]: ", formatC(res[1:min(5, N)], digits = 15, format = "g"), "\n")
cat("  y[1:5]:   ", formatC(y[1:min(5, N)], digits = 15, format = "g"), "\n")
cat("  V[1:5]:   ", formatC(V[1:min(5, N)], digits = 15, format = "g"), "\n")
cat("  S_a:      ", formatC(S_a, digits = 15, format = "g"), "\n")

cat("\n--- XVX matrix ---\n")
print(XVX)

cat("\n--- XVX_inv matrix ---\n")
print(XVX_inv)

# Extract alpha (coefficients) if available
alpha <- modglmm$coefficients
if (is.null(alpha)) {
  # Try to get from LOCOResult
  if (is_loco && !is.null(loco_result$alpha0)) {
    alpha <- loco_result$alpha0
  }
}
if (is.null(alpha)) {
  cat("WARNING: coefficients (alpha) not found in model, using zeros\n")
  alpha <- rep(0.0, p)
}
cat("\nalpha: ", formatC(alpha, digits = 15, format = "g"), "\n")

# --- Step 6: Write .arma binary files ---------------------------------------

cat("\n--- Writing .arma files ---\n")

write_arma_vec(file.path(output_dir, "mu.arma"), mu)
cat("  Wrote mu.arma        [", length(mu), "x 1]\n")

write_arma_vec(file.path(output_dir, "res.arma"), res)
cat("  Wrote res.arma       [", length(res), "x 1]\n")

write_arma_vec(file.path(output_dir, "y.arma"), y)
cat("  Wrote y.arma         [", length(y), "x 1]\n")

write_arma_vec(file.path(output_dir, "V.arma"), V)
cat("  Wrote V.arma         [", length(V), "x 1]\n")

write_arma_vec(file.path(output_dir, "S_a.arma"), S_a)
cat("  Wrote S_a.arma       [", length(S_a), "x 1]\n")

write_arma_mat(file.path(output_dir, "X.arma"), X)
cat("  Wrote X.arma         [", nrow(X), "x", ncol(X), "]\n")

write_arma_mat(file.path(output_dir, "XVX.arma"), XVX)
cat("  Wrote XVX.arma       [", nrow(XVX), "x", ncol(XVX), "]\n")

write_arma_mat(file.path(output_dir, "XVX_inv.arma"), XVX_inv)
cat("  Wrote XVX_inv.arma   [", nrow(XVX_inv), "x", ncol(XVX_inv), "]\n")

write_arma_mat(file.path(output_dir, "XXVX_inv.arma"), XXVX_inv)
cat("  Wrote XXVX_inv.arma  [", nrow(XXVX_inv), "x", ncol(XXVX_inv), "]\n")

write_arma_mat(file.path(output_dir, "XV.arma"), XV)
cat("  Wrote XV.arma        [", nrow(XV), "x", ncol(XV), "]\n")

write_arma_mat(file.path(output_dir, "XVX_inv_XV.arma"), XVX_inv_XV)
cat("  Wrote XVX_inv_XV.arma [", nrow(XVX_inv_XV), "x", ncol(XVX_inv_XV), "]\n")

# --- Step 7: Write nullmodel.json -------------------------------------------

cat("\n--- Writing nullmodel.json ---\n")

# Format sampleIDs as JSON array
if (!is.null(sampleIDs)) {
  sample_ids_json <- paste0('"', sampleIDs, '"', collapse = ", ")
  sample_ids_line <- paste0('  "sampleIDs": [', sample_ids_json, ']')
} else {
  sample_ids_line <- '  "sampleIDs": []'
  cat("WARNING: sampleIDs not found in model\n")
}

# Format tau as JSON array with full precision
tau_json <- paste(formatC(tau, digits = 15, format = "g"), collapse = ", ")

# Format alpha as JSON array with full precision
alpha_json <- paste(formatC(alpha, digits = 15, format = "g"), collapse = ", ")

json_lines <- c(
  '{',
  paste0('  "traitType": "', traitType, '",'),
  paste0('  "n": ', N, ','),
  paste0('  "p": ', p, ','),
  paste0('  "tau": [', tau_json, '],'),
  paste0('  "alpha": [', alpha_json, '],'),
  '  "SPA_Cutoff": 2.0,',
  '  "impute_method": "mean",',
  '  "flagSparseGRM": false,',
  '  "isFastTest": true,',
  '  "isnoadjCov": false,',
  '  "pval_cutoff_for_fastTest": 0.05,',
  '  "isCondition": false,',
  '  "is_Firth_beta": true,',
  '  "pCutoffforFirth": 0.01,',
  paste0(sample_ids_line),
  '}'
)

json_path <- file.path(output_dir, "nullmodel.json")
writeLines(json_lines, json_path)
cat("  Wrote nullmodel.json\n")

# --- Step 8: Verification ---------------------------------------------------

cat("\n--- Verification: reading back mu.arma ---\n")

mu_readback <- load_arma_bin(file.path(output_dir, "mu.arma"))

if (length(mu_readback) != length(mu)) {
  cat("  ERROR: length mismatch! wrote", length(mu), "read back", length(mu_readback), "\n")
} else {
  max_diff <- max(abs(mu_readback - mu))
  cat("  Length: ", length(mu_readback), " (expected ", length(mu), ")\n")
  cat("  Max absolute difference: ", formatC(max_diff, digits = 3, format = "e"), "\n")
  cat("  mu_readback[1:5]: ", formatC(mu_readback[1:min(5, length(mu_readback))], digits = 15, format = "g"), "\n")
  if (max_diff < 1e-15) {
    cat("  PASS: mu.arma roundtrip is exact\n")
  } else if (max_diff < 1e-10) {
    cat("  PASS: mu.arma roundtrip within floating-point tolerance\n")
  } else {
    cat("  FAIL: mu.arma roundtrip difference too large\n")
  }
}

cat("\n--- Verification: reading back XVX.arma ---\n")

XVX_readback <- load_arma_bin(file.path(output_dir, "XVX.arma"))

if (!is.matrix(XVX_readback)) {
  cat("  WARNING: XVX_readback is not a matrix (dim = ", dim(XVX_readback), ")\n")
} else {
  max_diff <- max(abs(XVX_readback - XVX))
  cat("  Dims:  ", nrow(XVX_readback), "x", ncol(XVX_readback),
      " (expected ", nrow(XVX), "x", ncol(XVX), ")\n")
  cat("  Max absolute difference: ", formatC(max_diff, digits = 3, format = "e"), "\n")
  if (max_diff < 1e-10) {
    cat("  PASS: XVX.arma roundtrip OK\n")
  } else {
    cat("  FAIL: XVX.arma roundtrip difference too large\n")
  }
}

cat("\n--- Verification: reading back XVX_inv_XV.arma ---\n")

XVX_inv_XV_readback <- load_arma_bin(file.path(output_dir, "XVX_inv_XV.arma"))

if (is.matrix(XVX_inv_XV_readback)) {
  cat("  Dims:  ", nrow(XVX_inv_XV_readback), "x", ncol(XVX_inv_XV_readback),
      " (expected ", nrow(XVX_inv_XV), "x", ncol(XVX_inv_XV), ")\n")
  max_diff <- max(abs(XVX_inv_XV_readback - XVX_inv_XV))
  cat("  Max absolute difference: ", formatC(max_diff, digits = 3, format = "e"), "\n")
  if (max_diff < 1e-10) {
    cat("  PASS: XVX_inv_XV.arma roundtrip OK\n")
  } else {
    cat("  FAIL: XVX_inv_XV.arma roundtrip difference too large\n")
  }
} else {
  cat("  XVX_inv_XV read back as vector of length", length(XVX_inv_XV_readback), "\n")
}

# --- Summary -----------------------------------------------------------------

cat("\n=== Summary ===\n")
cat("Output directory: ", output_dir, "\n")
output_files <- list.files(output_dir, full.names = FALSE)
for (f in output_files) {
  fpath <- file.path(output_dir, f)
  fsize <- file.info(fpath)$size
  cat(sprintf("  %-20s  %8d bytes\n", f, fsize))
}

cat("\nnullmodel.json contents:\n")
cat(paste(readLines(json_path), collapse = "\n"), "\n")

cat("\nDone.\n")
