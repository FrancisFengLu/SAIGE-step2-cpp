#!/usr/bin/env Rscript
# ==============================================================================
# convert_rda_to_arma_sparse.R
# ==============================================================================
#
# Purpose:
#   Convert a SAIGE .rda null model (fitted with sparseGRM) to Armadillo binary
#   files (.arma) and a nullmodel.json, so that the standalone C++ Step 2 binary
#   can load the exact same null model that R SAIGE Step 2 uses.
#
#   This script extends convert_rda_to_arma.R by additionally:
#   1. Reading the sparse GRM .mtx file
#   2. Applying the setSparseSigma_new transformation (subset, scale, diagonal)
#   3. Saving sparseGRM_locationMat.arma and sparseGRM_valueVec.arma
#   4. Setting flagSparseGRM=true in nullmodel.json
#
# Usage:
#   Rscript convert_rda_to_arma_sparse.R
#
#   Or via pixi:
#   ~/.pixi/bin/pixi run --manifest-path=.../pixi.toml Rscript convert_rda_to_arma_sparse.R
#
# Input:
#   SAIGE extdata example_quantitative_sparseGRM.rda
#   Sparse GRM .mtx file + sample ID file
#
# Output:
#   Directory: test/data/nullmodel_from_rda_sparse/
#     mu.arma, res.arma, y.arma, V.arma, S_a.arma,
#     X.arma, XVX.arma, XVX_inv.arma, XXVX_inv.arma, XV.arma, XVX_inv_XV.arma,
#     sparseGRM_locationMat.arma, sparseGRM_valueVec.arma,
#     nullmodel.json
#
# ==============================================================================

options(digits = 15)

cat("=== convert_rda_to_arma_sparse.R ===\n\n")

# --- Paths -------------------------------------------------------------------

rda_path <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/SAIGE/extdata/output/example_quantitative_sparseGRM.rda"
sparseGRM_file <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/SAIGE/extdata/output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
sparseGRM_sampleID_file <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/SAIGE/extdata/output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
varianceRatio_file <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/SAIGE/extdata/output/example_quantitative_sparseGRM.varianceRatio.txt"
output_dir <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Step_2_Feb_11/test/data/nullmodel_from_rda_sparse"

# Relatedness cutoff (same default as SAIGE Step 2)
relatednessCutoff <- 0.125

cat("Input .rda:         ", rda_path, "\n")
cat("Sparse GRM .mtx:    ", sparseGRM_file, "\n")
cat("Sparse GRM IDs:     ", sparseGRM_sampleID_file, "\n")
cat("Variance ratio:     ", varianceRatio_file, "\n")
cat("Output dir:         ", output_dir, "\n")
cat("relatednessCutoff:  ", relatednessCutoff, "\n\n")

for (f in c(rda_path, sparseGRM_file, sparseGRM_sampleID_file, varianceRatio_file)) {
  if (!file.exists(f)) {
    stop("ERROR: File not found: ", f)
  }
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Helper: write Armadillo binary format -----------------------------------

# Write a vector as ARMA_MAT_BIN_FN008 format (column vector)
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
  writeBin(as.double(mat), con, size = 8, endian = "little")
}

# Write a uword matrix as ARMA_MAT_BIN_IU008 format (column-major, 8-byte unsigned integers)
# This matches arma::umat::save(arma_binary)
write_arma_umat <- function(filepath, mat) {
  mat <- as.matrix(mat)
  nr <- nrow(mat)
  nc <- ncol(mat)
  con <- file(filepath, "wb")
  on.exit(close(con))
  # ARMA_MAT_BIN_IU008 is armadillo's header for unsigned 64-bit integer matrices
  writeLines(c("ARMA_MAT_BIN_IU008", paste(nr, nc)), con)
  # Write as 8-byte unsigned integers, column-major
  # R doesn't have unsigned 64-bit, but we can write as raw bytes
  for (c_idx in 1:nc) {
    for (r_idx in 1:nr) {
      val <- as.integer(mat[r_idx, c_idx])
      # Pack as 8-byte little-endian unsigned integer
      raw_bytes <- writeBin(as.double(val), raw(), size = 8, endian = "little")
      # Actually, we need proper unsigned int64 encoding
      # Use raw byte packing: val as 8-byte little-endian
      bytes <- raw(8)
      v <- as.integer(mat[r_idx, c_idx])
      for (b in 1:8) {
        bytes[b] <- as.raw(bitwAnd(v, 0xFF))
        v <- bitwShiftR(v, 8)
      }
      writeBin(bytes, con)
    }
  }
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
  if (grepl("IU008", header)) {
    # Unsigned integer matrix: read as raw bytes and convert
    values <- integer(n_elements)
    for (i in 1:n_elements) {
      raw_bytes <- readBin(con, what = "raw", n = 8)
      # Convert 8 little-endian bytes to integer (only works for values < 2^31)
      val <- 0
      for (b in 8:1) {
        val <- val * 256 + as.integer(raw_bytes[b])
      }
      values[i] <- val
    }
    if (ncol == 1) return(as.integer(values))
    return(matrix(values, nrow = nrow, ncol = ncol, byrow = FALSE))
  }
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

# Find the main model object
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

  loco_result <- modglmm$LOCOResult[[target_chr]]
  mu <- as.vector(loco_result$fitted.values)
  res <- as.vector(loco_result$residuals)

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

if (!is.null(obj.noK)) {
  cat("\nUsing obj.noK from model (pre-computed matrices)\n")
  XV <- obj.noK$XV
  XVX <- obj.noK$XVX
  XVX_inv <- obj.noK$XVX_inv
  XXVX_inv <- obj.noK$XXVX_inv
  XVX_inv_XV <- obj.noK$XVX_inv_XV
  S_a <- obj.noK$S_a

  if (is.null(XVX_inv)) {
    cat("  XVX_inv not in obj.noK, computing as solve(XVX)\n")
    XVX_inv <- solve(XVX)
  }
} else {
  cat("\nRecomputing obj.noK matrices (ScoreTest_NULL_Model equivalent)\n")
  XV <- t(X * V)
  XVX <- t(X) %*% t(XV)
  XVX_inv <- solve(XVX)
  XXVX_inv <- X %*% XVX_inv
  XVX_inv_XV <- XXVX_inv * V
  S_a <- as.vector(colSums(X * res))
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
  if (is_loco && !is.null(loco_result$alpha0)) {
    alpha <- loco_result$alpha0
  }
}
if (is.null(alpha)) {
  cat("WARNING: coefficients (alpha) not found in model, using zeros\n")
  alpha <- rep(0.0, p)
} else {
  alpha <- as.vector(alpha)
}
cat("\nalpha: ", formatC(alpha, digits = 15, format = "g"), "\n")

# --- Step 6: Process sparse GRM (setSparseSigma_new equivalent) -------------

cat("\n--- Processing sparse GRM ---\n")

# Load required library for sparse matrix reading
library(Matrix)

# Read sparse GRM in MatrixMarket format
cat("Reading sparse GRM from: ", sparseGRM_file, "\n")
sparseGRM <- readMM(sparseGRM_file)
cat("  sparseGRM dimensions: ", dim(sparseGRM), "\n")
cat("  sparseGRM non-zeros:  ", length(sparseGRM@x), "\n")

# Read sample IDs for the GRM
sparseGRMSampleID <- read.table(sparseGRM_sampleID_file, header = FALSE,
                                 stringsAsFactors = FALSE, colClasses = "character")
colnames(sparseGRMSampleID) <- "sampleID"
sparseGRMSampleID$IndexGRM <- 1:nrow(sparseGRMSampleID)
cat("  GRM sample IDs: ", nrow(sparseGRMSampleID), "\n")

# Build mapping: model sample IDs -> GRM indices
sampleInModel <- data.frame(IID = sampleIDs, stringsAsFactors = FALSE)
sampleInModel$IndexInModel <- seq_len(nrow(sampleInModel))
cat("  Model sample IDs: ", nrow(sampleInModel), "\n")

mergeID <- merge(sampleInModel, sparseGRMSampleID, by.x = "IID", by.y = "sampleID")
mergeID <- mergeID[order(mergeID$IndexInModel), ]
cat("  Matched samples: ", nrow(mergeID), "\n")

if (nrow(mergeID) < nrow(sampleInModel)) {
  stop("ERROR: ", nrow(sampleInModel) - nrow(mergeID),
       " samples from null model not found in sparse GRM")
}

# Subset GRM to model samples (reorder to match model sample order)
indexIDofGRM <- mergeID$IndexGRM
sparseGRM <- sparseGRM[indexIDofGRM, indexIDofGRM]
cat("  Subsetted sparseGRM: ", dim(sparseGRM), " with ", length(sparseGRM@x), " non-zeros\n")

# Remove elements below relatedness cutoff (but keep diagonal)
# Note: relatednessCutoff applies to off-diagonal elements
removeIndex <- which(sparseGRM@x < relatednessCutoff)
if (length(removeIndex) > 0) {
  cat("  Removing ", length(removeIndex), " elements < relatednessCutoff (", relatednessCutoff, ")\n")
  sparseGRM@x <- sparseGRM@x[-removeIndex]
  sparseGRM@i <- sparseGRM@i[-removeIndex]
  sparseGRM@j <- sparseGRM@j[-removeIndex]
}
cat("  After filtering: ", length(sparseGRM@x), " non-zeros\n")

# Transform to sparseSigma:
# 1. Scale all elements by tau[2] (the genetic variance component)
sparseGRM@x <- sparseGRM@x * tau[2]
cat("  Scaled by tau[2] = ", tau[2], "\n")

# 2. Add to diagonal based on trait type
# For quantitative: diag += tau[1]
# For binary:       diag += 1/W where W = mu*(1-mu) = mu2
diag_indices <- which(sparseGRM@i == sparseGRM@j)
cat("  Number of diagonal elements: ", length(diag_indices), "\n")

if (traitType == "binary" || traitType == "survival") {
  # W = mu2 = mu*(1-mu), diagonal gets 1/W added
  # The sparseGRM@i indices are 0-based
  diag_sample_indices <- sparseGRM@i[diag_indices] + 1  # convert to 1-based
  sparseGRM@x[diag_indices] <- sparseGRM@x[diag_indices] + 1.0 / mu2[diag_sample_indices]
  cat("  Binary/survival: added 1/mu2 to diagonal\n")
} else if (traitType == "quantitative") {
  sparseGRM@x[diag_indices] <- sparseGRM@x[diag_indices] + tau[1]
  cat("  Quantitative: added tau[1] = ", tau[1], " to diagonal\n")
}

# Extract location matrix and value vector
# locations: 2 x nnz matrix with row 1 = i (0-indexed), row 2 = j (0-indexed)
# values: nnz-length vector
locations <- rbind(sparseGRM@i, sparseGRM@j)  # Already 0-indexed from Matrix package
values <- sparseGRM@x
nSubj <- dim(sparseGRM)[1]

cat("\n  sparseSigma summary:\n")
cat("    nSubj (dimNum):    ", nSubj, "\n")
cat("    locationMat shape: ", nrow(locations), " x ", ncol(locations), "\n")
cat("    valueVec length:   ", length(values), "\n")
cat("    value range:       [", min(values), ", ", max(values), "]\n")
cat("    first 5 locations: i=", locations[1, 1:min(5, ncol(locations))],
    " j=", locations[2, 1:min(5, ncol(locations))], "\n")
cat("    first 5 values:    ", values[1:min(5, length(values))], "\n")

# --- Step 7: Write .arma binary files ---------------------------------------

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

# Write offset if present
if (!is.null(modglmm$offset)) {
  model_offset <- as.vector(modglmm$offset)
  write_arma_vec(file.path(output_dir, "offset.arma"), model_offset)
  cat("  Wrote offset.arma    [", length(model_offset), "x 1]\n")
}

# Write sparse GRM files
write_arma_umat(file.path(output_dir, "sparseGRM_locationMat.arma"), locations)
cat("  Wrote sparseGRM_locationMat.arma [", nrow(locations), "x", ncol(locations), "]\n")

write_arma_vec(file.path(output_dir, "sparseGRM_valueVec.arma"), values)
cat("  Wrote sparseGRM_valueVec.arma    [", length(values), "x 1]\n")

# Copy the variance ratio file
file.copy(varianceRatio_file, file.path(output_dir, "varianceRatio.txt"), overwrite = TRUE)
cat("  Copied varianceRatio.txt\n")

# Also copy the raw sparse GRM files for reference
file.copy(sparseGRM_file, file.path(output_dir, "sparseGRM.mtx"), overwrite = TRUE)
file.copy(sparseGRM_sampleID_file, file.path(output_dir, "sparseGRM.mtx.sampleIDs.txt"), overwrite = TRUE)
cat("  Copied raw sparseGRM.mtx and sampleIDs.txt for reference\n")

# --- Step 8: Write nullmodel.json -------------------------------------------

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
  '  "flagSparseGRM": true,',
  '  "isFastTest": true,',
  '  "isnoadjCov": false,',
  '  "pval_cutoff_for_fastTest": 0.05,',
  '  "isCondition": false,',
  '  "is_Firth_beta": false,',
  '  "pCutoffforFirth": 0.01,',
  paste0(sample_ids_line),
  '}'
)

json_path <- file.path(output_dir, "nullmodel.json")
writeLines(json_lines, json_path)
cat("  Wrote nullmodel.json\n")

# --- Step 9: Verification ---------------------------------------------------

cat("\n--- Verification: reading back mu.arma ---\n")

mu_readback <- load_arma_bin(file.path(output_dir, "mu.arma"))

if (length(mu_readback) != length(mu)) {
  cat("  ERROR: length mismatch! wrote", length(mu), "read back", length(mu_readback), "\n")
} else {
  max_diff <- max(abs(mu_readback - mu))
  cat("  Length: ", length(mu_readback), " (expected ", length(mu), ")\n")
  cat("  Max absolute difference: ", formatC(max_diff, digits = 3, format = "e"), "\n")
  if (max_diff < 1e-15) {
    cat("  PASS: mu.arma roundtrip is exact\n")
  } else if (max_diff < 1e-10) {
    cat("  PASS: mu.arma roundtrip within floating-point tolerance\n")
  } else {
    cat("  FAIL: mu.arma roundtrip difference too large\n")
  }
}

cat("\n--- Verification: reading back sparseGRM_valueVec.arma ---\n")

values_readback <- load_arma_bin(file.path(output_dir, "sparseGRM_valueVec.arma"))

if (length(values_readback) != length(values)) {
  cat("  ERROR: length mismatch! wrote", length(values), "read back", length(values_readback), "\n")
} else {
  max_diff <- max(abs(values_readback - values))
  cat("  Length: ", length(values_readback), " (expected ", length(values), ")\n")
  cat("  Max absolute difference: ", formatC(max_diff, digits = 3, format = "e"), "\n")
  if (max_diff < 1e-15) {
    cat("  PASS: sparseGRM_valueVec.arma roundtrip is exact\n")
  } else if (max_diff < 1e-10) {
    cat("  PASS: sparseGRM_valueVec.arma roundtrip within floating-point tolerance\n")
  } else {
    cat("  FAIL: sparseGRM_valueVec.arma roundtrip difference too large\n")
  }
}

cat("\n--- Verification: reading back sparseGRM_locationMat.arma ---\n")

locations_readback <- load_arma_bin(file.path(output_dir, "sparseGRM_locationMat.arma"))

if (is.matrix(locations_readback)) {
  cat("  Dims:  ", nrow(locations_readback), "x", ncol(locations_readback),
      " (expected ", nrow(locations), "x", ncol(locations), ")\n")
  max_diff <- max(abs(locations_readback - locations))
  cat("  Max absolute difference: ", formatC(max_diff, digits = 3, format = "e"), "\n")
  if (max_diff == 0) {
    cat("  PASS: sparseGRM_locationMat.arma roundtrip is exact\n")
  } else {
    cat("  FAIL: sparseGRM_locationMat.arma roundtrip mismatch\n")
  }
} else {
  cat("  Read back as vector length", length(locations_readback), " (may be OK if only 1 column)\n")
}

# --- Summary -----------------------------------------------------------------

cat("\n=== Summary ===\n")
cat("Output directory: ", output_dir, "\n")
cat("Trait type: ", traitType, "\n")
cat("LOCO: ", is_loco, "\n")
cat("Sparse GRM: YES (flagSparseGRM=true)\n")
cat("  dimNum (nSubj): ", nSubj, "\n")
cat("  locationMat: ", nrow(locations), "x", ncol(locations), "\n")
cat("  valueVec: ", length(values), " elements\n")

output_files <- list.files(output_dir, full.names = FALSE)
for (f in sort(output_files)) {
  fpath <- file.path(output_dir, f)
  fsize <- file.info(fpath)$size
  cat(sprintf("  %-40s  %8d bytes\n", f, fsize))
}

cat("\nnullmodel.json contents:\n")
cat(paste(readLines(json_path), collapse = "\n"), "\n")

cat("\nDone.\n")
