#!/usr/bin/env Rscript
# ==============================================================================
# compare_nullmodels.R
# ==============================================================================
# Compare C++ Step 1 null model output vs R Step 1 null model (converted from .rda)
#
# C++ Step 1 outputs files as: {prefix}.{field}.arma (flat naming in one directory)
# R converted model: {dir}/{field}.arma (directory-based naming)
#
# Also compares JSON scalars (tau/theta, trait, n, p, etc.)
# ==============================================================================

options(digits = 15)

cat("====================================================================\n")
cat("Null Model Comparison: C++ Step 1 vs R Step 1 (converted from .rda)\n")
cat("====================================================================\n\n")

# --- Helper: read Armadillo binary file (.arma) ---
load_arma_bin <- function(filepath) {
  if (!file.exists(filepath)) {
    cat("  WARNING: File not found:", filepath, "\n")
    return(NULL)
  }
  con <- file(filepath, "rb")
  on.exit(close(con))
  header <- readLines(con, n = 1)
  dims_line <- readLines(con, n = 1)
  dims <- as.integer(strsplit(trimws(dims_line), "\\s+")[[1]])
  nr <- dims[1]; nc <- dims[2]
  n_elements <- nr * nc
  values <- readBin(con, what = "double", n = n_elements, size = 8, endian = "little")
  if (nc == 1) return(as.numeric(values))
  matrix(values, nrow = nr, ncol = nc, byrow = FALSE)
}

# --- Paths ---
# R converted model (from .rda with LOCO=TRUE, chr 1)
r_dir <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Step_2_Feb_11/test/data/nullmodel_from_rda"

# C++ Step 1 output (multiple configs available)
cpp_out_dir <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/output"

# We'll compare multiple C++ configs
cpp_configs <- list(
  sparse_x1x2 = "cpp_sparse_x1x2_output",
  dense_x1x2  = "cpp_dense_x1x2_output",
  sparse_x1   = "cpp_sparse_x1_output",
  dense_x1    = "cpp_dense_x1_output"
)

# Fields to compare (vectors and matrices)
arma_fields <- c("mu", "res", "y", "V", "S_a", "X", "XV", "XVX", "XVX_inv", "XXVX_inv", "XVX_inv_XV")

# --- Helper: compare two arrays element-wise ---
compare_field <- function(r_data, cpp_data, field_name) {
  if (is.null(r_data) || is.null(cpp_data)) {
    cat(sprintf("  %-15s SKIP (missing data)\n", field_name))
    return(NULL)
  }

  # Check dimensions
  r_dim <- if (is.matrix(r_data)) dim(r_data) else c(length(r_data), 1)
  cpp_dim <- if (is.matrix(cpp_data)) dim(cpp_data) else c(length(cpp_data), 1)

  if (!identical(r_dim, cpp_dim)) {
    cat(sprintf("  %-15s DIM MISMATCH: R=%dx%d  C++=%dx%d\n",
                field_name, r_dim[1], r_dim[2], cpp_dim[1], cpp_dim[2]))
    return(NULL)
  }

  r_vec <- as.vector(r_data)
  cpp_vec <- as.vector(cpp_data)

  abs_diff <- abs(r_vec - cpp_vec)
  max_abs_diff <- max(abs_diff)
  mean_abs_diff <- mean(abs_diff)

  # Relative difference (avoiding division by zero)
  denom <- pmax(abs(r_vec), abs(cpp_vec), 1e-300)
  rel_diff <- abs_diff / denom
  max_rel_diff <- max(rel_diff)
  mean_rel_diff <- mean(rel_diff)

  # Correlation
  corr <- if (length(r_vec) > 1) cor(r_vec, cpp_vec) else NA

  verdict <- if (max_rel_diff < 1e-10) {
    "EXACT"
  } else if (max_rel_diff < 1e-6) {
    "VERY_CLOSE"
  } else if (max_rel_diff < 1e-3) {
    "CLOSE"
  } else if (max_rel_diff < 0.01) {
    "SIMILAR"
  } else {
    "DIFFERENT"
  }

  cat(sprintf("  %-15s [%dx%d] max_abs=%.3e  max_rel=%.3e  mean_rel=%.3e  corr=%.15f  %s\n",
              field_name, r_dim[1], r_dim[2], max_abs_diff, max_rel_diff, mean_rel_diff,
              ifelse(is.na(corr), 0, corr), verdict))

  return(list(field = field_name, max_abs = max_abs_diff, max_rel = max_rel_diff,
              mean_rel = mean_rel_diff, corr = corr, verdict = verdict))
}

# --- Load R model data ---
cat("Loading R converted model from:\n  ", r_dir, "\n\n")

r_data <- list()
for (field in arma_fields) {
  filepath <- file.path(r_dir, paste0(field, ".arma"))
  r_data[[field]] <- load_arma_bin(filepath)
}

# Load R JSON
r_json <- jsonlite::fromJSON(file.path(r_dir, "nullmodel.json"))
cat("R model JSON:\n")
cat("  traitType: ", r_json$traitType, "\n")
cat("  n:         ", r_json$n, "\n")
cat("  p:         ", r_json$p, "\n")
cat("  tau:       ", paste(formatC(r_json$tau, digits = 15, format = "g"), collapse = ", "), "\n")
cat("  alpha:     ", paste(formatC(r_json$alpha, digits = 15, format = "g"), collapse = ", "), "\n")
cat("  LOCO:      TRUE (extracted chr 1 from .rda)\n\n")

# --- Compare with each C++ config ---
all_results <- list()

for (cfg_name in names(cpp_configs)) {
  prefix <- cpp_configs[[cfg_name]]

  cat("====================================================================\n")
  cat("Comparing R model vs C++ Step 1:", cfg_name, "\n")
  cat("  C++ prefix:", prefix, "\n")
  cat("====================================================================\n")

  # Load C++ JSON
  cpp_json_path <- file.path(cpp_out_dir, paste0(prefix, ".nullmodel.json"))
  if (!file.exists(cpp_json_path)) {
    cat("  WARNING: JSON not found:", cpp_json_path, "\n\n")
    next
  }
  cpp_json <- jsonlite::fromJSON(cpp_json_path)

  cat("\nC++ JSON:\n")
  cat("  trait:  ", cpp_json$trait, "\n")
  cat("  n:      ", cpp_json$n, "\n")
  cat("  p:      ", cpp_json$p, "\n")
  cat("  theta:  ", paste(formatC(cpp_json$theta, digits = 15, format = "g"), collapse = ", "), "\n")
  cat("  alpha:  ", paste(formatC(cpp_json$alpha, digits = 15, format = "g"), collapse = ", "), "\n")
  cat("  loco:   ", cpp_json$loco, "\n\n")

  # Check dimension compatibility
  if (r_json$p != cpp_json$p) {
    cat("  *** DIMENSION MISMATCH: R p=", r_json$p, " vs C++ p=", cpp_json$p, " ***\n")
    cat("  (Skipping .arma field comparison for incompatible dimensions)\n\n")

    # Still compare scalars
    cat("  Scalar comparison:\n")
    cat("  tau/theta:\n")
    cat("    R  tau  = [", paste(formatC(r_json$tau, digits = 15, format = "g"), collapse = ", "), "]\n")
    cat("    C++ theta = [", paste(formatC(cpp_json$theta, digits = 15, format = "g"), collapse = ", "), "]\n")

    next
  }

  # Compare scalars
  cat("  Scalar comparison (tau vs theta):\n")
  for (i in 1:2) {
    r_val <- r_json$tau[i]
    cpp_val <- cpp_json$theta[i]
    abs_d <- abs(r_val - cpp_val)
    rel_d <- abs_d / max(abs(r_val), abs(cpp_val), 1e-300)
    cat(sprintf("    tau/theta[%d]: R=%.15g  C++=%.15g  abs_diff=%.3e  rel_diff=%.3e\n",
                i, r_val, cpp_val, abs_d, rel_d))
  }
  cat("\n")

  # Compare alpha (coefficients)
  cat("  Alpha comparison:\n")
  cat("    R  alpha = [", paste(formatC(r_json$alpha, digits = 15, format = "g"), collapse = ", "), "]\n")
  cat("    C++ alpha = [", paste(formatC(cpp_json$alpha, digits = 15, format = "g"), collapse = ", "), "]\n\n")

  # Compare .arma fields
  cat("  Field-by-field .arma comparison:\n")
  cfg_results <- list()

  for (field in arma_fields) {
    cpp_filepath <- file.path(cpp_out_dir, paste0(prefix, ".", field, ".arma"))
    cpp_field <- load_arma_bin(cpp_filepath)

    result <- compare_field(r_data[[field]], cpp_field, field)
    if (!is.null(result)) {
      cfg_results[[field]] <- result
    }
  }

  all_results[[cfg_name]] <- cfg_results
  cat("\n")
}

# --- Summary table ---
cat("====================================================================\n")
cat("SUMMARY TABLE: Max relative differences by config and field\n")
cat("====================================================================\n\n")

# Only configs with p=3 (compatible dimensions)
for (cfg_name in names(all_results)) {
  cat("Config:", cfg_name, "\n")
  cat(sprintf("  %-15s  %12s  %12s  %s\n", "Field", "Max_RelDiff", "Correlation", "Verdict"))
  cat(sprintf("  %-15s  %12s  %12s  %s\n", "-----", "-----------", "-----------", "-------"))
  for (field in names(all_results[[cfg_name]])) {
    r <- all_results[[cfg_name]][[field]]
    cat(sprintf("  %-15s  %12.3e  %12.10f  %s\n",
                r$field, r$max_rel, ifelse(is.na(r$corr), 0, r$corr), r$verdict))
  }
  cat("\n")
}

cat("====================================================================\n")
cat("KEY OBSERVATIONS:\n")
cat("====================================================================\n")
cat("\n")
cat("CRITICAL CONFIGURATION DIFFERENCES:\n")
cat("  - R model: LOCO=TRUE (per-chromosome model, extracted chr 1)\n")
cat("  - C++ Step 1: loco=false (global model, all chromosomes)\n")
cat("  - This means the fitted values (mu), residuals (res), and derived\n")
cat("    matrices will differ because LOCO removes the effect of the\n")
cat("    test chromosome's variants from the null model.\n")
cat("  - y (phenotype) and X (design matrix) should be identical.\n")
cat("  - tau/theta values will differ because LOCO refits without chr 1.\n")
cat("\n")
cat("Done.\n")
