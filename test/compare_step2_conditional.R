#!/usr/bin/env Rscript
# ==============================================================================
# Compare R vs C++ SAIGE Step 2 Conditional Analysis Results
# ==============================================================================
#
# Compares:
#   1. Unconditional columns: BETA, SE, Tstat, var, p.value
#   2. Conditional columns: BETA_c, SE_c, Tstat_c, var_c, p.value_c
#
# Usage:
#   Rscript compare_step2_conditional.R
#
# Output:
#   test/RESULTS_step2_conditional.txt
#
# ==============================================================================

options(digits = 15)

# --- Paths ---
base_dir <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Step_2_Feb_11"

r_file   <- file.path(base_dir, "test/R/r_step2_conditional_results.txt")
cpp_file <- file.path(base_dir, "test/output/cpp_compare_conditional_results.txt")
out_file <- file.path(base_dir, "test/RESULTS_step2_conditional.txt")

cat("=== Step 2 Conditional Analysis: R vs C++ Comparison ===\n\n")
cat("R file:   ", r_file, "\n")
cat("C++ file: ", cpp_file, "\n")
cat("Output:   ", out_file, "\n\n")

# --- Load data ---
r_data   <- read.table(r_file,   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cpp_data <- read.table(cpp_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("R rows:   ", nrow(r_data), "\n")
cat("C++ rows: ", nrow(cpp_data), "\n")
cat("R cols:   ", paste(names(r_data), collapse=", "), "\n")
cat("C++ cols: ", paste(names(cpp_data), collapse=", "), "\n\n")

# --- Check row counts ---
if (nrow(r_data) != nrow(cpp_data)) {
  cat("WARNING: Row count mismatch!\n")
}

# --- Verify marker IDs match ---
if (!all(r_data$MarkerID == cpp_data$MarkerID)) {
  cat("WARNING: MarkerID mismatch!\n")
  mismatches <- which(r_data$MarkerID != cpp_data$MarkerID)
  cat("  First 5 mismatches at rows:", head(mismatches, 5), "\n")
} else {
  cat("All MarkerIDs match.\n\n")
}

# --- Identify conditioning markers ---
# These markers are conditioned on themselves, so their conditional
# values are essentially 0/infinity. We track them separately.
cond_markers <- c("rs1e+05", "rs8", "rs9")
is_cond <- r_data$MarkerID %in% cond_markers
n_cond <- sum(is_cond)
cat("Conditioning markers in output: ", n_cond, " (",
    paste(cond_markers[cond_markers %in% r_data$MarkerID], collapse=", "), ")\n\n")

# --- Comparison function ---
compare_column <- function(r_vec, cpp_vec, col_name, exclude_cond = FALSE) {
  if (exclude_cond) {
    r_vec <- r_vec[!is_cond]
    cpp_vec <- cpp_vec[!is_cond]
    label <- paste0(col_name, " (excl. cond. markers)")
  } else {
    label <- col_name
  }

  # Convert to numeric (handles "NA" strings and scientific notation)
  r_num <- suppressWarnings(as.numeric(r_vec))
  cpp_num <- suppressWarnings(as.numeric(cpp_vec))

  # Both NA
  both_na <- is.na(r_num) & is.na(cpp_num)
  # One NA
  one_na <- xor(is.na(r_num), is.na(cpp_num))
  # Both finite
  both_finite <- !is.na(r_num) & !is.na(cpp_num) & is.finite(r_num) & is.finite(cpp_num)
  # Both zero
  both_zero <- both_finite & (r_num == 0) & (cpp_num == 0)

  n_total <- length(r_num)
  n_both_na <- sum(both_na)
  n_one_na <- sum(one_na)
  n_both_finite <- sum(both_finite)
  n_both_zero <- sum(both_zero)

  # Compute differences for finite, non-zero values
  compare_idx <- both_finite & !both_zero

  if (sum(compare_idx) > 0) {
    abs_diff <- abs(r_num[compare_idx] - cpp_num[compare_idx])
    # Relative difference (avoid division by very small numbers)
    denom <- pmax(abs(r_num[compare_idx]), abs(cpp_num[compare_idx]), 1e-300)
    rel_diff <- abs_diff / denom

    exact_match <- sum(abs_diff == 0)
    max_abs_diff <- max(abs_diff)
    max_rel_diff <- max(rel_diff)
    mean_rel_diff <- mean(rel_diff)

    # P-value comparison: how many match to 6 significant figures?
    if (grepl("p.value", col_name)) {
      r_str <- sprintf("%.6E", r_num[compare_idx])
      cpp_str <- sprintf("%.6E", cpp_num[compare_idx])
      n_pval_exact <- sum(r_str == cpp_str)
    } else {
      n_pval_exact <- NA
    }

    result <- list(
      label = label,
      n_total = n_total,
      n_both_na = n_both_na,
      n_one_na = n_one_na,
      n_both_finite = n_both_finite,
      n_compared = sum(compare_idx),
      n_exact = exact_match,
      max_abs = max_abs_diff,
      max_rel = max_rel_diff,
      mean_rel = mean_rel_diff,
      n_pval_exact = n_pval_exact,
      PASS = (max_rel_diff < 1e-4) & (n_one_na == 0)
    )
  } else {
    result <- list(
      label = label,
      n_total = n_total,
      n_both_na = n_both_na,
      n_one_na = n_one_na,
      n_both_finite = n_both_finite,
      n_compared = 0,
      n_exact = 0,
      max_abs = 0,
      max_rel = 0,
      mean_rel = 0,
      n_pval_exact = NA,
      PASS = (n_one_na == 0)
    )
  }

  return(result)
}

# --- Compare unconditional columns ---
cat("=========================================\n")
cat("UNCONDITIONAL COLUMNS (should match Test 1)\n")
cat("=========================================\n\n")

uncond_cols <- c("BETA", "SE", "Tstat", "var", "p.value")
uncond_results <- list()

for (col in uncond_cols) {
  res <- compare_column(r_data[[col]], cpp_data[[col]], col)
  uncond_results[[col]] <- res

  status <- if (res$PASS) "PASS" else "FAIL"
  cat(sprintf("%-12s: %s  (n=%d, exact=%d, max_rel=%.2e)\n",
              col, status, res$n_compared, res$n_exact, res$max_rel))
}

# --- Compare conditional columns ---
cat("\n=========================================\n")
cat("CONDITIONAL COLUMNS (all markers)\n")
cat("=========================================\n\n")

cond_cols <- c("BETA_c", "SE_c", "Tstat_c", "var_c", "p.value_c")
cond_results_all <- list()

for (col in cond_cols) {
  res <- compare_column(r_data[[col]], cpp_data[[col]], col)
  cond_results_all[[col]] <- res

  status <- if (res$PASS) "PASS" else "FAIL"
  cat(sprintf("%-12s: %s  (n=%d, exact=%d, max_rel=%.2e, max_abs=%.2e)\n",
              col, status, res$n_compared, res$n_exact, res$max_rel, res$max_abs))
}

# --- Compare conditional columns EXCLUDING conditioning markers ---
cat("\n=========================================\n")
cat("CONDITIONAL COLUMNS (excl. conditioning markers)\n")
cat("=========================================\n\n")

cond_results_excl <- list()

for (col in cond_cols) {
  res <- compare_column(r_data[[col]], cpp_data[[col]], col, exclude_cond = TRUE)
  cond_results_excl[[col]] <- res

  status <- if (res$PASS) "PASS" else "FAIL"
  cat(sprintf("%-12s: %s  (n=%d, exact=%d, max_rel=%.2e, max_abs=%.2e)\n",
              col, status, res$n_compared, res$n_exact, res$max_rel, res$max_abs))
}

# --- Show worst mismatches for conditional columns ---
cat("\n=========================================\n")
cat("WORST MISMATCHES (conditional, excl. cond markers)\n")
cat("=========================================\n\n")

for (col in cond_cols) {
  r_num <- suppressWarnings(as.numeric(r_data[[col]][!is_cond]))
  cpp_num <- suppressWarnings(as.numeric(cpp_data[[col]][!is_cond]))
  ids <- r_data$MarkerID[!is_cond]

  both_finite <- !is.na(r_num) & !is.na(cpp_num) & is.finite(r_num) & is.finite(cpp_num)
  both_nonzero <- both_finite & (r_num != 0) & (cpp_num != 0)

  if (sum(both_nonzero) > 0) {
    denom <- pmax(abs(r_num[both_nonzero]), abs(cpp_num[both_nonzero]))
    rel_diff <- abs(r_num[both_nonzero] - cpp_num[both_nonzero]) / denom

    worst_idx <- which(both_nonzero)[order(rel_diff, decreasing = TRUE)][1:min(5, sum(both_nonzero))]

    cat(sprintf("%-12s: Top %d worst relative differences:\n", col, length(worst_idx)))
    for (i in worst_idx) {
      rd <- abs(r_num[i] - cpp_num[i]) / max(abs(r_num[i]), abs(cpp_num[i]))
      cat(sprintf("  %s: R=%.10g  C++=%.10g  rel_diff=%.4e\n",
                  ids[i], r_num[i], cpp_num[i], rd))
    }
    cat("\n")
  }
}

# --- Show conditioning markers specifically ---
cat("=========================================\n")
cat("CONDITIONING MARKERS (conditioned on themselves)\n")
cat("=========================================\n\n")

for (m in cond_markers) {
  r_idx <- which(r_data$MarkerID == m)
  cpp_idx <- which(cpp_data$MarkerID == m)

  if (length(r_idx) > 0 && length(cpp_idx) > 0) {
    cat(sprintf("Marker: %s (index R=%d, C++=%d)\n", m, r_idx, cpp_idx))
    for (col in cond_cols) {
      cat(sprintf("  %-12s: R=%-20s  C++=%-20s\n",
                  col, as.character(r_data[[col]][r_idx]), as.character(cpp_data[[col]][cpp_idx])))
    }
    cat("\n")
  }
}

# --- Overall verdict ---
cat("=========================================\n")
cat("OVERALL VERDICT\n")
cat("=========================================\n\n")

# Unconditional: strict match (same as Test 1)
uncond_pass <- all(sapply(uncond_results, function(x) x$PASS))
cat("Unconditional columns (BETA, SE, Tstat, var, p.value):\n")
cat("  Status: ", if (uncond_pass) "PASS" else "FAIL", "\n\n")

# Conditional: strict match for non-conditioning markers
cond_excl_pass <- all(sapply(cond_results_excl, function(x) x$PASS))
cat("Conditional columns (excl. conditioning markers):\n")
cat("  Status: ", if (cond_excl_pass) "PASS" else "FAIL", "\n\n")

overall_pass <- uncond_pass && cond_excl_pass
cat("OVERALL: ", if (overall_pass) "PASS" else "FAIL", "\n")

# --- Write results file ---
sink(out_file)
cat("SAIGE Step 2 Conditional Analysis: R vs C++ Comparison Results\n")
cat("=============================================================\n")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("R file:   %s\n", r_file))
cat(sprintf("C++ file: %s\n", cpp_file))
cat(sprintf("R rows:   %d\n", nrow(r_data)))
cat(sprintf("C++ rows: %d\n", nrow(cpp_data)))
cat(sprintf("Conditioning markers: %s\n", paste(cond_markers, collapse=", ")))
cat("\n")

cat("UNCONDITIONAL COLUMNS:\n")
for (col in uncond_cols) {
  res <- uncond_results[[col]]
  cat(sprintf("  %-12s: %s  (n=%d, exact=%d, max_rel=%.2e)\n",
              col, if (res$PASS) "PASS" else "FAIL",
              res$n_compared, res$n_exact, res$max_rel))
}

cat("\nCONDITIONAL COLUMNS (excl. conditioning markers):\n")
for (col in cond_cols) {
  res <- cond_results_excl[[col]]
  cat(sprintf("  %-12s: %s  (n=%d, exact=%d, max_rel=%.2e, max_abs=%.2e)\n",
              col, if (res$PASS) "PASS" else "FAIL",
              res$n_compared, res$n_exact, res$max_rel, res$max_abs))
}

cat("\nCONDITIONAL COLUMNS (all markers, incl. conditioning markers):\n")
for (col in cond_cols) {
  res <- cond_results_all[[col]]
  cat(sprintf("  %-12s: %s  (n=%d, exact=%d, max_rel=%.2e, max_abs=%.2e)\n",
              col, if (res$PASS) "PASS" else "FAIL",
              res$n_compared, res$n_exact, res$max_rel, res$max_abs))
}

cat(sprintf("\nOVERALL VERDICT: %s\n", if (overall_pass) "PASS" else "FAIL"))
sink()

cat("\nResults written to:", out_file, "\n")
cat("\n=== Done ===\n")
