#!/usr/bin/env Rscript
# ==============================================================================
# SAIGE Step 2: R vs C++ Result Comparison (Sparse GRM + noadjCov)
# ==============================================================================
# Compares R SPAGMMATtest output with C++ saige-step2 output
# for the quantitative sparse GRM example model with is_noadjCov = TRUE.
#
# Usage: Rscript compare_step2_sparse_noadjcov.R
# ==============================================================================

options(digits = 15)

# --- Paths ------------------------------------------------------------------
base_dir <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Step_2_Feb_11"

r_file   <- file.path(base_dir, "test/R/r_step2_sparse_noadjcov_results.txt")
cpp_file <- file.path(base_dir, "test/output/cpp_compare_sparse_noadjcov_results.txt")
out_file <- file.path(base_dir, "test/RESULTS_step2_sparse_noadjcov.txt")

cat("=== SAIGE Step 2: R vs C++ Comparison (Sparse GRM + noadjCov) ===\n\n")
cat("R results:   ", r_file, "\n")
cat("C++ results: ", cpp_file, "\n")
cat("Output:      ", out_file, "\n\n")

# --- Load data ---------------------------------------------------------------
if (!file.exists(r_file)) stop(paste("R results file not found:", r_file))
if (!file.exists(cpp_file)) stop(paste("C++ results file not found:", cpp_file))

r_data   <- read.table(r_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cpp_data <- read.table(cpp_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("R results:   ", nrow(r_data), "markers\n")
cat("C++ results: ", nrow(cpp_data), "markers\n")

# --- Merge on MarkerID -------------------------------------------------------
merged <- merge(r_data, cpp_data, by = "MarkerID", suffixes = c(".R", ".Cpp"))
cat("Overlapping:  ", nrow(merged), "markers\n\n")

if (nrow(merged) == 0) {
  cat("ERROR: No overlapping markers found. Check MarkerID formats.\n")
  cat("R MarkerIDs (first 5):  ", head(r_data$MarkerID, 5), "\n")
  cat("C++ MarkerIDs (first 5):", head(cpp_data$MarkerID, 5), "\n")
  stop("No overlapping markers")
}

# Check for markers only in one output
r_only <- setdiff(r_data$MarkerID, cpp_data$MarkerID)
cpp_only <- setdiff(cpp_data$MarkerID, r_data$MarkerID)
if (length(r_only) > 0) {
  cat("  Markers only in R:   ", length(r_only), "\n")
  cat("    First 5:", head(r_only, 5), "\n")
}
if (length(cpp_only) > 0) {
  cat("  Markers only in C++: ", length(cpp_only), "\n")
  cat("    First 5:", head(cpp_only, 5), "\n")
}

# --- Also load Test 5 (sparse, no noadjCov) to compare difference -----------
test5_r_file <- file.path(base_dir, "test/R/r_step2_sparse_results.txt")
test5_cpp_file <- file.path(base_dir, "test/output/cpp_compare_sparse_results.txt")
has_test5 <- file.exists(test5_r_file) && file.exists(test5_cpp_file)

if (has_test5) {
  test5_r <- read.table(test5_r_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  test5_merged <- merge(r_data, test5_r, by = "MarkerID", suffixes = c(".noadj", ".orig"))
  cat("\n--- Comparing noadjCov vs original R results ---\n")
  cat("Markers in both:", nrow(test5_merged), "\n")
  if (nrow(test5_merged) > 0 && "p.value.noadj" %in% names(test5_merged) && "p.value.orig" %in% names(test5_merged)) {
    pv_noadj <- as.numeric(test5_merged$p.value.noadj)
    pv_orig <- as.numeric(test5_merged$p.value.orig)
    valid <- !is.na(pv_noadj) & !is.na(pv_orig)
    n_differ <- sum(abs(pv_noadj[valid] - pv_orig[valid]) > 1e-12)
    cat("  Markers with different p-values (noadjCov vs original):", n_differ, "out of", sum(valid), "\n")
    if (n_differ > 0) {
      cat("  This confirms the scoreTestFast_noadjCov path was exercised.\n")
    } else {
      cat("  WARNING: No markers differ. scoreTestFast_noadjCov may not have been exercised.\n")
    }
  }
  cat("\n")
}

# --- Compare numeric columns -------------------------------------------------
compare_cols <- c("BETA", "SE", "Tstat", "var", "p.value")

classify <- function(rel_diff) {
  ifelse(rel_diff < 1e-12, "EXACT",
  ifelse(rel_diff < 1e-8,  "MATCH",
  ifelse(rel_diff < 1e-4,  "CLOSE",
                            "DIFFER")))
}

results <- data.frame(
  Column = character(),
  N_compared = integer(),
  N_EXACT = integer(),
  N_MATCH = integer(),
  N_CLOSE = integer(),
  N_DIFFER = integer(),
  Max_RelDiff = numeric(),
  Max_AbsDiff = numeric(),
  Median_RelDiff = numeric(),
  stringsAsFactors = FALSE
)

all_details <- list()

for (col in compare_cols) {
  r_col <- paste0(col, ".R")
  cpp_col <- paste0(col, ".Cpp")

  if (!(r_col %in% names(merged)) || !(cpp_col %in% names(merged))) {
    cat("  WARNING: Column", col, "not found in both outputs. Skipping.\n")
    next
  }

  r_vals <- as.numeric(merged[[r_col]])
  cpp_vals <- as.numeric(merged[[cpp_col]])

  # Skip NAs
  valid <- !is.na(r_vals) & !is.na(cpp_vals)
  r_v <- r_vals[valid]
  cpp_v <- cpp_vals[valid]

  if (length(r_v) == 0) {
    cat("  WARNING: No valid values for column", col, "\n")
    next
  }

  abs_diff <- abs(r_v - cpp_v)
  denom <- pmax(abs(r_v), abs(cpp_v), 1e-300)
  rel_diff <- abs_diff / denom

  # For p-values, handle very small values
  if (col == "p.value") {
    tiny <- (abs(r_v) < 1e-100) & (abs(cpp_v) < 1e-100)
    if (any(tiny)) {
      rel_diff[tiny] <- 0
    }
  }

  cats <- classify(rel_diff)

  results <- rbind(results, data.frame(
    Column = col,
    N_compared = length(r_v),
    N_EXACT = sum(cats == "EXACT"),
    N_MATCH = sum(cats == "MATCH"),
    N_CLOSE = sum(cats == "CLOSE"),
    N_DIFFER = sum(cats == "DIFFER"),
    Max_RelDiff = max(rel_diff),
    Max_AbsDiff = max(abs_diff),
    Median_RelDiff = median(rel_diff),
    stringsAsFactors = FALSE
  ))

  # Store worst cases for detailed output
  if (any(rel_diff > 1e-8)) {
    worst_idx <- order(rel_diff, decreasing = TRUE)
    worst_idx <- head(worst_idx, 10)  # top 10 worst
    all_details[[col]] <- data.frame(
      MarkerID = merged$MarkerID[valid][worst_idx],
      R_value = r_v[worst_idx],
      Cpp_value = cpp_v[worst_idx],
      RelDiff = rel_diff[worst_idx],
      AbsDiff = abs_diff[worst_idx]
    )
  }
}

# --- Print summary -----------------------------------------------------------
cat("\n--- Column-by-Column Summary ---\n\n")
cat(sprintf("%-10s | %7s | %7s | %7s | %7s | %7s | %12s | %12s | %12s\n",
            "Column", "Total", "EXACT", "MATCH", "CLOSE", "DIFFER", "MaxRelDiff", "MaxAbsDiff", "MedRelDiff"))
cat(paste(rep("-", 105), collapse = ""), "\n")
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("%-10s | %7d | %7d | %7d | %7d | %7d | %12.3e | %12.3e | %12.3e\n",
              r$Column, r$N_compared, r$N_EXACT, r$N_MATCH, r$N_CLOSE, r$N_DIFFER,
              r$Max_RelDiff, r$Max_AbsDiff, r$Median_RelDiff))
}

# Print details for worst cases
if (length(all_details) > 0) {
  cat("\n--- Worst Disagreements (top 10 per column) ---\n")
  for (col in names(all_details)) {
    cat("\n", col, ":\n")
    d <- all_details[[col]]
    for (j in seq_len(nrow(d))) {
      cat(sprintf("  %s: R=%.10e  C++=%.10e  relDiff=%.3e  absDiff=%.3e\n",
                  d$MarkerID[j], d$R_value[j], d$Cpp_value[j], d$RelDiff[j], d$AbsDiff[j]))
    }
  }
}

# --- Overall verdict ---------------------------------------------------------
total_differ <- sum(results$N_DIFFER)
total_compared <- sum(results$N_compared)

cat("\n--- Overall ---\n")
cat(sprintf("Total comparisons: %d across %d columns\n", total_compared, nrow(results)))
cat(sprintf("EXACT  (rel diff < 1e-12): %d\n", sum(results$N_EXACT)))
cat(sprintf("MATCH  (rel diff < 1e-8):  %d\n", sum(results$N_MATCH)))
cat(sprintf("CLOSE  (rel diff < 1e-4):  %d\n", sum(results$N_CLOSE)))
cat(sprintf("DIFFER (rel diff >= 1e-4): %d\n", total_differ))

if (total_differ == 0) {
  cat("\nOVERALL: PASS -- All values match within tolerance.\n")
} else {
  pct_differ <- 100.0 * total_differ / total_compared
  cat(sprintf("\nOVERALL: FAIL -- %d values differ beyond tolerance (%.4f%%).\n",
              total_differ, pct_differ))
}

# --- Write to file -----------------------------------------------------------
sink(out_file)
cat("SAIGE Step 2: R vs C++ Comparison Results (Sparse GRM + noadjCov)\n")
cat(paste("Date:", Sys.time(), "\n"))
cat(paste("R results:", r_file, "\n"))
cat(paste("C++ results:", cpp_file, "\n"))
cat(paste("Model: example_quantitative_sparseGRM.rda (LOCO=FALSE, isFastTest=TRUE, is_noadjCov=TRUE)\n"))
cat(paste(rep("=", 105), collapse = ""), "\n\n")

cat(sprintf("R markers:   %d\n", nrow(r_data)))
cat(sprintf("C++ markers: %d\n", nrow(cpp_data)))
cat(sprintf("Overlapping: %d\n\n", nrow(merged)))

cat(sprintf("%-10s | %7s | %7s | %7s | %7s | %7s | %12s | %12s | %12s\n",
            "Column", "Total", "EXACT", "MATCH", "CLOSE", "DIFFER", "MaxRelDiff", "MaxAbsDiff", "MedRelDiff"))
cat(paste(rep("-", 105), collapse = ""), "\n")
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("%-10s | %7d | %7d | %7d | %7d | %7d | %12.3e | %12.3e | %12.3e\n",
              r$Column, r$N_compared, r$N_EXACT, r$N_MATCH, r$N_CLOSE, r$N_DIFFER,
              r$Max_RelDiff, r$Max_AbsDiff, r$Median_RelDiff))
}

if (length(all_details) > 0) {
  cat("\nWorst Disagreements:\n")
  for (col in names(all_details)) {
    cat("\n", col, ":\n")
    d <- all_details[[col]]
    for (j in seq_len(nrow(d))) {
      cat(sprintf("  %s: R=%.10e  C++=%.10e  relDiff=%.3e\n",
                  d$MarkerID[j], d$R_value[j], d$Cpp_value[j], d$RelDiff[j]))
    }
  }
}

cat(sprintf("\nTotal comparisons: %d\n", total_compared))
cat(sprintf("DIFFER: %d\n", total_differ))
if (total_differ == 0) {
  cat("OVERALL: PASS\n")
} else {
  pct_differ <- 100.0 * total_differ / total_compared
  cat(sprintf("OVERALL: FAIL (%d differ, %.4f%%)\n", total_differ, pct_differ))
}
sink()

cat("\nResults written to:", out_file, "\n")
