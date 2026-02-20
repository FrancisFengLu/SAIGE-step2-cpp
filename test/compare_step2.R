#!/usr/bin/env Rscript
# ==============================================================================
# SAIGE Step 2: R vs C++ Result Comparison
# ==============================================================================
# Compares R SPAGMMATtest output with C++ saige-step2 output.
# Both must have been run with the SAME null model and genotype data.
#
# Usage: Rscript compare_step2.R [r_results] [cpp_results] [output_file]
#   Defaults:
#     r_results:   test/R/r_step2_results.txt
#     cpp_results: test/output/cpp_compare_results.txt
#     output_file: test/RESULTS_step2.txt
# ==============================================================================

options(digits = 15)

# --- Parse arguments ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

script_dir <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) getwd())

base_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
if (!dir.exists(base_dir)) base_dir <- getwd()

r_file   <- if (length(args) >= 1) args[1] else file.path(base_dir, "test", "R", "r_step2_results.txt")
cpp_file <- if (length(args) >= 2) args[2] else file.path(base_dir, "test", "output", "cpp_compare_results.txt")
out_file <- if (length(args) >= 3) args[3] else file.path(base_dir, "test", "RESULTS_step2.txt")

cat("=== SAIGE Step 2: R vs C++ Comparison ===\n\n")
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

# --- Compare numeric columns -------------------------------------------------
# Columns to compare (present in both outputs)
compare_cols <- c("BETA", "SE", "Tstat", "var", "p.value")

# Find matching column names (R output may use slightly different names)
find_col <- function(data, target, suffix) {
  candidates <- c(
    paste0(target, suffix),
    paste0(target, ".", suffix),
    target
  )
  for (c in candidates) {
    if (c %in% names(data)) return(c)
  }
  return(NA)
}

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

  # Relative difference: use max(|r|, |cpp|, 1e-300) as denominator
  denom <- pmax(abs(r_v), abs(cpp_v), 1e-300)
  rel_diff <- abs_diff / denom

  # For p-values, also handle the case where both are very small
  if (col == "p.value") {
    # For very small p-values (< 1e-100), use log-scale comparison
    tiny <- (abs(r_v) < 1e-100) & (abs(cpp_v) < 1e-100)
    if (any(tiny)) {
      # Both essentially zero — mark as EXACT
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
  if (any(cats == "DIFFER")) {
    worst_idx <- which(cats == "DIFFER")
    worst_idx <- worst_idx[order(rel_diff[worst_idx], decreasing = TRUE)]
    worst_idx <- head(worst_idx, 5)  # top 5 worst
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
cat(sprintf("%-10s | %7s | %7s | %7s | %7s | %7s | %12s | %12s\n",
            "Column", "Total", "EXACT", "MATCH", "CLOSE", "DIFFER", "MaxRelDiff", "MaxAbsDiff"))
cat(paste(rep("-", 90), collapse = ""), "\n")
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("%-10s | %7d | %7d | %7d | %7d | %7d | %12.3e | %12.3e\n",
              r$Column, r$N_compared, r$N_EXACT, r$N_MATCH, r$N_CLOSE, r$N_DIFFER,
              r$Max_RelDiff, r$Max_AbsDiff))
}

# Print details for DIFFER cases
if (length(all_details) > 0) {
  cat("\n--- Worst Disagreements (DIFFER cases) ---\n")
  for (col in names(all_details)) {
    cat("\n", col, ":\n")
    d <- all_details[[col]]
    for (j in seq_len(nrow(d))) {
      cat(sprintf("  %s: R=%.10e  C++=%.10e  relDiff=%.3e\n",
                  d$MarkerID[j], d$R_value[j], d$Cpp_value[j], d$RelDiff[j]))
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
  cat("\nOVERALL: PASS — All values match within tolerance.\n")
} else {
  cat(sprintf("\nOVERALL: FAIL — %d values differ beyond tolerance.\n", total_differ))
}

# --- Write to file -----------------------------------------------------------
sink(out_file)
cat("SAIGE Step 2: R vs C++ Comparison Results\n")
cat(paste("Date:", Sys.time(), "\n"))
cat(paste("R results:", r_file, "\n"))
cat(paste("C++ results:", cpp_file, "\n"))
cat(paste(rep("=", 90), collapse = ""), "\n\n")

cat(sprintf("%-10s | %7s | %7s | %7s | %7s | %7s | %12s | %12s\n",
            "Column", "Total", "EXACT", "MATCH", "CLOSE", "DIFFER", "MaxRelDiff", "MaxAbsDiff"))
cat(paste(rep("-", 90), collapse = ""), "\n")
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("%-10s | %7d | %7d | %7d | %7d | %7d | %12.3e | %12.3e\n",
              r$Column, r$N_compared, r$N_EXACT, r$N_MATCH, r$N_CLOSE, r$N_DIFFER,
              r$Max_RelDiff, r$Max_AbsDiff))
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
  cat(sprintf("OVERALL: FAIL (%d differ)\n", total_differ))
}
sink()

cat("\nResults written to:", out_file, "\n")
