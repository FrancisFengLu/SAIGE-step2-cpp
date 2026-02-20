#!/usr/bin/env Rscript
# ==============================================================================
# SAIGE Step 2: R vs C++ Region-Based Result Comparison
# ==============================================================================
# Compares R SPAGMMATtest region output with C++ saige-step2 region output.
# Both must have been run with the SAME null model, genotype data, and group file.
#
# Usage: Rscript compare_step2_region.R [r_results] [cpp_results] [output_file]
#   Defaults:
#     r_results:   test/R/r_step2_region_results.txt
#     cpp_results: test/output/cpp_compare_region_results.txt
#     output_file: test/RESULTS_step2_region.txt
# ==============================================================================

options(digits = 15)

# --- Parse arguments ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

script_dir <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) getwd())

base_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
if (!dir.exists(base_dir)) base_dir <- getwd()

r_file   <- if (length(args) >= 1) args[1] else file.path(base_dir, "test", "R", "r_step2_region_results.txt")
cpp_file <- if (length(args) >= 2) args[2] else file.path(base_dir, "test", "output", "cpp_compare_region_results.txt")
out_file <- if (length(args) >= 3) args[3] else file.path(base_dir, "test", "RESULTS_step2_region.txt")

cat("=== SAIGE Step 2: R vs C++ Region Comparison ===\n\n")
cat("R results:   ", r_file, "\n")
cat("C++ results: ", cpp_file, "\n")
cat("Output:      ", out_file, "\n\n")

# --- Load data ---------------------------------------------------------------
if (!file.exists(r_file)) stop(paste("R results file not found:", r_file))
if (!file.exists(cpp_file)) stop(paste("C++ results file not found:", cpp_file))

r_data   <- read.table(r_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cpp_data <- read.table(cpp_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("R results:   ", nrow(r_data), "rows\n")
cat("C++ results: ", nrow(cpp_data), "rows\n")

# Show column names
cat("R columns:   ", paste(names(r_data), collapse = ", "), "\n")
cat("C++ columns: ", paste(names(cpp_data), collapse = ", "), "\n\n")

# --- Merge on Region + Group + max_MAF (the unique key for region results) ---
# Both should have: Region, Group, max_MAF as key columns

# Create a merge key
r_data$merge_key <- paste(r_data$Region, r_data$Group, r_data$max_MAF, sep = "___")
cpp_data$merge_key <- paste(cpp_data$Region, cpp_data$Group, cpp_data$max_MAF, sep = "___")

merged <- merge(r_data, cpp_data, by = "merge_key", suffixes = c(".R", ".Cpp"))
cat("Overlapping:  ", nrow(merged), "rows\n\n")

if (nrow(merged) == 0) {
  cat("ERROR: No overlapping rows found. Check Region/Group/max_MAF formats.\n")
  cat("R keys (first 5):  ", head(r_data$merge_key, 5), "\n")
  cat("C++ keys (first 5):", head(cpp_data$merge_key, 5), "\n")
  stop("No overlapping rows")
}

# --- Compare numeric columns -------------------------------------------------
# Region output columns to compare
compare_cols <- c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT", "BETA_Burden", "SE_Burden",
                  "MAC", "Number_rare", "Number_ultra_rare")

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

  r_vals <- suppressWarnings(as.numeric(merged[[r_col]]))
  cpp_vals <- suppressWarnings(as.numeric(merged[[cpp_col]]))

  # Skip NAs (Cauchy rows have NA for some columns)
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

  # For p-values, handle very small values
  if (grepl("Pvalue|pval", col, ignore.case = TRUE)) {
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
  if (any(cats %in% c("DIFFER", "CLOSE"))) {
    worst_idx <- which(cats %in% c("DIFFER", "CLOSE"))
    worst_idx <- worst_idx[order(rel_diff[worst_idx], decreasing = TRUE)]
    worst_idx <- head(worst_idx, 5)  # top 5 worst
    all_details[[col]] <- data.frame(
      Key = merged$merge_key[valid][worst_idx],
      R_value = r_v[worst_idx],
      Cpp_value = cpp_v[worst_idx],
      RelDiff = rel_diff[worst_idx],
      AbsDiff = abs_diff[worst_idx]
    )
  }
}

# --- Print summary -----------------------------------------------------------
cat("\n--- Column-by-Column Summary ---\n\n")
cat(sprintf("%-16s | %7s | %7s | %7s | %7s | %7s | %12s | %12s\n",
            "Column", "Total", "EXACT", "MATCH", "CLOSE", "DIFFER", "MaxRelDiff", "MaxAbsDiff"))
cat(paste(rep("-", 100), collapse = ""), "\n")
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("%-16s | %7d | %7d | %7d | %7d | %7d | %12.3e | %12.3e\n",
              r$Column, r$N_compared, r$N_EXACT, r$N_MATCH, r$N_CLOSE, r$N_DIFFER,
              r$Max_RelDiff, r$Max_AbsDiff))
}

# Print details for DIFFER/CLOSE cases
if (length(all_details) > 0) {
  cat("\n--- Worst Disagreements ---\n")
  for (col in names(all_details)) {
    cat("\n", col, ":\n")
    d <- all_details[[col]]
    for (j in seq_len(nrow(d))) {
      cat(sprintf("  %s: R=%.10e  C++=%.10e  relDiff=%.3e\n",
                  d$Key[j], d$R_value[j], d$Cpp_value[j], d$RelDiff[j]))
    }
  }
}

# --- Side-by-side comparison --------------------------------------------------
cat("\n--- Side-by-Side Results ---\n\n")
# Print each row showing R and C++ values
for (i in seq_len(nrow(merged))) {
  row <- merged[i, ]
  cat(sprintf("Row %d: %s\n", i, row$merge_key))
  for (col in compare_cols) {
    r_col <- paste0(col, ".R")
    cpp_col <- paste0(col, ".Cpp")
    if (r_col %in% names(merged) && cpp_col %in% names(merged)) {
      r_val <- row[[r_col]]
      cpp_val <- row[[cpp_col]]
      if (!is.na(r_val) && !is.na(cpp_val)) {
        r_num <- suppressWarnings(as.numeric(r_val))
        cpp_num <- suppressWarnings(as.numeric(cpp_val))
        if (!is.na(r_num) && !is.na(cpp_num)) {
          denom <- max(abs(r_num), abs(cpp_num), 1e-300)
          rd <- abs(r_num - cpp_num) / denom
          cat(sprintf("  %-16s: R=%-20.10e  C++=%-20.10e  relDiff=%.3e\n",
                      col, r_num, cpp_num, rd))
        }
      }
    }
  }
  cat("\n")
}

# --- Overall verdict ---------------------------------------------------------
total_differ <- sum(results$N_DIFFER)
total_compared <- sum(results$N_compared)

cat("--- Overall ---\n")
cat(sprintf("Total comparisons: %d across %d columns\n", total_compared, nrow(results)))
cat(sprintf("EXACT  (rel diff < 1e-12): %d\n", sum(results$N_EXACT)))
cat(sprintf("MATCH  (rel diff < 1e-8):  %d\n", sum(results$N_MATCH)))
cat(sprintf("CLOSE  (rel diff < 1e-4):  %d\n", sum(results$N_CLOSE)))
cat(sprintf("DIFFER (rel diff >= 1e-4): %d\n", total_differ))

if (total_differ == 0) {
  cat("\nOVERALL: PASS -- All values match within tolerance.\n")
} else {
  cat(sprintf("\nOVERALL: FAIL -- %d values differ beyond tolerance.\n", total_differ))
}

# --- Write to file -----------------------------------------------------------
sink(out_file)
cat("SAIGE Step 2: R vs C++ Region Comparison Results\n")
cat(paste("Date:", Sys.time(), "\n"))
cat(paste("R results:", r_file, "\n"))
cat(paste("C++ results:", cpp_file, "\n"))
cat(paste(rep("=", 100), collapse = ""), "\n\n")

cat(sprintf("%-16s | %7s | %7s | %7s | %7s | %7s | %12s | %12s\n",
            "Column", "Total", "EXACT", "MATCH", "CLOSE", "DIFFER", "MaxRelDiff", "MaxAbsDiff"))
cat(paste(rep("-", 100), collapse = ""), "\n")
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("%-16s | %7d | %7d | %7d | %7d | %7d | %12.3e | %12.3e\n",
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
                  d$Key[j], d$R_value[j], d$Cpp_value[j], d$RelDiff[j]))
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
