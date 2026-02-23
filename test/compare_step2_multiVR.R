#!/usr/bin/env Rscript
# ==============================================================================
# Test 8: Multi-Category VR -- R vs C++ Comparison
# ==============================================================================
# Compares R and C++ output for multi-category variance ratio test.
# Additionally verifies that results differ from single-VR run (proving
# multi-VR lookup is actually exercised).
#
# Usage:
#   Rscript compare_step2_multiVR.R
#
# Output:
#   test/RESULTS_step2_multiVR.txt
# ==============================================================================

options(digits = 15)

base_dir <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Step_2_Feb_11"

r_file      <- file.path(base_dir, "test/output/r_multiVR_results.txt")
cpp_file    <- file.path(base_dir, "test/output/cpp_multiVR_results.txt")
singleVR_file <- file.path(base_dir, "test/output/cpp_compare_results.txt")  # single-VR baseline
out_file    <- file.path(base_dir, "test/RESULTS_step2_multiVR.txt")

cat("=== Test 8: Multi-Category VR -- R vs C++ Comparison ===\n\n")
cat("R results (multi-VR):   ", r_file, "\n")
cat("C++ results (multi-VR): ", cpp_file, "\n")
cat("C++ baseline (single-VR):", singleVR_file, "\n")
cat("Output:                  ", out_file, "\n\n")

# --- Load data ---------------------------------------------------------------
if (!file.exists(r_file)) stop("R results not found: ", r_file)
if (!file.exists(cpp_file)) stop("C++ results not found: ", cpp_file)

r_data   <- read.table(r_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cpp_data <- read.table(cpp_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("R results:   ", nrow(r_data), "markers\n")
cat("C++ results: ", nrow(cpp_data), "markers\n")

# Merge on MarkerID
merged <- merge(r_data, cpp_data, by = "MarkerID", suffixes = c(".R", ".Cpp"))
cat("Overlapping: ", nrow(merged), "markers\n\n")

if (nrow(merged) == 0) {
  stop("No overlapping markers found.")
}

# --- MAC category analysis ----------------------------------------------------
# Compute MAC for each marker
mac_r <- pmin(merged$AC_Allele2.R, 2 * merged$N.R - merged$AC_Allele2.R)

n_lowMAC  <- sum(mac_r <= 10.5, na.rm = TRUE)
n_midMAC  <- sum(mac_r > 10.5 & mac_r <= 20.5, na.rm = TRUE)
n_highMAC <- sum(mac_r > 20.5, na.rm = TRUE)

cat("MAC category counts:\n")
cat("  MAC <= 10.5        (fallback, VR=1.15):  ", n_lowMAC, "markers\n")
cat("  10.5 < MAC <= 20.5 (cat 1, VR=1.15):     ", n_midMAC, "markers\n")
cat("  MAC > 20.5          (cat 2, VR=1.22420): ", n_highMAC, "markers\n\n")

# --- Column-by-column comparison ----------------------------------------------
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

  valid <- !is.na(r_vals) & !is.na(cpp_vals)
  r_v <- r_vals[valid]
  cpp_v <- cpp_vals[valid]

  if (length(r_v) == 0) next

  abs_diff <- abs(r_v - cpp_v)
  denom <- pmax(abs(r_v), abs(cpp_v), 1e-300)
  rel_diff <- abs_diff / denom

  if (col == "p.value") {
    tiny <- (abs(r_v) < 1e-100) & (abs(cpp_v) < 1e-100)
    if (any(tiny)) rel_diff[tiny] <- 0
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

  if (any(cats == "DIFFER")) {
    worst_idx <- which(cats == "DIFFER")
    worst_idx <- worst_idx[order(rel_diff[worst_idx], decreasing = TRUE)]
    worst_idx <- head(worst_idx, 5)
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
cat("\n--- Column-by-Column Summary (R multiVR vs C++ multiVR) ---\n\n")
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

# --- Verify multi-VR differs from single-VR ----------------------------------
multiVR_differs_from_single <- FALSE

if (file.exists(singleVR_file)) {
  cat("\n--- Comparing multi-VR vs single-VR (sanity check) ---\n")
  single_data <- read.table(singleVR_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Merge multi-VR C++ with single-VR C++
  merged_vs <- merge(cpp_data, single_data, by = "MarkerID", suffixes = c(".Multi", ".Single"))

  if (nrow(merged_vs) > 0) {
    # Compare var column (variance is directly scaled by VR)
    var_multi  <- as.numeric(merged_vs$var.Multi)
    var_single <- as.numeric(merged_vs$var.Single)
    valid2 <- !is.na(var_multi) & !is.na(var_single)

    if (sum(valid2) > 0) {
      var_diff <- abs(var_multi[valid2] - var_single[valid2])
      n_differ <- sum(var_diff > 1e-10)

      # Also check by MAC category
      mac_vs <- pmin(merged_vs$AC_Allele2.Multi, 2 * merged_vs$N.Multi - merged_vs$AC_Allele2.Multi)

      low_mac_idx   <- which(valid2 & mac_vs <= 10.5)
      mid_mac_idx   <- which(valid2 & mac_vs > 10.5 & mac_vs <= 20.5)
      high_mac_idx  <- which(valid2 & mac_vs > 20.5)

      cat("  Multi-VR vs Single-VR 'var' differences:\n")
      cat("    Total markers compared:", sum(valid2), "\n")
      cat("    Markers with var diff > 1e-10:", n_differ, "\n")

      if (length(low_mac_idx) > 0) {
        n_diff_low <- sum(var_diff[low_mac_idx - which(valid2)[1] + 1] > 1e-10, na.rm = TRUE)
        cat("    MAC <= 10.5 with diff:         ", n_diff_low, "/", length(low_mac_idx), "\n")
      }
      if (length(mid_mac_idx) > 0) {
        n_diff_mid <- sum(abs(var_multi[mid_mac_idx] - var_single[mid_mac_idx]) > 1e-10, na.rm = TRUE)
        cat("    10.5 < MAC <= 20.5 with diff:  ", n_diff_mid, "/", length(mid_mac_idx), "\n")
      }
      if (length(high_mac_idx) > 0) {
        n_diff_high <- sum(abs(var_multi[high_mac_idx] - var_single[high_mac_idx]) > 1e-10, na.rm = TRUE)
        cat("    MAC > 20.5 with diff:          ", n_diff_high, "/", length(high_mac_idx), "\n")
      }

      if (n_differ > 0) {
        cat("  CONFIRMED: Multi-VR produces different results from single-VR.\n")
        multiVR_differs_from_single <- TRUE
      } else {
        cat("  WARNING: Multi-VR results identical to single-VR. MAC lookup may not be exercised.\n")
      }
    }
  }
} else {
  cat("\n  (Single-VR baseline file not found, skipping sanity check)\n")
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
  cat("\nOVERALL: PASS -- All R vs C++ values match within tolerance.\n")
  if (multiVR_differs_from_single) {
    cat("SANITY:  PASS -- Multi-VR results differ from single-VR (MAC lookup exercised).\n")
  }
} else {
  cat(sprintf("\nOVERALL: FAIL -- %d values differ beyond tolerance.\n", total_differ))
}

# --- Write to file -----------------------------------------------------------
sink(out_file)
cat("==========================================================================\n")
cat("Test 8: Multi-Category Variance Ratio -- R vs C++ Comparison Results\n")
cat("==========================================================================\n")
cat(paste("Date:", Sys.time(), "\n"))
cat(paste("R results (multi-VR):   ", r_file, "\n"))
cat(paste("C++ results (multi-VR): ", cpp_file, "\n"))
cat(paste("C++ baseline (single-VR):", singleVR_file, "\n"))
cat(paste(rep("=", 75), collapse = ""), "\n\n")

cat("VR file: test/data/multiVR_varianceRatio.txt\n")
cat("  Category 1: MAC in (10.5, 20.5] -> VR = 1.15\n")
cat("  Category 2: MAC in (20.5, 1000] -> VR = 1.22420115246319\n")
cat("  Fallback:   MAC <= 10.5         -> VR = 1.15\n\n")

cat("MAC category counts:\n")
cat(sprintf("  MAC <= 10.5:         %d markers\n", n_lowMAC))
cat(sprintf("  10.5 < MAC <= 20.5:  %d markers\n", n_midMAC))
cat(sprintf("  MAC > 20.5:          %d markers\n\n", n_highMAC))

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
  cat("\nOVERALL: PASS\n")
} else {
  cat(sprintf("\nOVERALL: FAIL (%d differ)\n", total_differ))
}

if (multiVR_differs_from_single) {
  cat("SANITY:  PASS (multi-VR differs from single-VR)\n")
} else if (file.exists(singleVR_file)) {
  cat("SANITY:  WARNING (multi-VR identical to single-VR -- MAC lookup may not be exercised)\n")
}
sink()

cat("\nResults written to:", out_file, "\n")
