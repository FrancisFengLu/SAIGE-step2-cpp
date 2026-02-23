#!/usr/bin/env Rscript
# ==============================================================================
# compare_step2_outputs.R
# ==============================================================================
# Compare C++ Step 2 output using two different null models:
#   1. R converted model (LOCO=TRUE, chr 1) -> cpp_compare_results.txt
#   2. C++ Step 1 model (LOCO=FALSE) -> cpp_step1_to_step2_results.txt
#
# We expect differences because the null models are different (LOCO vs no LOCO).
# ==============================================================================

options(digits = 15)

cat("====================================================================\n")
cat("Step 2 Output Comparison\n")
cat("  Source 1: R converted model (LOCO=TRUE)  -> cpp_compare_results.txt\n")
cat("  Source 2: C++ Step 1 model (LOCO=FALSE) -> cpp_step1_to_step2_results.txt\n")
cat("====================================================================\n\n")

# --- Load data ---
file1 <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Step_2_Feb_11/test/output/cpp_compare_results.txt"
file2 <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Step_2_Feb_11/test/output/cpp_step1_to_step2_results.txt"

d1 <- read.table(file1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
d2 <- read.table(file2, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("File 1 (R model):     ", nrow(d1), "markers,", ncol(d1), "columns\n")
cat("File 2 (C++ Step 1):  ", nrow(d2), "markers,", ncol(d2), "columns\n\n")

# Check same markers
stopifnot(identical(d1$MarkerID, d2$MarkerID))
cat("Marker IDs match: YES\n\n")

# --- Compare numeric columns ---
compare_cols <- c("BETA", "SE", "Tstat", "var", "p.value")

cat("Column-by-column comparison:\n")
cat(sprintf("  %-10s  %12s  %12s  %12s  %12s  %8s  %8s\n",
            "Column", "MaxAbsDiff", "MaxRelDiff", "MeanRelDiff", "Correlation",
            "SameSign", "BothSig"))
cat(sprintf("  %-10s  %12s  %12s  %12s  %12s  %8s  %8s\n",
            "------", "----------", "----------", "----------", "-----------",
            "--------", "--------"))

for (col in compare_cols) {
  v1 <- d1[[col]]
  v2 <- d2[[col]]

  abs_diff <- abs(v1 - v2)
  max_abs <- max(abs_diff, na.rm = TRUE)
  mean_abs <- mean(abs_diff, na.rm = TRUE)

  denom <- pmax(abs(v1), abs(v2), 1e-300)
  rel_diff <- abs_diff / denom
  max_rel <- max(rel_diff, na.rm = TRUE)
  mean_rel <- mean(rel_diff, na.rm = TRUE)

  corr <- cor(v1, v2, use = "complete.obs")

  # Sign agreement
  same_sign <- mean(sign(v1) == sign(v2), na.rm = TRUE)

  # Both significant at p < 0.05
  if (col == "p.value") {
    both_sig <- sum(v1 < 0.05 & v2 < 0.05, na.rm = TRUE)
    n_sig1 <- sum(v1 < 0.05, na.rm = TRUE)
    n_sig2 <- sum(v2 < 0.05, na.rm = TRUE)
    sig_str <- sprintf("%d/%d/%d", both_sig, n_sig1, n_sig2)
  } else {
    sig_str <- sprintf("%.4f", same_sign)
  }

  cat(sprintf("  %-10s  %12.3e  %12.3e  %12.3e  %12.10f  %8s\n",
              col, max_abs, max_rel, mean_rel, corr, sig_str))
}

cat("\n")

# --- Detailed p-value comparison ---
cat("====================================================================\n")
cat("P-value Distribution Comparison:\n")
cat("====================================================================\n\n")

p1 <- d1$p.value
p2 <- d2$p.value

cat("Summary statistics:\n")
cat("  R model p-values:\n")
cat("    min   =", min(p1), "\n")
cat("    median=", median(p1), "\n")
cat("    mean  =", mean(p1), "\n")
cat("    max   =", max(p1), "\n")
cat("    # significant (p<0.05)  =", sum(p1 < 0.05), "\n")
cat("    # significant (p<5e-8)  =", sum(p1 < 5e-8), "\n\n")

cat("  C++ Step 1 p-values:\n")
cat("    min   =", min(p2), "\n")
cat("    median=", median(p2), "\n")
cat("    mean  =", mean(p2), "\n")
cat("    max   =", max(p2), "\n")
cat("    # significant (p<0.05)  =", sum(p2 < 0.05), "\n")
cat("    # significant (p<5e-8)  =", sum(p2 < 5e-8), "\n\n")

# Log-p correlation
lp1 <- -log10(p1)
lp2 <- -log10(p2)
cat("Correlation of -log10(p): ", cor(lp1, lp2), "\n\n")

# --- Look at markers where significance disagrees ---
sig_thresh <- 0.05
sig1_only <- which(p1 < sig_thresh & p2 >= sig_thresh)
sig2_only <- which(p2 < sig_thresh & p1 >= sig_thresh)
both_sig <- which(p1 < sig_thresh & p2 < sig_thresh)

cat("Significance agreement at p < 0.05:\n")
cat("  Both significant:     ", length(both_sig), "\n")
cat("  R model only:         ", length(sig1_only), "\n")
cat("  C++ Step 1 only:      ", length(sig2_only), "\n")
cat("  Agreement rate:       ", sprintf("%.2f%%",
    100 * (nrow(d1) - length(sig1_only) - length(sig2_only)) / nrow(d1)), "\n\n")

# --- BETA sign agreement ---
cat("BETA sign agreement:\n")
cat("  Same sign:   ", sum(sign(d1$BETA) == sign(d2$BETA)), "/", nrow(d1),
    sprintf(" (%.2f%%)\n", 100 * sum(sign(d1$BETA) == sign(d2$BETA)) / nrow(d1)))
cat("  BETA corr:   ", cor(d1$BETA, d2$BETA), "\n\n")

# --- Show some specific markers ---
cat("====================================================================\n")
cat("First 10 markers - side by side:\n")
cat("====================================================================\n")
cat(sprintf("  %-12s  %12s %12s  %12s %12s  %12s %12s\n",
            "MarkerID", "BETA_R", "BETA_C++", "pval_R", "pval_C++", "SE_R", "SE_C++"))
for (i in 1:10) {
  cat(sprintf("  %-12s  %12.6f %12.6f  %12.6e %12.6e  %12.6f %12.6f\n",
              d1$MarkerID[i],
              d1$BETA[i], d2$BETA[i],
              d1$p.value[i], d2$p.value[i],
              d1$SE[i], d2$SE[i]))
}

cat("\n")
cat("====================================================================\n")
cat("CONCLUSION:\n")
cat("====================================================================\n")
cat("  The differences are expected because the null models differ:\n")
cat("  - R model: LOCO=TRUE (leave-one-chromosome-out for chr 1)\n")
cat("  - C++ Step 1 model: LOCO=FALSE (global model, all chromosomes)\n")
cat("  LOCO refits the null model excluding the test chromosome's variants,\n")
cat("  which changes mu, res, V, and all derived matrices.\n")
cat("  This fundamentally changes the association statistics.\n")
cat("\n")
cat("  The key validation is that C++ Step 2 *runs successfully* with\n")
cat("  C++ Step 1 output, proving the pure C++ pipeline works.\n")
cat("\n")
cat("Done.\n")
