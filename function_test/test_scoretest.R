# R-side score test computation for comparison with C++ SAIGEClass::scoreTestFast
# This script manually computes the score test using the same .arma matrices
# that the C++ code uses, following the exact same mathematical formulas.

options(digits = 15)

cat("=== R Score Test Comparison ===\n\n")

# --- Helper: Load Armadillo binary file ---
# Armadillo binary format: ARMA_MAT_BIN_FN008 or ARMA_COL_BIN_FN008
# Line 1 (text): header type
# Line 2 (text): "rows cols"
# Remaining: raw binary doubles (column-major)
load_arma_bin <- function(filepath) {
  con <- file(filepath, "rb")
  on.exit(close(con))

  # Read header line (text, newline-terminated)
  header <- readLines(con, n = 1)
  dims_line <- readLines(con, n = 1)
  dims <- as.integer(strsplit(trimws(dims_line), "\\s+")[[1]])

  nrow <- dims[1]
  ncol <- dims[2]

  # Read remaining binary data (doubles, 8 bytes each)
  n_elements <- nrow * ncol
  values <- readBin(con, what = "double", n = n_elements, size = 8, endian = "little")

  if (ncol == 1) {
    return(as.numeric(values))
  } else {
    # Armadillo stores column-major
    mat <- matrix(values, nrow = nrow, ncol = ncol, byrow = FALSE)
    return(mat)
  }
}

# --- Load .arma files ---
prefix <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/output/cpp_dense_x1_output"

cat("Loading .arma files...\n")
mu_vec <- load_arma_bin(paste0(prefix, ".mu.arma"))
res_vec <- load_arma_bin(paste0(prefix, ".res.arma"))
y_vec <- load_arma_bin(paste0(prefix, ".y.arma"))
V_vec <- load_arma_bin(paste0(prefix, ".V.arma"))
S_a <- load_arma_bin(paste0(prefix, ".S_a.arma"))
X_mat <- load_arma_bin(paste0(prefix, ".X.arma"))
XVX_mat <- load_arma_bin(paste0(prefix, ".XVX.arma"))
XVX_inv_mat <- load_arma_bin(paste0(prefix, ".XVX_inv.arma"))
XXVX_inv_mat <- load_arma_bin(paste0(prefix, ".XXVX_inv.arma"))
XV_mat <- load_arma_bin(paste0(prefix, ".XV.arma"))
XVX_inv_XV_mat <- load_arma_bin(paste0(prefix, ".XVX_inv_XV.arma"))

n <- length(y_vec)

cat("n =", n, "\n")
cat("X dims:", dim(X_mat), "\n")
cat("XV dims:", dim(XV_mat), "\n")
cat("XVX dims:", dim(XVX_mat), "\n")
cat("XXVX_inv dims:", dim(XXVX_inv_mat), "\n")
cat("XVX_inv_XV dims:", dim(XVX_inv_XV_mat), "\n")
cat("S_a:", S_a, "\n")
cat("mu[1:5]:", mu_vec[1:5], "\n")
cat("res[1:5]:", res_vec[1:5], "\n")
cat("y[1:5]:", y_vec[1:5], "\n")

# --- Trait parameters ---
# From nullmodel.json: trait=quantitative, theta=[0.286165, 0.458186]
tau0 <- 0.286165
tau1 <- 0.458186
tauvec <- c(tau0, tau1)
traitType <- "quantitative"

# For quantitative trait: mu2 = 1/tau[1] * ones(N)
mu2 <- rep(1.0 / tau1, n)

# Variance ratio
varRatio <- 1.1105

cat("tau0 =", tau0, "\n")
cat("tau1 =", tau1, "\n")
cat("varRatio =", varRatio, "\n\n")

# --- Create test genotype vectors (same as C++) ---
# C++ uses 0-indexed: G1(0)=1, G1(5)=1, G1(10)=1
# R is 1-indexed: G1[1]=1, G1[6]=1, G1[11]=1
G1 <- rep(0.0, n)
G1[1] <- 1.0; G1[6] <- 1.0; G1[11] <- 1.0

G2 <- rep(0.0, n)
G2[1:100] <- 1.0

G3 <- rep(0.0, n)
G3[seq(2, n, by = 2)] <- 1.0  # indices 2,4,6,...,1000 = C++ indices 1,3,5,...,999

gnames <- c("G1_sparse", "G2_block100", "G3_alternating")
gvecs <- list(G1, G2, G3)

# --- scoreTestFast implementation in R ---
# Following the exact C++ code in saige_test.cpp lines 205-286
# For quantitative trait path (lines 234-237):
#   var2 = ZtXVXZ(0,0)*m_tauvec[0] + dot(g1,g1) - 2*Bmu2
#   where Bmu2 = dot(g1, B) for quantitative
#
# Note: m_XVX_inv_XV in C++ is .rows(t_indexForNonZero) - it's N x p
# So A1 = XVX_inv_XV[nonzero, ] which is |nz| x p
# Z = t(A1) %*% g1 = p x 1
# B = X1 %*% Z = |nz| x 1

scoreTestFast_R <- function(GVec, indexNonZero, m_X, m_XVX_inv_XV, m_XVX, m_res, m_S_a, m_tauvec, m_varRatioVal, traitType, mu2) {
  # indexNonZero is 1-indexed in R
  g1 <- GVec[indexNonZero]
  X1 <- m_X[indexNonZero, , drop = FALSE]
  A1 <- m_XVX_inv_XV[indexNonZero, , drop = FALSE]
  res1 <- m_res[indexNonZero]

  Z <- as.vector(t(A1) %*% g1)   # p x 1
  B <- as.vector(X1 %*% Z)       # |nz| x 1
  g1_tilde <- g1 - B

  ZtXVXZ <- as.numeric(t(Z) %*% m_XVX %*% Z)

  if (traitType == "quantitative") {
    Bmu2 <- sum(g1 * B)
    var2 <- ZtXVXZ * m_tauvec[1] + sum(g1^2) - 2 * Bmu2
  } else if (traitType == "binary" || traitType == "survival") {
    mu21 <- mu2[indexNonZero]
    g1tildemu2 <- sum(g1_tilde^2 * mu21)
    Bmu2_val <- sum(B^2 * mu21)
    var2 <- ZtXVXZ - Bmu2_val + g1tildemu2
  }

  var1 <- var2 * m_varRatioVal

  # Score statistic
  S1 <- sum(res1 * g1_tilde)
  res1X1 <- as.vector(t(res1) %*% X1)
  S_a2 <- m_S_a - res1X1
  S2 <- -sum(S_a2 * Z)
  S <- S1 + S2
  S <- S / m_tauvec[1]

  stat <- S^2 / var1

  if (var1 <= .Machine$double.xmin) {
    pval <- 1.0
  } else {
    if (!is.nan(stat) && is.finite(stat)) {
      pval <- pchisq(stat, df = 1, lower.tail = FALSE)
    } else {
      pval <- 1.0
      stat <- 0.0
    }
  }

  Beta <- S / var1
  seBeta <- abs(Beta) / sqrt(abs(stat))

  return(list(
    Tstat = S,
    var1 = var1,
    var2 = var2,
    Beta = Beta,
    seBeta = seBeta,
    pval = pval,
    stat = stat
  ))
}

# --- Also implement the non-fast scoreTest for comparison ---
# Following saige_test.cpp lines 125-202
scoreTest_R <- function(GVec, indexNonZero, m_X, m_XXVX_inv, m_XV, m_XVX_inv_XV, m_res, m_mu2, m_tauvec, m_varRatioVal, traitType) {
  p <- nrow(m_XV)

  # getadjGFast: compute adjusted genotype
  m_XVG <- rep(0.0, p)
  for (i in seq_along(indexNonZero)) {
    idx <- indexNonZero[i]
    m_XVG <- m_XVG + m_XV[, idx] * GVec[idx]
  }
  gtilde <- GVec - as.vector(m_XXVX_inv %*% m_XVG)

  S <- sum(gtilde * m_res)
  S <- S / m_tauvec[1]

  # For non-sparse GRM path:
  # t_P2Vec = t_gtilde % m_mu2 *m_tauvec[0]
  P2Vec <- gtilde * m_mu2 * m_tauvec[1]
  var2 <- sum(P2Vec * gtilde)
  var1 <- var2 * m_varRatioVal

  stat <- S^2 / var1

  if (var1 <= .Machine$double.xmin) {
    pval <- 1.0
  } else {
    if (!is.nan(stat) && is.finite(stat)) {
      pval <- pchisq(stat, df = 1, lower.tail = FALSE)
    } else {
      pval <- 1.0
      stat <- 0.0
    }
  }

  Beta <- S / var1
  seBeta <- abs(Beta) / sqrt(abs(stat))

  return(list(
    Tstat = S,
    var1 = var1,
    var2 = var2,
    Beta = Beta,
    seBeta = seBeta,
    pval = pval,
    stat = stat
  ))
}


# --- Run tests ---
cat("=== scoreTestFast results ===\n\n")

results_fast <- list()
for (gi in seq_along(gvecs)) {
  G <- gvecs[[gi]]
  indexNonZero <- which(G != 0)  # 1-indexed in R

  res <- scoreTestFast_R(G, indexNonZero, X_mat, XVX_inv_XV_mat, XVX_mat, res_vec, S_a, tauvec, varRatio, traitType, mu2)

  cat(sprintf("=== %s ===\n", gnames[gi]))
  cat(sprintf("  Tstat (S)   = %.15f\n", res$Tstat))
  cat(sprintf("  var1        = %.15f\n", res$var1))
  cat(sprintf("  var2        = %.15f\n", res$var2))
  cat(sprintf("  Beta        = %.15f\n", res$Beta))
  cat(sprintf("  seBeta      = %.15f\n", res$seBeta))
  cat(sprintf("  pval        = %.15e\n", res$pval))
  cat("\n")

  results_fast[[gi]] <- res
}

# Write results to file
outfile <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Step_2_Feb_11/function_test/r_scoretest_results.txt"
con <- file(outfile, "w")
for (gi in seq_along(gvecs)) {
  res <- results_fast[[gi]]
  writeLines(sprintf("%s\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15e\tNA\t0",
                     gnames[gi], res$Tstat, res$var1, res$var2, res$Beta, res$seBeta, res$pval), con)
}
close(con)
cat("R fast results written to:", outfile, "\n\n")


cat("=== scoreTest (non-fast) results ===\n\n")

results_nonfast <- list()
for (gi in seq_along(gvecs)) {
  G <- gvecs[[gi]]
  indexNonZero <- which(G != 0)

  res <- scoreTest_R(G, indexNonZero, X_mat, XXVX_inv_mat, XV_mat, XVX_inv_XV_mat, res_vec, mu2, tauvec, varRatio, traitType)

  cat(sprintf("=== %s (non-fast) ===\n", gnames[gi]))
  cat(sprintf("  Tstat (S)   = %.15f\n", res$Tstat))
  cat(sprintf("  var1        = %.15f\n", res$var1))
  cat(sprintf("  var2        = %.15f\n", res$var2))
  cat(sprintf("  Beta        = %.15f\n", res$Beta))
  cat(sprintf("  seBeta      = %.15f\n", res$seBeta))
  cat(sprintf("  pval        = %.15e\n", res$pval))
  cat("\n")

  results_nonfast[[gi]] <- res
}

# Write non-fast results
outfile2 <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Step_2_Feb_11/function_test/r_scoretest_nonfastresults.txt"
con2 <- file(outfile2, "w")
for (gi in seq_along(gvecs)) {
  res <- results_nonfast[[gi]]
  writeLines(sprintf("%s\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15e\tNA\t0",
                     gnames[gi], res$Tstat, res$var1, res$var2, res$Beta, res$seBeta, res$pval), con2)
}
close(con2)
cat("R non-fast results written to:", outfile2, "\n\n")

# --- Sanity check: verify fast vs non-fast agree in R ---
cat("=== R internal consistency: fast vs non-fast ===\n")
for (gi in seq_along(gvecs)) {
  rf <- results_fast[[gi]]
  rn <- results_nonfast[[gi]]
  diff_S <- abs(rf$Tstat - rn$Tstat)
  diff_var1 <- abs(rf$var1 - rn$var1)
  diff_var2 <- abs(rf$var2 - rn$var2)
  cat(sprintf("%s: |dS|=%.2e, |dvar1|=%.2e, |dvar2|=%.2e\n", gnames[gi], diff_S, diff_var1, diff_var2))
}
cat("\nDone.\n")
