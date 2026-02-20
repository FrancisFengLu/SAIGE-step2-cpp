#!/usr/bin/env Rscript
# Debug script: extract Score, Phi, eigenvalues, Q statistic for each gene group
# so we can compare R vs C++ SKAT inputs

library(SKAT)

# Load the R region results to understand which groups are being tested
# We need to run the actual region test with debug output

# First, let's understand what SKAT's internal functions look like
cat("=== Get_Davies_PVal source ===\n")
print(SKAT:::Get_Davies_PVal)

cat("\n=== Get_Lambda source ===\n")
print(SKAT:::Get_Lambda)

cat("\n=== SKAT_liu.MOD.Lambda source ===\n")
tryCatch(print(SKAT:::SKAT_liu.MOD.Lambda), error=function(e) cat("Not found:", conditionMessage(e), "\n"))

cat("\n=== Met_SKAT_Get_Pvalue source ===\n")
print(SKAT:::Met_SKAT_Get_Pvalue)

cat("\n=== SKAT_Optimal_Each_Q source ===\n")
print(SKAT:::SKAT_Optimal_Each_Q)

cat("\n=== SKAT_Optimal_PValue_Davies source ===\n")
print(SKAT:::SKAT_Optimal_PValue_Davies)

cat("\n=== SKAT_Optimal_Param source ===\n")
tryCatch(print(SKAT:::SKAT_Optimal_Param), error=function(e) cat("Not found:", conditionMessage(e), "\n"))

cat("\n=== SKAT_Optimal_Integrate_Func_Davies source ===\n")
tryCatch(print(SKAT:::SKAT_Optimal_Integrate_Func_Davies), error=function(e) cat("Not found:", conditionMessage(e), "\n"))

# Try a simple test to verify our understanding of the Davies method
cat("\n\n=== Simple Davies test ===\n")
# Create a simple test case: Q = sum(lambda * chi2(1))
# lambda = c(1, 2, 3), Q = 5
# Get_Davies_PVal expects: Q, lambda
tryCatch({
  result <- SKAT:::Get_Davies_PVal(5, c(1, 2, 3))
  cat("Get_Davies_PVal(5, c(1,2,3)) =", result, "\n")
}, error=function(e) cat("Error:", conditionMessage(e), "\n"))

# Another test
tryCatch({
  result2 <- SKAT:::Get_Davies_PVal(1.5, c(0.5, 0.3, 0.1))
  cat("Get_Davies_PVal(1.5, c(0.5, 0.3, 0.1)) =", result2, "\n")
}, error=function(e) cat("Error:", conditionMessage(e), "\n"))
