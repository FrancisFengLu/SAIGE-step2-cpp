#!/usr/bin/env Rscript
# Extract all SKAT internal function sources for porting to C++

library(SKAT)

outfile <- "skat_internals.txt"
sink(outfile)

cat("=== SKAT_Optimal_PValue_Davies ===\n")
tryCatch(print(SKAT:::SKAT_Optimal_PValue_Davies), error=function(e) cat("Not found:", conditionMessage(e), "\n"))

cat("\n\n=== SKAT_Optimal_Integrate_Func_Davies ===\n")
tryCatch(print(SKAT:::SKAT_Optimal_Integrate_Func_Davies), error=function(e) cat("Not found:", conditionMessage(e), "\n"))

cat("\n\n=== SKAT_Optimal_Param ===\n")
tryCatch(print(SKAT:::SKAT_Optimal_Param), error=function(e) cat("Not found:", conditionMessage(e), "\n"))

cat("\n\n=== Get_Davies_PVal ===\n")
tryCatch(print(SKAT:::Get_Davies_PVal), error=function(e) cat("Not found:", conditionMessage(e), "\n"))

cat("\n\n=== Get_Lambda ===\n")
tryCatch(print(SKAT:::Get_Lambda), error=function(e) cat("Not found:", conditionMessage(e), "\n"))

cat("\n\n=== SKAT_liu.MOD.Lambda ===\n")
tryCatch(print(SKAT:::SKAT_liu.MOD.Lambda), error=function(e) cat("Not found:", conditionMessage(e), "\n"))

cat("\n\n=== Met_SKAT_Get_Pvalue ===\n")
tryCatch(print(SKAT:::Met_SKAT_Get_Pvalue), error=function(e) cat("Not found:", conditionMessage(e), "\n"))

cat("\n\n=== SKAT_Optimal_Each_Q ===\n")
tryCatch(print(SKAT:::SKAT_Optimal_Each_Q), error=function(e) cat("Not found:", conditionMessage(e), "\n"))

cat("\n\n=== SKAT_Optimal_PValue ===\n")
tryCatch(print(SKAT:::SKAT_Optimal_PValue), error=function(e) cat("Not found:", conditionMessage(e), "\n"))

sink()
cat("Wrote SKAT internal function sources to", outfile, "\n")

# Now run a simple test to get intermediate values
cat("\n=== Running SKAT test cases ===\n")

# Test case 1: simple 2-eigenvalue case
lambda1 <- c(3.5, 1.2)
Q1 <- 5.0
cat("Test 1: lambda=c(3.5, 1.2), Q=5.0\n")
cat("  Davies:", tryCatch(SKAT:::Get_Davies_PVal(Q1, lambda1), error=function(e) NA), "\n")
cat("  Liu:", tryCatch({
    lp <- SKAT:::SKAT_liu.MOD.Lambda(lambda1)
    cat(sprintf("    muQ=%g, sigmaQ=%g, l=%g, delta=%g\n", lp$muQ, lp$sigmaQ, lp$l, lp$delta))
    Qnorm <- (Q1 - lp$muQ) / lp$sigmaQ * sqrt(2*lp$l + 4*lp$delta) + lp$l + lp$delta
    cat(sprintf("    Qnorm=%g\n", Qnorm))
    pchisq(Qnorm, lp$l, ncp=lp$delta, lower.tail=FALSE)
}, error=function(e) {cat("Error:", conditionMessage(e), "\n"); NA}), "\n")

# Test case 2: 4-eigenvalue case (similar to GENE1 missense;lof 0.01)
lambda2 <- c(10000, 7000, 5000, 3000)
Q2 <- 20000
cat("\nTest 2: lambda=c(10000, 7000, 5000, 3000), Q=20000\n")
cat("  Davies:", tryCatch(SKAT:::Get_Davies_PVal(Q2, lambda2), error=function(e) NA), "\n")
cat("  Liu:", tryCatch({
    lp <- SKAT:::SKAT_liu.MOD.Lambda(lambda2)
    cat(sprintf("    muQ=%g, sigmaQ=%g, l=%g, delta=%g\n", lp$muQ, lp$sigmaQ, lp$l, lp$delta))
    Qnorm <- (Q2 - lp$muQ) / lp$sigmaQ * sqrt(2*lp$l + 4*lp$delta) + lp$l + lp$delta
    cat(sprintf("    Qnorm=%g\n", Qnorm))
    pchisq(Qnorm, lp$l, ncp=lp$delta, lower.tail=FALSE)
}, error=function(e) {cat("Error:", conditionMessage(e), "\n"); NA}), "\n")
