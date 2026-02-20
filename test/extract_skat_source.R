#!/usr/bin/env Rscript
# Extract SKAT internal function sources for reference

library(SKAT)

outfile <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Step_2_Feb_11/test/skat_source_dump.txt"

sink(outfile)

cat("=== Get_Davies_PVal ===\n")
print(SKAT:::Get_Davies_PVal)

cat("\n\n=== Met_SKAT_Get_Pvalue ===\n")
print(SKAT:::Met_SKAT_Get_Pvalue)

cat("\n\n=== SKAT_Optimal_Each_Q ===\n")
print(SKAT:::SKAT_Optimal_Each_Q)

cat("\n\n=== SKAT_Optimal_PValue_Davies ===\n")
print(SKAT:::SKAT_Optimal_PValue_Davies)

cat("\n\n=== SKAT_Optimal_Param ===\n")
tryCatch(print(SKAT:::SKAT_Optimal_Param), error=function(e) cat("Not found\n"))

cat("\n\n=== SKAT_liu.MOD ===\n")
tryCatch(print(SKAT:::SKAT_liu.MOD), error=function(e) cat("Not found\n"))

cat("\n\n=== SKAT_liu.MOD.Lambda ===\n")
tryCatch(print(SKAT:::SKAT_liu.MOD.Lambda), error=function(e) cat("Not found\n"))

cat("\n\n=== Get_Lambda ===\n")
tryCatch(print(SKAT:::Get_Lambda), error=function(e) cat("Not found\n"))

cat("\n\n=== SKAT_GET_pvalue ===\n")
tryCatch(print(SKAT:::SKAT_GET_pvalue), error=function(e) cat("Not found\n"))

sink()

cat("Written to:", outfile, "\n")
