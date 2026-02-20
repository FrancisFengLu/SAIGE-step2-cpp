#!/bin/bash
# Master test runner: compiles C++ tests, runs them, runs R tests, and compares outputs
# Usage: cd function_test && bash run_all.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

PIXI_CMD="$HOME/.pixi/bin/pixi run --manifest-path=/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/code_copy/SAIGE_isolated/pixi.toml Rscript"

echo "============================================"
echo "  SAIGE C++ vs R Function Comparison Tests"
echo "============================================"
echo ""

# ---- Step 1: Compile all C++ tests ----
echo "--- Compiling C++ tests ---"
make -f Makefile.tests clean 2>/dev/null || true
make -f Makefile.tests all
echo "Compilation successful."
echo ""

# ---- Step 2: Run C++ tests ----
echo "--- Running C++ test_cct ---"
./test_cct > cpp_cct_output.txt 2>&1
cat cpp_cct_output.txt
echo ""

echo "--- Running C++ test_util ---"
./test_util > cpp_util_output.txt 2>&1
cat cpp_util_output.txt
echo ""

echo "--- Running C++ test_spa ---"
./test_spa > cpp_spa_output.txt 2>&1
cat cpp_spa_output.txt
echo ""

# ---- Step 3: Run R tests ----
echo "--- Running R test_cct.R ---"
$PIXI_CMD "$SCRIPT_DIR/test_cct.R" > r_cct_output.txt 2>&1
cat r_cct_output.txt
echo ""

echo "--- Running R test_util.R ---"
$PIXI_CMD "$SCRIPT_DIR/test_util.R" > r_util_output.txt 2>&1
cat r_util_output.txt
echo ""

echo "--- Running R test_spa.R ---"
$PIXI_CMD "$SCRIPT_DIR/test_spa.R" > r_spa_output.txt 2>&1
cat r_spa_output.txt
echo ""

# ---- Step 4: Compare outputs ----
echo "============================================"
echo "  COMPARISON RESULTS"
echo "============================================"
echo ""

# Write comparison script
cat > compare_outputs.R << 'REOF'
# Compare C++ and R outputs line by line
# Reads pairs of output files and compares numeric values

compare_files <- function(cpp_file, r_file, test_name) {
    cpp_lines <- readLines(cpp_file)
    r_lines <- readLines(r_file)

    results <- data.frame(
        Test = character(),
        Label = character(),
        Cpp_Value = character(),
        R_Value = character(),
        Match = character(),
        Difference = character(),
        stringsAsFactors = FALSE
    )

    # Parse lines into label:value pairs
    parse_line <- function(line) {
        parts <- strsplit(line, ": ", fixed = TRUE)[[1]]
        if (length(parts) >= 2) {
            label <- parts[1]
            value <- paste(parts[-1], collapse = ": ")
            return(list(label = trimws(label), value = trimws(value)))
        }
        return(NULL)
    }

    cpp_parsed <- lapply(cpp_lines, parse_line)
    r_parsed <- lapply(r_lines, parse_line)

    # Remove NULLs and header lines
    cpp_parsed <- cpp_parsed[!sapply(cpp_parsed, is.null)]
    r_parsed <- r_parsed[!sapply(r_parsed, is.null)]

    # Match by label
    for (cp in cpp_parsed) {
        matching_r <- NULL
        for (rp in r_parsed) {
            if (cp$label == rp$label) {
                matching_r <- rp
                break
            }
        }

        if (!is.null(matching_r)) {
            cpp_num <- suppressWarnings(as.numeric(cp$value))
            r_num <- suppressWarnings(as.numeric(matching_r$value))

            if (!is.na(cpp_num) && !is.na(r_num)) {
                if (r_num == 0 && cpp_num == 0) {
                    diff_str <- "0"
                    match_str <- "EXACT"
                } else if (r_num == 0) {
                    diff_str <- sprintf("%.2e", abs(cpp_num))
                    match_str <- if (abs(cpp_num) < 1e-12) "YES" else "NO"
                } else {
                    rel_diff <- abs(cpp_num - r_num) / max(abs(r_num), 1e-300)
                    abs_diff <- abs(cpp_num - r_num)
                    diff_str <- sprintf("rel=%.2e abs=%.2e", rel_diff, abs_diff)
                    if (rel_diff < 1e-12) {
                        match_str <- "EXACT"
                    } else if (rel_diff < 1e-8) {
                        match_str <- "YES"
                    } else if (rel_diff < 1e-4) {
                        match_str <- "CLOSE"
                    } else {
                        match_str <- "NO"
                    }
                }
            } else {
                diff_str <- "N/A (non-numeric)"
                match_str <- if (cp$value == matching_r$value) "EXACT" else "MISMATCH"
            }

            results <- rbind(results, data.frame(
                Test = test_name,
                Label = cp$label,
                Cpp_Value = cp$value,
                R_Value = matching_r$value,
                Match = match_str,
                Difference = diff_str,
                stringsAsFactors = FALSE
            ))
        } else {
            # Skip header lines like "=== ... ==="
            if (!grepl("^===", cp$label)) {
                results <- rbind(results, data.frame(
                    Test = test_name,
                    Label = cp$label,
                    Cpp_Value = cp$value,
                    R_Value = "(no R match)",
                    Match = "N/A",
                    Difference = "N/A",
                    stringsAsFactors = FALSE
                ))
            }
        }
    }

    return(results)
}

all_results <- data.frame()
all_results <- rbind(all_results, compare_files("cpp_cct_output.txt", "r_cct_output.txt", "CCT"))
all_results <- rbind(all_results, compare_files("cpp_util_output.txt", "r_util_output.txt", "UTIL"))
all_results <- rbind(all_results, compare_files("cpp_spa_output.txt", "r_spa_output.txt", "SPA"))

# Print summary table
cat(sprintf("%-8s | %-35s | %-22s | %-22s | %-7s | %s\n",
            "Test", "Label", "C++ Value", "R Value", "Match?", "Difference"))
cat(paste(rep("-", 130), collapse = ""), "\n")

for (i in seq_len(nrow(all_results))) {
    r <- all_results[i, ]
    cat(sprintf("%-8s | %-35s | %-22s | %-22s | %-7s | %s\n",
                r$Test, r$Label,
                substr(r$Cpp_Value, 1, 22),
                substr(r$R_Value, 1, 22),
                r$Match, r$Difference))
}

# Write to RESULTS.txt
sink("RESULTS.txt")
cat("SAIGE C++ vs R Function Comparison Test Results\n")
cat(paste("Date:", Sys.time(), "\n"))
cat(paste(rep("=", 130), collapse = ""), "\n\n")

cat(sprintf("%-8s | %-35s | %-22s | %-22s | %-7s | %s\n",
            "Test", "Label", "C++ Value", "R Value", "Match?", "Difference"))
cat(paste(rep("-", 130), collapse = ""), "\n")

for (i in seq_len(nrow(all_results))) {
    r <- all_results[i, ]
    cat(sprintf("%-8s | %-35s | %-22s | %-22s | %-7s | %s\n",
                r$Test, r$Label,
                substr(r$Cpp_Value, 1, 22),
                substr(r$R_Value, 1, 22),
                r$Match, r$Difference))
}

cat("\n")
n_exact <- sum(all_results$Match == "EXACT")
n_yes <- sum(all_results$Match == "YES")
n_close <- sum(all_results$Match == "CLOSE")
n_no <- sum(all_results$Match == "NO")
n_na <- sum(all_results$Match %in% c("N/A", "MISMATCH"))
n_total <- nrow(all_results)

cat(sprintf("\nSummary: %d total comparisons\n", n_total))
cat(sprintf("  EXACT  (rel diff < 1e-12): %d\n", n_exact))
cat(sprintf("  YES    (rel diff < 1e-8):  %d\n", n_yes))
cat(sprintf("  CLOSE  (rel diff < 1e-4):  %d\n", n_close))
cat(sprintf("  NO     (rel diff >= 1e-4): %d\n", n_no))
cat(sprintf("  N/A    (non-numeric/missing): %d\n", n_na))

if (n_no == 0) {
    cat("\nOVERALL: PASS -- All numeric values match within tolerance.\n")
} else {
    cat("\nOVERALL: FAIL -- Some values differ beyond tolerance.\n")
}
sink()

cat("\nResults written to RESULTS.txt\n")

# Also print summary
n_exact <- sum(all_results$Match == "EXACT")
n_yes <- sum(all_results$Match == "YES")
n_close <- sum(all_results$Match == "CLOSE")
n_no <- sum(all_results$Match == "NO")
n_total <- nrow(all_results)

cat(sprintf("\nSummary: %d total comparisons\n", n_total))
cat(sprintf("  EXACT  (rel diff < 1e-12): %d\n", n_exact))
cat(sprintf("  YES    (rel diff < 1e-8):  %d\n", n_yes))
cat(sprintf("  CLOSE  (rel diff < 1e-4):  %d\n", n_close))
cat(sprintf("  NO     (rel diff >= 1e-4): %d\n", n_no))

if (n_no == 0) {
    cat("\nOVERALL: PASS\n")
} else {
    cat("\nOVERALL: FAIL\n")
}
REOF

echo "--- Comparing C++ vs R outputs ---"
$PIXI_CMD "$SCRIPT_DIR/compare_outputs.R"

echo ""
echo "============================================"
echo "  Done. See RESULTS.txt for full report."
echo "============================================"
