# Test harness for CCT (Cauchy Combination Test) - R side
# Implements the same CCT function as SAIGE to compare against C++ standalone

# CCT function matching SAIGE's implementation
CCT <- function(pvals, weights = rep(1/length(pvals), length(pvals))) {
    # Check for NAs
    if (any(is.na(pvals))) stop("Cannot have NAs in the p-values!")
    # Check range
    if (any(pvals < 0) || any(pvals > 1)) stop("P-values must be between 0 and 1")
    # If any zero
    if (any(pvals == 0)) return(0)
    # If any one
    if (any(pvals == 1)) {
        minp <- min(pvals) * length(pvals)
        return(min(1, minp))
    }

    # For very small p-values, use approximation: tan((0.5 - p) * pi) ~ 1/(p*pi)
    is.small <- (pvals < 1e-16)
    if (all(!is.small)) {
        cctstat <- sum(weights * tan((0.5 - pvals) * pi))
    } else {
        cctstat <- sum((weights[is.small] / pvals[is.small]) / pi)
        cctstat <- cctstat + sum(weights[!is.small] * tan((0.5 - pvals[!is.small]) * pi))
    }

    # Convert to p-value
    if (cctstat > 1e+15) {
        pval <- (1/cctstat) / pi
    } else {
        pval <- pcauchy(cctstat, lower.tail = FALSE)
    }
    return(pval)
}

options(digits = 15)

# Test case 1: {0.01, 0.05, 0.1, 0.5}
pvals1 <- c(0.01, 0.05, 0.1, 0.5)
cat(sprintf("CCT_test1: %.15f\n", CCT(pvals1)))

# Test case 2: {1e-5, 0.3, 0.7, 0.9}
pvals2 <- c(1e-5, 0.3, 0.7, 0.9)
cat(sprintf("CCT_test2: %.15f\n", CCT(pvals2)))

# Test case 3: {0.5, 0.5, 0.5, 0.5}
pvals3 <- c(0.5, 0.5, 0.5, 0.5)
cat(sprintf("CCT_test3: %.15f\n", CCT(pvals3)))

# Test case 4: {1e-10, 1e-8, 1e-6} (extreme small)
pvals4 <- c(1e-10, 1e-8, 1e-6)
cat(sprintf("CCT_test4: %.15f\n", CCT(pvals4)))

# Test case 5: {0.001} (single value)
pvals5 <- c(0.001)
cat(sprintf("CCT_test5: %.15f\n", CCT(pvals5)))
