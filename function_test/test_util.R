# Test harness for UTIL::getWeights - R side
# Compares R's dbeta() against C++ standalone getWeights()

options(digits = 15)

mafs <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5)

# Beta(1, 25) -- SAIGE default
cat("=== getWeights with Beta(1, 25) ===\n")
for (maf in mafs) {
    w <- dbeta(maf, 1, 25)
    cat(sprintf("getWeights_beta1_25_maf%s: %.15f\n", maf, w))
}

# Beta(1, 1) -- uniform
cat("=== getWeights with Beta(1, 1) ===\n")
for (maf in mafs) {
    w <- dbeta(maf, 1, 1)
    cat(sprintf("getWeights_beta1_1_maf%s: %.15f\n", maf, w))
}

# Linear kernel (always 1)
cat("=== getWeights with linear kernel ===\n")
for (maf in mafs) {
    cat(sprintf("getWeights_linear_maf%s: %.15f\n", maf, 1.0))
}
