# Test harness for SPA Binary (Saddlepoint Approximation) - R side
# Implements the same SPA functions as SAIGE to compare against C++ standalone

options(digits = 15)

# ---- CGF and derivatives for Binomial SPA ----

Korg_Binom <- function(t1, mu, g) {
    temp <- log(1 - mu + mu * exp(g * t1))
    return(sum(temp))
}

K1_adj_Binom <- function(t1, mu, g, q) {
    temp1 <- (1 - mu) * exp(-g * t1) + mu
    temp2 <- mu * g
    temp3 <- temp2 / temp1
    return(sum(temp3) - q)
}

K2_Binom <- function(t1, mu, g) {
    temp0 <- exp(-g * t1)
    temp1 <- (1 - mu) * temp0 + mu
    temp1 <- temp1^2
    temp2 <- g^2 * temp0
    temp2 <- (1 - mu) * mu * temp2
    temp3 <- temp2 / temp1
    # sum_arma1: sum only finite values
    return(sum(temp3[is.finite(temp3)]))
}

getroot_K1_Binom <- function(init, mu, g, q, tol, maxiter = 1000) {
    gpos <- sum(g[g > 0])
    gneg <- sum(g[g < 0])

    if (q >= gpos || q <= gneg) {
        return(list(root = Inf, niter = 0, Isconverge = TRUE))
    }

    t <- init
    K1_eval <- K1_adj_Binom(t, mu, g, q)
    prevJump <- Inf
    conv <- TRUE
    rep <- 1

    while (rep <= maxiter) {
        K2_eval <- K2_Binom(t, mu, g)
        tnew <- t - K1_eval / K2_eval

        if (is.nan(tnew)) {
            conv <- FALSE
            break
        }

        if (abs(tnew - t) < tol) {
            conv <- TRUE
            break
        }

        if (rep == maxiter) {
            conv <- FALSE
            break
        }

        newK1 <- K1_adj_Binom(tnew, mu, g, q)
        if (sign(K1_eval) != sign(newK1)) {
            if (abs(tnew - t) > (prevJump - tol)) {
                tnew <- t + sign(newK1 - K1_eval) * prevJump / 2
                newK1 <- K1_adj_Binom(tnew, mu, g, q)
                prevJump <- prevJump / 2
            } else {
                prevJump <- abs(tnew - t)
            }
        }

        rep <- rep + 1
        t <- tnew
        K1_eval <- newK1
    }

    return(list(root = t, niter = rep, Isconverge = conv))
}

Get_Saddle_Prob_Binom <- function(zeta, mu, g, q, logp = FALSE) {
    k1 <- Korg_Binom(zeta, mu, g)
    k2 <- K2_Binom(zeta, mu, g)
    temp1 <- zeta * q - k1
    isSaddle <- FALSE

    flagrun <- FALSE
    if (is.finite(k1) && is.finite(k2) && temp1 >= 0 && k2 >= 0) {
        w <- sign(zeta) * sqrt(2 * temp1)
        v <- zeta * sqrt(k2)
        if (w != 0) {
            flagrun <- TRUE
        }
    }

    if (flagrun) {
        Ztest <- w + (1/w) * log(v/w)

        if (Ztest > 0) {
            if (logp) {
                pval0 <- pnorm(Ztest, 0, 1, lower.tail = FALSE, log.p = TRUE)
            } else {
                pval0 <- pnorm(Ztest, 0, 1, lower.tail = FALSE)
            }
            pval <- pval0
        } else {
            if (logp) {
                pval0 <- pnorm(Ztest, 0, 1, lower.tail = TRUE, log.p = TRUE)
            } else {
                pval0 <- pnorm(Ztest, 0, 1, lower.tail = TRUE)
            }
            pval <- -pval0
        }
        isSaddle <- TRUE
    } else {
        if (logp) {
            pval <- -Inf
        } else {
            pval <- 0
        }
    }

    return(list(pval = pval, isSaddle = isSaddle))
}

add_logp <- function(p1, p2) {
    p1 <- -abs(p1)
    p2 <- -abs(p2)
    maxp <- max(p1, p2)
    minp <- min(p1, p2)
    result <- maxp + log(1 + exp(minp - maxp))
    return(result)
}

SPA_binary_R <- function(mu, g, q, qinv, pval_noadj, tol, logp = FALSE) {
    outuni1 <- getroot_K1_Binom(0, mu, g, q, tol)
    outuni2 <- getroot_K1_Binom(0, mu, g, qinv, tol)
    Isconverge <- TRUE

    if (outuni1$Isconverge && outuni2$Isconverge) {
        getSaddle <- Get_Saddle_Prob_Binom(outuni1$root, mu, g, q, logp)
        if (getSaddle$isSaddle) {
            p1 <- getSaddle$pval
        } else {
            if (logp) { p1 <- pval_noadj - log(2) } else { p1 <- pval_noadj / 2 }
        }

        getSaddle2 <- Get_Saddle_Prob_Binom(outuni2$root, mu, g, qinv, logp)
        if (getSaddle2$isSaddle) {
            p2 <- getSaddle2$pval
        } else {
            if (logp) { p2 <- pval_noadj - log(2) } else { p2 <- pval_noadj / 2 }
        }

        if (logp) {
            pval <- add_logp(p1, p2)
        } else {
            pval <- abs(p1) + abs(p2)
        }
        Isconverge <- TRUE
    } else {
        pval <- pval_noadj
        Isconverge <- FALSE
    }

    return(list(pval = pval, Isconverge = Isconverge))
}

# The SPA() dispatcher in R matches SPA_binary for traitType="binary"
SPA_R <- function(mu, g, q, qinv, pval_noadj, tol, logp, traitType) {
    if (traitType == "binary") {
        outuni1 <- getroot_K1_Binom(0, mu, g, q, tol)
        outuni2 <- getroot_K1_Binom(0, mu, g, qinv, tol)
    }
    Isconverge <- TRUE

    if (outuni1$Isconverge && outuni2$Isconverge) {
        getSaddle <- Get_Saddle_Prob_Binom(outuni1$root, mu, g, q, logp)
        getSaddle2 <- Get_Saddle_Prob_Binom(outuni2$root, mu, g, qinv, logp)

        if (getSaddle$isSaddle) {
            p1 <- getSaddle$pval
        } else {
            Isconverge <- FALSE
            if (logp) { p1 <- pval_noadj - log(2) } else { p1 <- pval_noadj / 2 }
        }
        if (getSaddle2$isSaddle) {
            p2 <- getSaddle2$pval
        } else {
            Isconverge <- FALSE
            if (logp) { p2 <- pval_noadj - log(2) } else { p2 <- pval_noadj / 2 }
        }

        if (logp) {
            pval <- add_logp(p1, p2)
        } else {
            pval <- abs(p1) + abs(p2)
        }
    } else {
        pval <- pval_noadj
        Isconverge <- FALSE
    }
    return(list(pval = pval, Isconverge = Isconverge))
}

# ---- Run tests ----

mu <- c(0.3, 0.7, 0.2, 0.8, 0.5, 0.1, 0.9, 0.4, 0.6, 0.15)
g  <- c(1.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 2.0)

q_mean <- sum(mu * g)
q <- q_mean + 1.5
qinv <- 2 * q_mean - q

cat(sprintf("mu*g_mean: %.15f\n", q_mean))
cat(sprintf("q: %.15f\n", q))
cat(sprintf("qinv: %.15f\n", qinv))

# Test 1: Korg_Binom at t=0.1
cat(sprintf("Korg_Binom_t0.1: %.15f\n", Korg_Binom(0.1, mu, g)))

# Test 2: K1_adj_Binom at t=0.1
cat(sprintf("K1_adj_Binom_t0.1: %.15f\n", K1_adj_Binom(0.1, mu, g, q)))

# Test 3: K2_Binom at t=0.1
cat(sprintf("K2_Binom_t0.1: %.15f\n", K2_Binom(0.1, mu, g)))

# Test 4: Korg_Binom at t=0
cat(sprintf("Korg_Binom_t0: %.15f\n", Korg_Binom(0, mu, g)))

# Test 5: K1_adj_Binom at t=0
cat(sprintf("K1_adj_Binom_t0: %.15f\n", K1_adj_Binom(0, mu, g, q)))

# Test 6: K2_Binom at t=0
cat(sprintf("K2_Binom_t0: %.15f\n", K2_Binom(0, mu, g)))

# Test 7: getroot_K1_Binom for q
tol <- 0.0001220703125  # (0.5)^20, SAIGE default
rr <- getroot_K1_Binom(0, mu, g, q, tol)
cat(sprintf("getroot_root: %.15f\n", rr$root))
cat(sprintf("getroot_niter: %d\n", rr$niter))
cat(sprintf("getroot_converge: %d\n", as.integer(rr$Isconverge)))

# Test 8: getroot for qinv
rr2 <- getroot_K1_Binom(0, mu, g, qinv, tol)
cat(sprintf("getroot_qinv_root: %.15f\n", rr2$root))
cat(sprintf("getroot_qinv_niter: %d\n", rr2$niter))
cat(sprintf("getroot_qinv_converge: %d\n", as.integer(rr2$Isconverge)))

# Test 9: Get_Saddle_Prob_Binom
if (rr$Isconverge && is.finite(rr$root)) {
    sr <- Get_Saddle_Prob_Binom(rr$root, mu, g, q, logp = FALSE)
    cat(sprintf("SaddleProb_pval: %.15f\n", sr$pval))
    cat(sprintf("SaddleProb_isSaddle: %d\n", as.integer(sr$isSaddle)))
} else {
    cat("SaddleProb_pval: inf (root not finite)\n")
    cat("SaddleProb_isSaddle: 0\n")
}

# Test 10: SPA_binary full pipeline
pval_noadj <- 0.05
sr_binary <- SPA_binary_R(mu, g, q, qinv, pval_noadj, tol, logp = FALSE)
cat(sprintf("SPA_binary_pval: %.15f\n", sr_binary$pval))
cat(sprintf("SPA_binary_converge: %d\n", as.integer(sr_binary$Isconverge)))

# Test 11: SPA() dispatcher (binary trait)
sr_disp <- SPA_R(mu, g, q, qinv, pval_noadj, tol, logp = FALSE, traitType = "binary")
cat(sprintf("SPA_dispatcher_pval: %.15f\n", sr_disp$pval))
cat(sprintf("SPA_dispatcher_converge: %d\n", as.integer(sr_disp$Isconverge)))

# Test 12: SPA_pval function (same as dispatcher for binary)
sr_pval <- SPA_R(mu, g, q, qinv, pval_noadj, tol, logp = FALSE, traitType = "binary")
cat(sprintf("SPA_pval_func: %.15f\n", sr_pval$pval))
cat(sprintf("SPA_pval_converge: %d\n", as.integer(sr_pval$Isconverge)))
