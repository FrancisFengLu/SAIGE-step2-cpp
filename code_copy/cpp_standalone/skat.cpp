// Standalone port of SKAT/Burden/SKAT-O tests for SAIGE Step 2
//
// Implements:
//   1. Davies method (characteristic function inversion) for mixture of chi-squared
//   2. Liu moment-matching method (fallback when Davies fails)
//   3. SKAT-O optimal.adj method
//   4. get_SKAT_pvalue() matching SAIGE's R wrapper
//
// Original R sources:
//   - SAIGE/R/SAIGE_SPATest_Region_Func.R :: get_SKAT_pvalue()
//   - SKAT:::Met_SKAT_Get_Pvalue (R/SKAT_Optimal_Adj.R)
//   - CompQuadForm::davies (qfc.c)
//   - SKAT:::SKAT_liu.MOD.Lambda
//
// Conversions from R/Rcpp:
//   - R::pchisq  -> boost::math::cdf(chi_squared)
//   - R::qnorm   -> boost::math::quantile(normal)
//   - R::qchisq  -> boost::math::quantile(chi_squared)
//   - Rcpp::List  -> C++ structs
//   - Rcpp::stop  -> throw std::runtime_error

#include "skat.hpp"

#include <armadillo>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <limits>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>


// ============================================================
// Davies method: P(Q > q) where Q = sum(lambda_j * chi2(1))
//
// Uses characteristic function inversion via the trapezoidal rule.
// The characteristic function of Q = sum(lambda_j * X_j^2) where X_j ~ N(0,1) is:
//   phi(t) = prod_j (1 - 2*i*t*lambda_j)^(-1/2)
//
// The survival function is:
//   P(Q > q) = 0.5 - (1/pi) * integral_0^inf Im[phi(t)*exp(-i*t*q)] / t dt
//
// The integrand in polar form:
//   f(t) = exp(log_amplitude) * sin(phase) / t
// where:
//   log_amplitude = sum_j (-1/4) * log(1 + 4*t^2*lambda_j^2)
//   phase = sum_j (1/2) * atan(2*t*lambda_j) - t*q
//
// Extended for general chi-squared(n_j, delta_j) components:
//   log_amplitude += -delta_j * (2*t*lambda_j)^2 / (2*(1 + (2*t*lambda_j)^2))
//   phase += delta_j * (2*t*lambda_j) / (2*(1 + (2*t*lambda_j)^2))
//
// This implementation handles the general case:
//   Q = sum_j lambda_j * chi2(n_j, delta_j) + sigma*N(0,1)
// but SAIGE only needs: n_j=1, delta_j=0, sigma=0.
// ============================================================

namespace davies_impl {

struct DaviesParams {
    std::vector<double> lambda;   // eigenvalues
    std::vector<int> mult;        // multiplicities (degrees of freedom per term)
    std::vector<double> delta;    // noncentrality parameters
    double sigma;                 // std dev of additional normal component
    double c;                     // threshold: compute P(Q > c)
    int lim;                      // max integration terms
    double acc;                   // target accuracy
};

// Compute the integrand: Im[phi(t)*exp(-i*t*c)] / t
static double compute_integrand(double t, const DaviesParams& p) {
    double log_amp = 0.0;
    double phase = -t * p.c;
    int n = (int)p.lambda.size();

    for (int j = 0; j < n; j++) {
        double x = 2.0 * t * p.lambda[j];
        double x2 = x * x;
        double denom = 1.0 + x2;

        // Amplitude: (1 + x^2)^(-mult_j/4)
        log_amp -= 0.25 * p.mult[j] * std::log(denom);

        // Phase: (mult_j/2) * atan(x)
        phase += 0.5 * p.mult[j] * std::atan(x);

        // Noncentrality contributions
        if (p.delta[j] > 0.0) {
            log_amp -= 0.5 * p.delta[j] * x2 / denom;
            phase += 0.5 * p.delta[j] * x / denom;
        }
    }

    // Additional normal component
    if (p.sigma > 0.0) {
        log_amp -= 0.5 * p.sigma * p.sigma * t * t;
    }

    // Protect against very small amplitudes
    if (log_amp < -50.0) return 0.0;

    return std::exp(log_amp) * std::sin(phase) / t;
}

// Main Davies computation
// Returns {P(Q > c), ifault}
// ifault: 0=success, 1=not converged, 2=round-off error
static std::pair<double, int> compute(const DaviesParams& p) {
    int n = (int)p.lambda.size();
    if (n == 0) {
        double pval = (p.c > 0.0) ? 0.0 : 1.0;
        return {pval, 0};
    }

    // Find max |eigenvalue| for step size selection
    double lmax = 0.0;
    for (int j = 0; j < n; j++) {
        double al = std::abs(p.lambda[j]);
        if (al > lmax) lmax = al;
    }

    if (lmax == 0.0) {
        return {(p.c > 0.0) ? 0.0 : 1.0, 0};
    }

    // Compute mean of Q for determining integration range
    double mean_q = 0.0;
    double var_q = p.sigma * p.sigma;
    for (int j = 0; j < n; j++) {
        mean_q += p.lambda[j] * (p.mult[j] + p.delta[j]);
        var_q += 2.0 * p.lambda[j] * p.lambda[j] * (p.mult[j] + 2.0 * p.delta[j]);
    }
    double sd_q = std::sqrt(var_q);

    // Step size for trapezoidal rule
    // Need: 2 * lmax * h < pi to avoid aliasing
    // Also want h small enough to capture the integrand shape
    // A good heuristic: h = pi / (4 * max(lmax, |c-mean|/sd))
    double almx = lmax;
    if (sd_q > 0.0) {
        almx = std::max(lmax, std::abs(p.c - mean_q) / sd_q);
    }
    double h = 1.0 / (2.0 * almx);
    h = std::min(h, M_PI / (4.0 * lmax));  // Nyquist condition

    // Trapezoidal integration
    double integral = 0.0;
    int count = 0;
    int consecutive_small = 0;
    int ifault = 0;

    for (int k = 1; k <= p.lim; k++) {
        double t = k * h;
        double val = compute_integrand(t, p);
        integral += val;
        count++;

        // Convergence check: several consecutive small contributions
        double contrib = std::abs(val * h);
        if (k > 20 && contrib < p.acc * 1e-4) {
            consecutive_small++;
            if (consecutive_small > 10) break;
        } else {
            consecutive_small = 0;
        }

        // Also check if amplitude has decayed substantially
        if (k > 200 && contrib < p.acc * 1e-8) break;
    }

    if (count >= p.lim) ifault = 1;

    double pvalue = 0.5 - h * integral / M_PI;

    // Clamp to [0, 1]
    if (pvalue < 0.0) {
        if (pvalue > -4.0 * p.acc) {
            pvalue = 0.0;
        } else {
            ifault = 2;
            pvalue = 0.0;
        }
    }
    if (pvalue > 1.0) pvalue = 1.0;

    return {pvalue, ifault};
}

}  // namespace davies_impl


// ============================================================
// Public Davies interface
// ============================================================

double davies_pvalue(double q, const arma::vec& lambda) {
    // Filter out zero/near-zero eigenvalues
    arma::vec lam = lambda(arma::find(arma::abs(lambda) > 1e-10));
    int n = (int)lam.n_elem;

    if (n == 0) {
        return (q > 0.0) ? 0.0 : 1.0;
    }

    // Special case: single eigenvalue -> simple scaled chi-squared(1)
    if (n == 1) {
        double l1 = lam(0);
        if (l1 <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        double q_scaled = q / l1;
        if (q_scaled <= 0.0) return 1.0;
        try {
            boost::math::chi_squared chi2_1(1.0);
            return boost::math::cdf(boost::math::complement(chi2_1, q_scaled));
        } catch (...) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    // Build parameters
    davies_impl::DaviesParams params;
    params.lambda.resize(n);
    params.mult.resize(n, 1);
    params.delta.resize(n, 0.0);
    params.sigma = 0.0;
    params.c = q;
    params.lim = 10000;
    params.acc = 1e-6;

    for (int i = 0; i < n; i++) {
        params.lambda[i] = lam(i);
    }

    auto [pvalue, ifault] = davies_impl::compute(params);

    // If Davies failed, fall back to Liu's method
    if (ifault != 0 || pvalue < 0.0 || !std::isfinite(pvalue)) {
        return liu_pvalue(q, lambda);
    }

    return pvalue;
}


// ============================================================
// Liu moment-matching method
// ============================================================
// Approximates P(Q > q) where Q = sum(lambda_i * chi2(1))
// by matching the first 4 cumulants to a scaled noncentral chi-squared.
//
// From SKAT:::SKAT_liu.MOD.Lambda:
//   c_k = sum(lambda^k) for k=1,2,3,4
//   s1 = c3/c2^1.5, s2 = c4/c2^2
//   if s1^2 > s2: use noncentral chi-squared
//   else: use central chi-squared
//
// The normalized statistic is mapped to the chi-squared distribution.

double liu_pvalue(double q, const arma::vec& lambda) {
    arma::vec lam = lambda(arma::find(arma::abs(lambda) > 1e-10));
    int n = (int)lam.n_elem;

    if (n == 0) return (q > 0.0) ? 0.0 : 1.0;

    // Cumulant sums: c_k = sum(lambda^k)
    double c1 = arma::sum(lam);
    double c2 = arma::sum(lam % lam);
    double c3 = arma::sum(lam % lam % lam);
    double c4 = arma::sum(lam % lam % lam % lam);

    if (c2 <= 0.0) return (q > c1) ? 0.0 : 1.0;

    // Skewness and kurtosis ratios
    double s1 = c3 / std::pow(c2, 1.5);
    double s2 = c4 / (c2 * c2);

    double muQ = c1;
    double sigmaQ = std::sqrt(2.0 * c2);

    // Fit chi-squared approximation parameters
    double l, delta;  // df and noncentrality

    if (s1 * s1 > s2) {
        // Noncentral chi-squared approximation
        double a = 1.0 / (s1 - std::sqrt(s1 * s1 - s2));
        delta = s1 * a * a * a - a * a;
        l = a * a - 2.0 * delta;
    } else {
        // Central chi-squared approximation
        delta = 0.0;
        l = 1.0 / (s1 * s1);
    }

    // Ensure valid parameters
    if (l <= 0.0) l = 1.0;
    if (delta < 0.0) delta = 0.0;

    // Normalize Q to the chi-squared scale
    double Qnorm = (q - muQ) / sigmaQ * std::sqrt(2.0 * l + 4.0 * delta) + l + delta;

    if (Qnorm <= 0.0) return 1.0;

    // Compute p-value
    try {
        if (delta > 1e-6) {
            boost::math::non_central_chi_squared ncchi2(l, delta);
            return boost::math::cdf(boost::math::complement(ncchi2, Qnorm));
        } else {
            boost::math::chi_squared chi2(l);
            return boost::math::cdf(boost::math::complement(chi2, Qnorm));
        }
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}


// ============================================================
// Internal: Liu method returning all parameters (for SKAT-O)
// ============================================================

struct LiuParams {
    double pvalue;
    double muQ;
    double sigmaQ;
    double l;       // df
    double delta;   // noncentrality
};

static LiuParams liu_params(const arma::vec& lambda) {
    LiuParams res;
    res.pvalue = 1.0;
    res.muQ = 0.0;
    res.sigmaQ = 1.0;
    res.l = 1.0;
    res.delta = 0.0;

    arma::vec lam = lambda(arma::find(arma::abs(lambda) > 1e-10));
    if (lam.n_elem == 0) return res;

    double c1 = arma::sum(lam);
    double c2 = arma::sum(lam % lam);
    double c3 = arma::sum(lam % lam % lam);
    double c4 = arma::sum(lam % lam % lam % lam);

    if (c2 <= 0.0) return res;

    double s1 = c3 / std::pow(c2, 1.5);
    double s2 = c4 / (c2 * c2);

    res.muQ = c1;
    res.sigmaQ = std::sqrt(2.0 * c2);

    if (s1 * s1 > s2) {
        double a = 1.0 / (s1 - std::sqrt(s1 * s1 - s2));
        res.delta = s1 * a * a * a - a * a;
        res.l = a * a - 2.0 * res.delta;
    } else {
        res.delta = 0.0;
        res.l = 1.0 / (s1 * s1);
    }

    if (res.l <= 0.0) res.l = 1.0;
    if (res.delta < 0.0) res.delta = 0.0;

    return res;
}


// ============================================================
// Helper: compute eigenvalues of Phi_rho for a given rho
// ============================================================
// Q(rho) = Score' * R_rho * Score where R_rho = (1-rho)*I + rho*11'
// Under H0, Q(rho) ~ sum(lambda_j * chi2(1)) where lambda_j = eigenvalues
// of R_rho^{1/2} * Phi * R_rho^{1/2}

static arma::vec compute_phi_rho_eigenvalues(const arma::mat& Phi, double rho, int m) {
    if (std::abs(rho) < 1e-10) {
        // SKAT: eigenvalues of Phi
        arma::vec ev;
        arma::eig_sym(ev, Phi);
        return ev(arma::find(ev > 1e-10));
    }

    if (std::abs(rho - 1.0) < 1e-10) {
        // Burden: single eigenvalue = sum(Phi)
        double total = arma::accu(Phi);
        if (total > 1e-10) {
            arma::vec result(1);
            result(0) = total;
            return result;
        }
        return arma::vec();
    }

    // General case: R_rho = (1-rho)*I + rho*11'
    // Eigenvalues of R_rho: (1-rho) [multiplicity m-1] and (1-rho+m*rho) [multiplicity 1]
    // sqrt(R_rho) has same eigenvectors with sqrt eigenvalues
    double ev_small = 1.0 - rho;
    double ev_large = 1.0 - rho + m * rho;

    arma::vec ones_v = arma::ones<arma::vec>(m);
    arma::mat P_proj = ones_v * ones_v.t() / (double)m;  // projection onto 1
    arma::mat sqrtR = std::sqrt(ev_small) * (arma::eye<arma::mat>(m, m) - P_proj)
                    + std::sqrt(ev_large) * P_proj;

    arma::mat Phi_rho = sqrtR * Phi * sqrtR;
    Phi_rho = 0.5 * (Phi_rho + Phi_rho.t());  // ensure numerical symmetry

    arma::vec ev;
    arma::eig_sym(ev, Phi_rho);
    return ev(arma::find(ev > 1e-10));
}


// ============================================================
// SKAT-O: optimal.adj p-value
// ============================================================
// Implements SKAT:::Met_SKAT_Get_Pvalue with method="optimal.adj"
//
// Algorithm:
// 1. For each rho in grid, compute p-value of Q(rho) using Davies/Liu
// 2. Find T_min = min(p-values across rho)
// 3. Compute P(min_rho p(rho) < T_min) using numerical integration
//
// The integration factorizes Q(rho) into a leading eigenvalue component
// (integrated over chi2(1)) and a remainder (evaluated via Liu's method).

double SKATO_optimal_pvalue(const arma::vec& Score,
                             const arma::mat& Phi,
                             const arma::vec& rho_vec,
                             const arma::vec& pval_each,
                             const arma::vec& q_each) {
    int n_rho = (int)rho_vec.n_elem;
    int m = (int)Score.n_elem;

    // Find minimum p-value across the rho grid
    double min_pval = arma::min(pval_each);
    if (!std::isfinite(min_pval) || min_pval <= 0.0) {
        return min_pval;
    }

    // For each rho, get eigenvalues and Liu parameters
    std::vector<arma::vec> all_eigvals(n_rho);
    arma::vec lam_max(n_rho, arma::fill::zeros);
    std::vector<LiuParams> liu_vec(n_rho);
    arma::vec tau_rho(n_rho, arma::fill::zeros);

    for (int k = 0; k < n_rho; k++) {
        all_eigvals[k] = compute_phi_rho_eigenvalues(Phi, rho_vec(k), m);

        if (all_eigvals[k].n_elem == 0) {
            liu_vec[k] = LiuParams{1.0, 0.0, 1.0, 1.0, 0.0};
            tau_rho(k) = 0.0;
            continue;
        }

        lam_max(k) = arma::max(all_eigvals[k]);
        liu_vec[k] = liu_params(all_eigvals[k]);

        // Find tau_rho(k): the Q(rho) threshold such that
        // P(Q_approx > tau) = min_pval under the Liu approximation
        // tau_norm = qchisq(1 - min_pval, df=l, ncp=delta)
        // tau = (tau_norm - l - delta) / sqrt(2*l + 4*delta) * sigmaQ + muQ
        try {
            double q_norm;
            if (liu_vec[k].delta > 1e-6) {
                boost::math::non_central_chi_squared ncchi2(liu_vec[k].l, liu_vec[k].delta);
                q_norm = boost::math::quantile(boost::math::complement(ncchi2, min_pval));
            } else {
                boost::math::chi_squared chi2(liu_vec[k].l);
                q_norm = boost::math::quantile(boost::math::complement(chi2, min_pval));
            }
            double denom = std::sqrt(2.0 * liu_vec[k].l + 4.0 * liu_vec[k].delta);
            if (denom > 0.0) {
                tau_rho(k) = (q_norm - liu_vec[k].l - liu_vec[k].delta) / denom
                           * liu_vec[k].sigmaQ + liu_vec[k].muQ;
            } else {
                tau_rho(k) = q_each(k);
            }
        } catch (...) {
            tau_rho(k) = q_each(k);  // fallback
        }
    }

    // Numerical integration over chi2(1) to compute P(SKAT-O)
    //
    // P(min_rho p(rho) < min_pval)
    //   = P(exists k: Q(rho_k) > tau_k)
    //   = integral_0^inf max_k P(Q_remain_k > tau_k - lam_max_k * z) * f_chi2(z) dz
    //
    // where Q(rho_k) is decomposed as lam_max_k * Z + Q_remain_k,
    // Z ~ chi2(1) is the component from the largest eigenvalue,
    // and Q_remain_k = sum_{j>1} lam_j(rho_k) * chi2_1_j is independent of Z.

    boost::math::chi_squared chi2_1(1.0);
    double dx = 0.05;
    double x_max = 40.0;
    int n_grid = (int)(x_max / dx);
    double p_skato = 0.0;

    for (int i = 0; i < n_grid; i++) {
        double x_lo = i * dx;
        double x_hi = (i + 1) * dx;
        double x_mid = (x_lo + x_hi) / 2.0;

        // Weight from chi2(1) CDF
        double weight;
        try {
            weight = boost::math::cdf(chi2_1, x_hi) - boost::math::cdf(chi2_1, x_lo);
        } catch (...) {
            continue;
        }
        if (weight < 1e-20) continue;

        // For each rho, compute P(Q_remain > tau - lam_max * x_mid)
        double max_p_exceed = 0.0;
        for (int k = 0; k < n_rho; k++) {
            if (all_eigvals[k].n_elem == 0) continue;

            double remaining_threshold = tau_rho(k) - lam_max(k) * x_mid;

            if (remaining_threshold <= 0.0) {
                // Leading eigenvalue component alone exceeds threshold
                max_p_exceed = 1.0;
                break;
            }

            if (all_eigvals[k].n_elem <= 1) {
                // Only one eigenvalue; after removing it, no remaining terms
                // P(0 > positive_threshold) = 0
                continue;
            }

            // Get remaining eigenvalues (all except the largest)
            arma::vec sorted_ev = arma::sort(all_eigvals[k], "descend");
            arma::vec lam_remain = sorted_ev.subvec(1, sorted_ev.n_elem - 1);

            // P(Q_remain > remaining_threshold) via Liu
            double p_exceed = liu_pvalue(remaining_threshold, lam_remain);
            if (std::isnan(p_exceed) || p_exceed < 0.0) p_exceed = 0.0;

            max_p_exceed = std::max(max_p_exceed, p_exceed);
            if (max_p_exceed >= 1.0) break;
        }

        p_skato += max_p_exceed * weight;
    }

    // Clamp to valid range
    if (p_skato < 0.0) p_skato = 0.0;
    if (p_skato > 1.0) p_skato = 1.0;

    // SKAT-O p-value must be at least as large as the minimum per-rho p-value
    // (the correction can only make it larger)
    p_skato = std::max(p_skato, min_pval);

    // But should not exceed the Bonferroni bound
    double bonf = std::min(min_pval * (double)n_rho, 1.0);
    p_skato = std::min(p_skato, bonf);

    return p_skato;
}


// ============================================================
// get_SKAT_pvalue: Main entry point
// ============================================================
// Direct port of SAIGE's get_SKAT_pvalue R function.
// Calls SKAT/Burden/SKAT-O tests depending on regionTestType.
//
// Arguments:
//   Score: weighted score vector (m x 1), weights already applied
//   Phi: weighted variance-covariance matrix (m x m), weights already applied
//   r_corr: correlation parameter grid for SKAT-O
//   regionTestType: "SKAT", "BURDEN", or "SKATO"

SKATResult get_SKAT_pvalue(const arma::vec& Score,
                            const arma::mat& Phi,
                            const arma::vec& r_corr,
                            const std::string& regionTestType) {
    SKATResult result;
    result.pvalue_SKATO = std::numeric_limits<double>::quiet_NaN();
    result.pvalue_Burden = std::numeric_limits<double>::quiet_NaN();
    result.pvalue_SKAT = std::numeric_limits<double>::quiet_NaN();
    result.beta_Burden = std::numeric_limits<double>::quiet_NaN();
    result.se_Burden = std::numeric_limits<double>::quiet_NaN();

    int m = (int)Score.n_elem;
    if (m == 0) return result;

    // ----------------------------------------------------------------
    // Burden effect size: BETA_Burden = sum(Score) / sum(diag(Phi))
    // This matches: BETA_Burden = sum(Score)/(sum(diag(Phi))) in R
    // ----------------------------------------------------------------
    double sum_score = arma::sum(Score);
    double sum_diag_phi = arma::trace(Phi);
    if (sum_diag_phi > 0.0) {
        result.beta_Burden = sum_score / sum_diag_phi;
    }

    // ----------------------------------------------------------------
    // SKAT test (rho = 0): Q_SKAT = Score' * Score
    // Under H0, Q_SKAT ~ sum(lambda_j * chi2(1))
    // where lambda_j = eigenvalues of Phi
    // ----------------------------------------------------------------
    {
        double Q_SKAT = arma::dot(Score, Score);
        arma::vec eigvals;
        bool ok = arma::eig_sym(eigvals, Phi);
        if (ok) {
            arma::vec lam = eigvals(arma::find(eigvals > 1e-10));
            if (lam.n_elem > 0) {
                result.pvalue_SKAT = davies_pvalue(Q_SKAT, lam);
                if (std::isnan(result.pvalue_SKAT) || result.pvalue_SKAT < 0.0) {
                    result.pvalue_SKAT = liu_pvalue(Q_SKAT, lam);
                }
            } else {
                result.pvalue_SKAT = 1.0;
            }
        }
    }

    // ----------------------------------------------------------------
    // Burden test (rho = 1): Q_Burden = (sum(Score))^2
    // Under H0, Q_Burden / sum(Phi) ~ chi2(1)
    // ----------------------------------------------------------------
    {
        double Q_Burden = sum_score * sum_score;
        double Var_Burden = arma::accu(Phi);  // sum of all elements

        if (Var_Burden > 0.0) {
            try {
                boost::math::chi_squared chi2_1(1.0);
                result.pvalue_Burden = boost::math::cdf(
                    boost::math::complement(chi2_1, Q_Burden / Var_Burden));
            } catch (...) {
                result.pvalue_Burden = std::numeric_limits<double>::quiet_NaN();
            }
        } else {
            result.pvalue_Burden = 1.0;
        }
    }

    // ----------------------------------------------------------------
    // Burden SE: SE = |BETA| / |qnorm(p/2)|
    // Matches R: SE_Burden = abs(BETA_Burden/qnorm(p_burden/2))
    // ----------------------------------------------------------------
    if (std::isfinite(result.pvalue_Burden) && result.pvalue_Burden > 0.0 &&
        std::isfinite(result.beta_Burden)) {
        try {
            boost::math::normal norm(0, 1);
            double z = boost::math::quantile(norm, result.pvalue_Burden / 2.0);
            if (z != 0.0) {
                result.se_Burden = std::abs(result.beta_Burden / z);
            }
        } catch (...) {
            result.se_Burden = std::numeric_limits<double>::quiet_NaN();
        }
    }

    // ----------------------------------------------------------------
    // SKAT-O or single test dispatch
    // ----------------------------------------------------------------
    if (regionTestType == "BURDEN") {
        result.pvalue_SKATO = result.pvalue_Burden;
    } else if (regionTestType == "SKAT") {
        result.pvalue_SKATO = result.pvalue_SKAT;
    } else {
        // SKAT-O: optimal rho search

        // Default rho grid: {0, 0.1^2, 0.2^2, ..., 0.9^2, 1}
        arma::vec rho_vec;
        if (r_corr.n_elem > 0) {
            rho_vec = r_corr;
        } else {
            rho_vec = arma::vec({0, 0.01, 0.04, 0.09, 0.16, 0.25, 0.36, 0.49, 0.64, 0.81, 1.0});
        }

        int n_rho = (int)rho_vec.n_elem;
        arma::vec pval_each(n_rho);
        arma::vec q_each(n_rho);

        for (int k = 0; k < n_rho; k++) {
            double rho = rho_vec(k);

            // Q(rho) = (1-rho) * Score'*Score + rho * (sum(Score))^2
            double Q_rho = (1.0 - rho) * arma::dot(Score, Score) + rho * sum_score * sum_score;
            q_each(k) = Q_rho;

            if (std::abs(rho) < 1e-10) {
                pval_each(k) = result.pvalue_SKAT;
            } else if (std::abs(rho - 1.0) < 1e-10) {
                pval_each(k) = result.pvalue_Burden;
            } else {
                arma::vec ev_rho = compute_phi_rho_eigenvalues(Phi, rho, m);
                if (ev_rho.n_elem > 0) {
                    pval_each(k) = davies_pvalue(Q_rho, ev_rho);
                    if (std::isnan(pval_each(k)) || pval_each(k) < 0.0) {
                        pval_each(k) = liu_pvalue(Q_rho, ev_rho);
                    }
                } else {
                    pval_each(k) = 1.0;
                }
            }

            if (std::isnan(pval_each(k))) {
                pval_each(k) = 1.0;
            }
        }

        // Check if both rho=0 and rho=1 are in the grid
        bool has_rho0 = false, has_rho1 = false;
        for (int k = 0; k < n_rho; k++) {
            if (std::abs(rho_vec(k)) < 1e-10) has_rho0 = true;
            if (std::abs(rho_vec(k) - 1.0) < 1e-10) has_rho1 = true;
        }

        if (!has_rho0 || !has_rho1) {
            // Edge case: grid does not contain both endpoints
            double min_pval = arma::min(pval_each);
            result.pvalue_SKATO = min_pval;
            result.pvalue_SKAT = min_pval;
            result.pvalue_Burden = min_pval;
        } else {
            // Compute optimal SKAT-O p-value
            result.pvalue_SKATO = SKATO_optimal_pvalue(Score, Phi, rho_vec, pval_each, q_each);

            // Extract SKAT and Burden p-values from the grid
            for (int k = 0; k < n_rho; k++) {
                if (std::abs(rho_vec(k)) < 1e-10) {
                    result.pvalue_SKAT = pval_each(k);
                }
                if (std::abs(rho_vec(k) - 1.0) < 1e-10) {
                    result.pvalue_Burden = pval_each(k);
                }
            }
        }

        // Recompute SE_Burden with final burden p-value
        if (std::isfinite(result.pvalue_Burden) && result.pvalue_Burden > 0.0 &&
            std::isfinite(result.beta_Burden)) {
            try {
                boost::math::normal norm(0, 1);
                double z = boost::math::quantile(norm, result.pvalue_Burden / 2.0);
                if (z != 0.0) {
                    result.se_Burden = std::abs(result.beta_Burden / z);
                }
            } catch (...) {
                result.se_Burden = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }

    return result;
}
