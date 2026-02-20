// Test harness for SPA Binary (Saddlepoint Approximation)
// Tests Korg_Binom, K1_adj_Binom, K2_Binom, getroot_K1_Binom,
// Get_Saddle_Prob_Binom, SPA_binary, and the SPA() dispatcher

#include <iostream>
#include <iomanip>
#include <armadillo>
#include "spa_binary.hpp"
#include "spa.hpp"
#include "UTIL.hpp"

int main() {
    std::cout << std::setprecision(15) << std::fixed;

    // Synthetic test data (10 elements)
    arma::vec mu = {0.3, 0.7, 0.2, 0.8, 0.5, 0.1, 0.9, 0.4, 0.6, 0.15};
    arma::vec g  = {1.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 2.0};

    // q = sum(mu * g) + 1.5  (shifted to test upper tail)
    double q_mean = arma::dot(mu, g);
    double q = q_mean + 1.5;
    // qinv = 2 * sum(mu * g) - q  (mirror for two-sided test)
    double qinv = 2.0 * q_mean - q;

    std::cout << "mu*g_mean: " << q_mean << std::endl;
    std::cout << "q: " << q << std::endl;
    std::cout << "qinv: " << qinv << std::endl;

    // Test 1: Korg_Binom at t=0.1
    {
        double val = Korg_Binom(0.1, mu, g);
        std::cout << "Korg_Binom_t0.1: " << val << std::endl;
    }

    // Test 2: K1_adj_Binom at t=0.1
    {
        double val = K1_adj_Binom(0.1, mu, g, q);
        std::cout << "K1_adj_Binom_t0.1: " << val << std::endl;
    }

    // Test 3: K2_Binom at t=0.1
    {
        double val = K2_Binom(0.1, mu, g);
        std::cout << "K2_Binom_t0.1: " << val << std::endl;
    }

    // Test 4: Korg_Binom at t=0
    {
        double val = Korg_Binom(0.0, mu, g);
        std::cout << "Korg_Binom_t0: " << val << std::endl;
    }

    // Test 5: K1_adj_Binom at t=0 (should be sum(mu*g) - q = -1.5)
    {
        double val = K1_adj_Binom(0.0, mu, g, q);
        std::cout << "K1_adj_Binom_t0: " << val << std::endl;
    }

    // Test 6: K2_Binom at t=0
    {
        double val = K2_Binom(0.0, mu, g);
        std::cout << "K2_Binom_t0: " << val << std::endl;
    }

    // Test 7: getroot_K1_Binom
    {
        double tol = 0.0001220703125;  // SAIGE default: (0.5)^20
        RootResult rr = getroot_K1_Binom(0, mu, g, q, tol);
        std::cout << "getroot_root: " << rr.root << std::endl;
        std::cout << "getroot_niter: " << rr.niter << std::endl;
        std::cout << "getroot_converge: " << rr.Isconverge << std::endl;
    }

    // Test 8: getroot for qinv
    {
        double tol = 0.0001220703125;
        RootResult rr = getroot_K1_Binom(0, mu, g, qinv, tol);
        std::cout << "getroot_qinv_root: " << rr.root << std::endl;
        std::cout << "getroot_qinv_niter: " << rr.niter << std::endl;
        std::cout << "getroot_qinv_converge: " << rr.Isconverge << std::endl;
    }

    // Test 9: Get_Saddle_Prob_Binom
    {
        double tol = 0.0001220703125;
        RootResult rr = getroot_K1_Binom(0, mu, g, q, tol);
        if (rr.Isconverge && std::isfinite(rr.root)) {
            SaddleResult sr = Get_Saddle_Prob_Binom(rr.root, mu, g, q, false);
            std::cout << "SaddleProb_pval: " << sr.pval << std::endl;
            std::cout << "SaddleProb_isSaddle: " << sr.isSaddle << std::endl;
        } else {
            std::cout << "SaddleProb_pval: inf (root not finite)" << std::endl;
            std::cout << "SaddleProb_isSaddle: 0" << std::endl;
        }
    }

    // Test 10: SPA_binary full pipeline
    {
        double tol = 0.0001220703125;
        double pval_noadj = 0.05;  // dummy normal-approx p-value
        SPAResult sr = SPA_binary(mu, g, q, qinv, pval_noadj, tol, false);
        std::cout << "SPA_binary_pval: " << sr.pvalue << std::endl;
        std::cout << "SPA_binary_converge: " << sr.Isconverge << std::endl;
    }

    // Test 11: SPA() dispatcher (binary trait)
    {
        double tol = 0.0001220703125;
        double pval_noadj = 0.05;
        double pval_out = 0;
        bool converge_out = false;
        std::string traitType = "binary";
        SPA(mu, g, q, qinv, pval_noadj, tol, false, traitType, pval_out, converge_out);
        std::cout << "SPA_dispatcher_pval: " << pval_out << std::endl;
        std::cout << "SPA_dispatcher_converge: " << converge_out << std::endl;
    }

    // Test 12: SPA_pval function
    {
        double tol = 0.0001220703125;
        double pval_noadj = 0.05;
        bool converge_out = false;
        std::string traitType = "binary";
        double pval_out = SPA_pval(mu, g, q, qinv, pval_noadj, tol, false, traitType, converge_out);
        std::cout << "SPA_pval_func: " << pval_out << std::endl;
        std::cout << "SPA_pval_converge: " << converge_out << std::endl;
    }

    return 0;
}
