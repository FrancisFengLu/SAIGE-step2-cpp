// Test harness for UTIL::getWeights
// Compares C++ standalone dbeta output against R's dbeta

#include <iostream>
#include <iomanip>
#include <armadillo>
#include "UTIL.hpp"

int main() {
    std::cout << std::setprecision(15) << std::fixed;

    // MAF values to test
    std::vector<double> mafs = {0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5};

    // Default beta params: a1=1, a2=25
    arma::vec wBeta_default = {1.0, 25.0};

    std::cout << "=== getWeights with Beta(1, 25) ===" << std::endl;
    for (double maf : mafs) {
        double w = getWeights("linear.weighted", maf, wBeta_default);
        std::cout << "getWeights_beta1_25_maf" << maf << ": " << w << std::endl;
    }

    // Also test with Beta(1, 1) -- uniform
    arma::vec wBeta_uniform = {1.0, 1.0};

    std::cout << "=== getWeights with Beta(1, 1) ===" << std::endl;
    for (double maf : mafs) {
        double w = getWeights("linear.weighted", maf, wBeta_uniform);
        std::cout << "getWeights_beta1_1_maf" << maf << ": " << w << std::endl;
    }

    // Test linear kernel
    std::cout << "=== getWeights with linear kernel ===" << std::endl;
    for (double maf : mafs) {
        double w = getWeights("linear", maf, wBeta_default);
        std::cout << "getWeights_linear_maf" << maf << ": " << w << std::endl;
    }

    return 0;
}
