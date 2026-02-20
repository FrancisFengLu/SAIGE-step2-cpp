// Test harness for CCT (Cauchy Combination Test)
// Compares C++ standalone output against R's SAIGE::CCT

#include <iostream>
#include <iomanip>
#include <vector>
#include <armadillo>
#include "cct.hpp"

int main() {
    std::cout << std::setprecision(15) << std::fixed;

    // Test case 1: {0.01, 0.05, 0.1, 0.5}
    {
        arma::vec pvals = {0.01, 0.05, 0.1, 0.5};
        double result = CCT_cpp(pvals);
        std::cout << "CCT_test1: " << result << std::endl;
    }

    // Test case 2: {1e-5, 0.3, 0.7, 0.9}
    {
        arma::vec pvals = {1e-5, 0.3, 0.7, 0.9};
        double result = CCT_cpp(pvals);
        std::cout << "CCT_test2: " << result << std::endl;
    }

    // Test case 3: {0.5, 0.5, 0.5, 0.5}
    {
        arma::vec pvals = {0.5, 0.5, 0.5, 0.5};
        double result = CCT_cpp(pvals);
        std::cout << "CCT_test3: " << result << std::endl;
    }

    // Test case 4: {1e-10, 1e-8, 1e-6} (extreme small)
    {
        arma::vec pvals = {1e-10, 1e-8, 1e-6};
        double result = CCT_cpp(pvals);
        std::cout << "CCT_test4: " << result << std::endl;
    }

    // Test case 5: {0.001} (single value)
    {
        arma::vec pvals = {0.001};
        double result = CCT_cpp(pvals);
        std::cout << "CCT_test5: " << result << std::endl;
    }

    return 0;
}
