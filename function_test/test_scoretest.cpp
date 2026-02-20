// Test harness for SAIGEClass::scoreTestFast
// Loads Step 1 .arma files, constructs SAIGEClass, runs scoreTestFast on
// synthetic genotype vectors, and prints results for comparison with R.

#define ARMA_USE_SUPERLU 1

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "saige_test.hpp"

int main(int argc, char* argv[])
{
    // --- Paths to Step 1 output (dense, x1 covariate, quantitative trait) ---
    std::string prefix = "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/output/cpp_dense_x1_output";

    // --- Load .arma files ---
    arma::vec mu_vec, res_vec, y_vec, V_vec, S_a_vec;
    arma::mat X_mat, XVX_mat, XVX_inv_mat, XXVX_inv_mat, XV_mat, XVX_inv_XV_mat;

    std::cout << "Loading .arma files from: " << prefix << std::endl;

    bool ok = true;
    ok &= mu_vec.load(prefix + ".mu.arma");
    ok &= res_vec.load(prefix + ".res.arma");
    ok &= y_vec.load(prefix + ".y.arma");
    ok &= V_vec.load(prefix + ".V.arma");
    ok &= S_a_vec.load(prefix + ".S_a.arma");
    ok &= X_mat.load(prefix + ".X.arma");
    ok &= XVX_mat.load(prefix + ".XVX.arma");
    ok &= XVX_inv_mat.load(prefix + ".XVX_inv.arma");
    ok &= XXVX_inv_mat.load(prefix + ".XXVX_inv.arma");
    ok &= XV_mat.load(prefix + ".XV.arma");
    ok &= XVX_inv_XV_mat.load(prefix + ".XVX_inv_XV.arma");

    if (!ok) {
        std::cerr << "ERROR: Failed to load one or more .arma files" << std::endl;
        return 1;
    }

    int n = static_cast<int>(y_vec.n_elem);
    int p = static_cast<int>(XV_mat.n_rows);

    std::cout << "Loaded: n=" << n << ", p=" << p << std::endl;
    std::cout << "mu dims: " << mu_vec.n_elem << std::endl;
    std::cout << "res dims: " << res_vec.n_elem << std::endl;
    std::cout << "X dims: " << X_mat.n_rows << " x " << X_mat.n_cols << std::endl;
    std::cout << "XV dims: " << XV_mat.n_rows << " x " << XV_mat.n_cols << std::endl;
    std::cout << "XVX dims: " << XVX_mat.n_rows << " x " << XVX_mat.n_cols << std::endl;
    std::cout << "XXVX_inv dims: " << XXVX_inv_mat.n_rows << " x " << XXVX_inv_mat.n_cols << std::endl;
    std::cout << "XVX_inv_XV dims: " << XVX_inv_XV_mat.n_rows << " x " << XVX_inv_XV_mat.n_cols << std::endl;
    std::cout << "S_a dims: " << S_a_vec.n_elem << std::endl;

    // --- Trait parameters ---
    // From nullmodel.json: trait=quantitative, theta=[0.286165, 0.458186]
    double tau0 = 0.286165;
    double tau1 = 0.458186;
    arma::vec tauvec = {tau0, tau1};
    std::string traitType = "quantitative";

    // mu2 for quantitative = (1/tau[1]) * ones(N)
    arma::vec mu2_vec = arma::vec(n, arma::fill::ones) * (1.0 / tau1);

    // Variance ratio from cpp_dense_x1_vr.varianceRatio.txt: 1.1105
    double varRatio_null_val = 1.1105;
    double varRatio_null_noXadj_val = 1.10911;
    arma::vec varRatio_sparse = {varRatio_null_val};
    arma::vec varRatio_null = {varRatio_null_val};
    arma::vec varRatio_null_noXadj = {varRatio_null_noXadj_val};
    arma::vec cateVarRatioMinMACVecExclude = {0.0};
    arma::vec cateVarRatioMaxMACVecInclude = {1e10};

    double SPA_Cutoff = 2.0;
    std::string impute_method = "mean";
    bool flagSparseGRM = false;
    bool isFastTest = true;
    bool isnoadjCov = false;
    double pval_cutoff_for_fastTest = 0.05;

    // No sparse GRM
    arma::umat locationMat(2, 1, arma::fill::zeros);
    arma::vec valueVec(1, arma::fill::zeros);
    int dimNum = 0;

    bool isCondition = false;
    std::vector<uint32_t> condition_genoIndex;
    bool is_Firth_beta = false;
    double pCutoffforFirth = 0.01;
    arma::vec offset(n, arma::fill::zeros);
    arma::vec resout(n, arma::fill::zeros);

    // Sigma_iXXSigma_iX: 1x1 zeros -> m_isVarPsadj = false
    arma::mat Sigma_iXXSigma_iX(1, 1, arma::fill::zeros);

    // --- Construct SAIGEClass ---
    std::cout << "\nConstructing SAIGEClass..." << std::endl;

    SAIGE::SAIGEClass saige(
        XVX_mat,
        XXVX_inv_mat,
        XV_mat,
        XVX_inv_XV_mat,
        Sigma_iXXSigma_iX,
        X_mat,
        S_a_vec,
        res_vec,
        mu2_vec,
        mu_vec,
        varRatio_sparse,
        varRatio_null,
        varRatio_null_noXadj,
        cateVarRatioMinMACVecExclude,
        cateVarRatioMaxMACVecInclude,
        SPA_Cutoff,
        tauvec,
        traitType,
        y_vec,
        impute_method,
        flagSparseGRM,
        isFastTest,
        isnoadjCov,
        pval_cutoff_for_fastTest,
        locationMat,
        valueVec,
        dimNum,
        isCondition,
        condition_genoIndex,
        is_Firth_beta,
        pCutoffforFirth,
        offset,
        resout
    );

    // Set flagSparseGRM_cur = false (no sparse GRM path)
    saige.set_flagSparseGRM_cur(false);
    // Assign the variance ratio
    saige.assignSingleVarianceRatio(false);  // use null (not sparse) variance ratio

    std::cout << "SAIGEClass constructed. varRatioVal = " << saige.m_varRatioVal << std::endl;

    // --- Create test genotype vectors ---
    // G1: sparse -- zeros except indices 0,5,10 = 1.0
    arma::vec G1(n, arma::fill::zeros);
    G1(0) = 1.0; G1(5) = 1.0; G1(10) = 1.0;

    // G2: first 100 samples = 1.0, rest = 0.0
    arma::vec G2(n, arma::fill::zeros);
    for (int i = 0; i < 100; i++) G2(i) = 1.0;

    // G3: alternating 0, 1, 0, 1, ...
    arma::vec G3(n, arma::fill::zeros);
    for (int i = 0; i < n; i++) {
        if (i % 2 == 1) G3(i) = 1.0;
    }

    std::vector<std::string> gnames = {"G1_sparse", "G2_block100", "G3_alternating"};
    std::vector<arma::vec*> gvecs = {&G1, &G2, &G3};

    std::cout << std::setprecision(15) << std::fixed;

    // Open output file
    std::ofstream out("/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Step_2_Feb_11/function_test/cpp_scoretest_results.txt");
    out << std::setprecision(15) << std::fixed;

    for (size_t gi = 0; gi < gvecs.size(); gi++) {
        arma::vec& G = *gvecs[gi];

        // Find nonzero indices
        arma::uvec indexNonZero = arma::find(G != 0);

        double Beta = 0, seBeta = 0, pval = 0, Tstat = 0, var1 = 0, var2 = 0;
        std::string pval_str;
        bool islogp = false;
        double altFreq = arma::sum(G) / (2.0 * n);

        // --- scoreTestFast ---
        saige.scoreTestFast(G, indexNonZero, Beta, seBeta, pval_str, pval, islogp, altFreq, Tstat, var1, var2);

        std::cout << "\n=== " << gnames[gi] << " ===" << std::endl;
        std::cout << "  Tstat (S)   = " << Tstat << std::endl;
        std::cout << "  var1        = " << var1 << std::endl;
        std::cout << "  var2        = " << var2 << std::endl;
        std::cout << "  Beta        = " << Beta << std::endl;
        std::cout << "  seBeta      = " << seBeta << std::endl;
        std::cout << "  pval        = " << pval << std::endl;
        std::cout << "  pval_str    = " << pval_str << std::endl;
        std::cout << "  islogp      = " << islogp << std::endl;

        out << gnames[gi] << "\t" << Tstat << "\t" << var1 << "\t" << var2
            << "\t" << Beta << "\t" << seBeta << "\t" << pval << "\t" << pval_str
            << "\t" << islogp << std::endl;
    }

    // --- Also test scoreTest (non-fast, full version) for verification ---
    std::cout << "\n\n=== scoreTest (non-fast) for comparison ===" << std::endl;
    std::ofstream out2("/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Step_2_Feb_11/function_test/cpp_scoretest_nonfastresults.txt");
    out2 << std::setprecision(15) << std::fixed;

    for (size_t gi = 0; gi < gvecs.size(); gi++) {
        arma::vec& G = *gvecs[gi];
        arma::uvec indexNonZero = arma::find(G != 0);

        double Beta = 0, seBeta = 0, pval = 0, Tstat = 0, var1 = 0, var2 = 0, gy = 0;
        std::string pval_str;
        bool islogp = false;
        double altFreq = arma::sum(G) / (2.0 * n);
        arma::vec gtilde(n, arma::fill::zeros);
        arma::vec P2Vec(n, arma::fill::zeros);
        bool is_region = false;

        saige.scoreTest(G, Beta, seBeta, pval_str, pval, islogp, altFreq, Tstat, var1, var2, gtilde, P2Vec, gy, is_region, indexNonZero);

        std::cout << "\n=== " << gnames[gi] << " (non-fast) ===" << std::endl;
        std::cout << "  Tstat (S)   = " << Tstat << std::endl;
        std::cout << "  var1        = " << var1 << std::endl;
        std::cout << "  var2        = " << var2 << std::endl;
        std::cout << "  Beta        = " << Beta << std::endl;
        std::cout << "  seBeta      = " << seBeta << std::endl;
        std::cout << "  pval        = " << pval << std::endl;
        std::cout << "  pval_str    = " << pval_str << std::endl;

        out2 << gnames[gi] << "\t" << Tstat << "\t" << var1 << "\t" << var2
             << "\t" << Beta << "\t" << seBeta << "\t" << pval << "\t" << pval_str
             << "\t" << islogp << std::endl;
    }

    out.close();
    out2.close();

    std::cout << "\nC++ results written to cpp_scoretest_results.txt and cpp_scoretest_nonfastresults.txt" << std::endl;
    return 0;
}
