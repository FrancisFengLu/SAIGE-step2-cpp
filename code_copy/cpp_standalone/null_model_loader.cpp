// Null model loader for SAIGE Step 2
// Loads Step 1 output: JSON manifest + .arma binary files + variance ratios
//
// Phase 1 implementation: full loading of nullmodel.json, .arma files,
// variance ratios, sparse GRM (optional), and derived quantities.

#include "null_model_loader.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <filesystem>
#include <iomanip>
#include <yaml-cpp/yaml.h>

extern bool g_writeCheckpoints;
extern std::string g_checkpointDir;


arma::vec loadArmaVec(const std::string & filepath) {
    arma::vec v;
    bool ok = v.load(filepath, arma::arma_binary);
    if (!ok) {
        throw std::runtime_error("Failed to load arma::vec from: " + filepath);
    }
    return v;
}

arma::mat loadArmaMat(const std::string & filepath) {
    arma::mat m;
    bool ok = m.load(filepath, arma::arma_binary);
    if (!ok) {
        throw std::runtime_error("Failed to load arma::mat from: " + filepath);
    }
    return m;
}

// Try to load an arma::vec; return true if successful, false if file doesn't exist.
// Throws only if file exists but fails to load.
static bool tryLoadArmaVec(const std::string & filepath, arma::vec & out) {
    if (!std::filesystem::exists(filepath)) {
        return false;
    }
    bool ok = out.load(filepath, arma::arma_binary);
    if (!ok) {
        throw std::runtime_error("File exists but failed to load arma::vec from: " + filepath);
    }
    return true;
}

// Try to load an arma::mat; return true if successful, false if file doesn't exist.
static bool tryLoadArmaMat(const std::string & filepath, arma::mat & out) {
    if (!std::filesystem::exists(filepath)) {
        return false;
    }
    bool ok = out.load(filepath, arma::arma_binary);
    if (!ok) {
        throw std::runtime_error("File exists but failed to load arma::mat from: " + filepath);
    }
    return true;
}

// Try to load an arma::umat; return true if successful, false if file doesn't exist.
static bool tryLoadArmaUmat(const std::string & filepath, arma::umat & out) {
    if (!std::filesystem::exists(filepath)) {
        return false;
    }
    bool ok = out.load(filepath, arma::arma_binary);
    if (!ok) {
        throw std::runtime_error("File exists but failed to load arma::umat from: " + filepath);
    }
    return true;
}


void loadVarianceRatios(const std::string & filepath,
                        arma::vec & varRatio_null,
                        arma::vec & varRatio_sparse,
                        arma::vec & varRatio_null_noXadj,
                        arma::vec & cateVarRatioMinMACVecExclude,
                        arma::vec & cateVarRatioMaxMACVecInclude) {
    std::ifstream infile(filepath);
    if (!infile.is_open()) {
        throw std::runtime_error("loadVarianceRatios: cannot open file: " + filepath);
    }

    std::cout << "  Loading variance ratios from: " << filepath << std::endl;

    // Temporary vectors to accumulate values
    std::vector<double> vr_null_vec;
    std::vector<double> vr_sparse_vec;
    std::vector<double> mac_min_vec;
    std::vector<double> mac_max_vec;
    std::vector<double> vr_noXadj_vec;

    std::string line;
    bool headerSkipped = false;

    while (std::getline(infile, line)) {
        // Skip empty lines
        if (line.empty()) continue;

        // Skip header line: detect by checking if first token is non-numeric
        if (!headerSkipped) {
            // Try to parse the first token as a number
            std::istringstream testss(line);
            std::string firstToken;
            testss >> firstToken;
            // If the first token starts with a letter, it's a header
            if (!firstToken.empty() && (std::isalpha(firstToken[0]) || firstToken[0] == '#')) {
                headerSkipped = true;
                continue;
            }
            headerSkipped = true;
            // Fall through to parse this line as data
        }

        std::istringstream ss(line);
        double col1, col2, col3, col4, col5;

        // Column 1: varianceRatio (null)
        if (!(ss >> col1)) continue;

        // Column 2: varianceRatioSparse
        if (!(ss >> col2)) {
            // Only 1 column: use col1 for everything
            vr_null_vec.push_back(col1);
            vr_sparse_vec.push_back(col1);
            mac_min_vec.push_back(0.0);
            mac_max_vec.push_back(1e10);
            vr_noXadj_vec.push_back(col1);
            continue;
        }

        // Column 3: MAC_cate_min
        if (!(ss >> col3)) {
            // 2 columns only
            vr_null_vec.push_back(col1);
            vr_sparse_vec.push_back(col2);
            mac_min_vec.push_back(0.0);
            mac_max_vec.push_back(1e10);
            vr_noXadj_vec.push_back(col1);
            continue;
        }

        // Column 4: MAC_cate_max
        if (!(ss >> col4)) {
            // 3 columns
            vr_null_vec.push_back(col1);
            vr_sparse_vec.push_back(col2);
            mac_min_vec.push_back(col3);
            mac_max_vec.push_back(1e10);
            vr_noXadj_vec.push_back(col1);
            continue;
        }

        // Column 5: noXadj_varianceRatio (optional)
        if (!(ss >> col5)) {
            // 4 columns: no noXadj column
            col5 = col1; // default to same as col1
        }

        vr_null_vec.push_back(col1);
        vr_sparse_vec.push_back(col2);
        mac_min_vec.push_back(col3);
        mac_max_vec.push_back(col4);
        vr_noXadj_vec.push_back(col5);
    }

    infile.close();

    if (vr_null_vec.empty()) {
        throw std::runtime_error("loadVarianceRatios: no data rows found in: " + filepath);
    }

    // Convert std::vector<double> to arma::vec
    size_t nrows = vr_null_vec.size();
    varRatio_null.set_size(nrows);
    varRatio_sparse.set_size(nrows);
    cateVarRatioMinMACVecExclude.set_size(nrows);
    cateVarRatioMaxMACVecInclude.set_size(nrows);
    varRatio_null_noXadj.set_size(nrows);

    for (size_t i = 0; i < nrows; i++) {
        varRatio_null(i)                    = vr_null_vec[i];
        varRatio_sparse(i)                  = vr_sparse_vec[i];
        cateVarRatioMinMACVecExclude(i)     = mac_min_vec[i];
        cateVarRatioMaxMACVecInclude(i)     = mac_max_vec[i];
        varRatio_null_noXadj(i)             = vr_noXadj_vec[i];
    }

    std::cout << "  Loaded " << nrows << " variance ratio categories" << std::endl;
    for (size_t i = 0; i < nrows; i++) {
        std::cout << "    VR[" << i << "]: null=" << varRatio_null(i)
                  << "  sparse=" << varRatio_sparse(i)
                  << "  MAC=[" << cateVarRatioMinMACVecExclude(i)
                  << ", " << cateVarRatioMaxMACVecInclude(i)
                  << "]  noXadj=" << varRatio_null_noXadj(i)
                  << std::endl;
    }
}


NullModelData loadNullModel(const std::string & model_dir,
                            const std::string & varianceRatio_file) {

    NullModelData data;

    std::cout << "========================================" << std::endl;
    std::cout << "Loading null model from: " << model_dir << std::endl;
    std::cout << "========================================" << std::endl;

    // =========================================================================
    // 1. Parse nullmodel.json from model_dir using yaml-cpp (yaml-cpp can parse JSON)
    // =========================================================================
    std::string json_path = model_dir + "/nullmodel.json";
    std::cout << "\n[1/5] Parsing JSON manifest: " << json_path << std::endl;

    if (!std::filesystem::exists(json_path)) {
        throw std::runtime_error("nullmodel.json not found at: " + json_path);
    }

    YAML::Node config = YAML::LoadFile(json_path);

    // --- tau (array of 2) ---
    if (config["tau"]) {
        YAML::Node tau_node = config["tau"];
        if (tau_node.IsSequence() && tau_node.size() >= 2) {
            data.tau0 = tau_node[0].as<double>();
            data.tau1 = tau_node[1].as<double>();
        } else if (tau_node.IsScalar()) {
            data.tau0 = tau_node.as<double>();
            data.tau1 = 0.0;
        } else {
            throw std::runtime_error("Invalid tau format in nullmodel.json");
        }
    } else {
        throw std::runtime_error("Missing 'tau' in nullmodel.json");
    }
    std::cout << "  tau = [" << data.tau0 << ", " << data.tau1 << "]" << std::endl;

    // --- traitType ---
    if (config["traitType"]) {
        data.traitType = config["traitType"].as<std::string>();
    } else {
        // Fallback: check for "trait" key
        if (config["trait"]) {
            data.traitType = config["trait"].as<std::string>();
        } else {
            throw std::runtime_error("Missing 'traitType' in nullmodel.json");
        }
    }
    std::cout << "  traitType = " << data.traitType << std::endl;

    // --- n (sample size) ---
    if (config["n"]) {
        data.n = config["n"].as<int>();
    } else {
        // Will be set after loading vectors
        data.n = 0;
    }

    // --- p (number of covariates including intercept) ---
    if (config["p"]) {
        data.p = config["p"].as<int>();
    } else {
        // Will be set after loading matrices
        data.p = 0;
    }

    // --- SPA_Cutoff (default 2.0) ---
    if (config["SPA_Cutoff"]) {
        data.SPA_Cutoff = config["SPA_Cutoff"].as<double>();
    } else {
        data.SPA_Cutoff = 2.0;
    }

    // --- impute_method ---
    if (config["impute_method"]) {
        data.impute_method = config["impute_method"].as<std::string>();
    } else {
        data.impute_method = "mean";
    }

    // --- flagSparseGRM ---
    if (config["flagSparseGRM"]) {
        data.flagSparseGRM = config["flagSparseGRM"].as<bool>();
    } else {
        data.flagSparseGRM = false;
    }

    // --- isFastTest ---
    if (config["isFastTest"]) {
        data.isFastTest = config["isFastTest"].as<bool>();
    } else {
        data.isFastTest = true;
    }

    // --- isnoadjCov ---
    if (config["isnoadjCov"]) {
        data.isnoadjCov = config["isnoadjCov"].as<bool>();
    } else {
        data.isnoadjCov = false;
    }

    // --- pval_cutoff_for_fastTest ---
    if (config["pval_cutoff_for_fastTest"]) {
        data.pval_cutoff_for_fastTest = config["pval_cutoff_for_fastTest"].as<double>();
    } else {
        data.pval_cutoff_for_fastTest = 0.05;
    }

    // --- isCondition ---
    if (config["isCondition"]) {
        data.isCondition = config["isCondition"].as<bool>();
    } else {
        data.isCondition = false;
    }

    // --- is_Firth_beta ---
    if (config["is_Firth_beta"]) {
        data.is_Firth_beta = config["is_Firth_beta"].as<bool>();
    } else {
        data.is_Firth_beta = false;
    }

    // --- pCutoffforFirth ---
    if (config["pCutoffforFirth"]) {
        data.pCutoffforFirth = config["pCutoffforFirth"].as<double>();
    } else {
        data.pCutoffforFirth = 0.01;
    }

    // --- sampleIDs (array of strings) ---
    if (config["sampleIDs"]) {
        YAML::Node ids = config["sampleIDs"];
        if (ids.IsSequence()) {
            for (size_t i = 0; i < ids.size(); i++) {
                data.sampleIDs.push_back(ids[i].as<std::string>());
            }
        }
    }

    std::cout << "  n = " << data.n << ", p = " << data.p << std::endl;
    std::cout << "  SPA_Cutoff = " << data.SPA_Cutoff << std::endl;
    std::cout << "  impute_method = " << data.impute_method << std::endl;
    std::cout << "  flagSparseGRM = " << data.flagSparseGRM << std::endl;
    std::cout << "  isFastTest = " << data.isFastTest << std::endl;
    std::cout << "  isnoadjCov = " << data.isnoadjCov << std::endl;
    std::cout << "  isCondition = " << data.isCondition << std::endl;
    std::cout << "  is_Firth_beta = " << data.is_Firth_beta << std::endl;
    std::cout << "  pCutoffforFirth = " << data.pCutoffforFirth << std::endl;
    std::cout << "  sampleIDs: " << data.sampleIDs.size() << " entries" << std::endl;

    // =========================================================================
    // 2. Load .arma vectors
    // =========================================================================
    std::cout << "\n[2/5] Loading .arma vector files..." << std::endl;

    data.mu  = loadArmaVec(model_dir + "/mu.arma");
    std::cout << "  mu:  " << data.mu.n_elem << " elements" << std::endl;

    data.res = loadArmaVec(model_dir + "/res.arma");
    std::cout << "  res: " << data.res.n_elem << " elements" << std::endl;

    data.y   = loadArmaVec(model_dir + "/y.arma");
    std::cout << "  y:   " << data.y.n_elem << " elements" << std::endl;

    data.V   = loadArmaVec(model_dir + "/V.arma");
    std::cout << "  V:   " << data.V.n_elem << " elements" << std::endl;

    data.S_a = loadArmaVec(model_dir + "/S_a.arma");
    std::cout << "  S_a: " << data.S_a.n_elem << " elements" << std::endl;

    // Update n from vector size if not set from JSON
    if (data.n == 0) {
        data.n = static_cast<int>(data.y.n_elem);
        std::cout << "  [inferred n = " << data.n << " from y.arma]" << std::endl;
    }

    // =========================================================================
    // 3. Load .arma matrices
    // =========================================================================
    std::cout << "\n[3/5] Loading .arma matrix files..." << std::endl;

    data.X          = loadArmaMat(model_dir + "/X.arma");
    std::cout << "  X:          " << data.X.n_rows << " x " << data.X.n_cols << std::endl;

    data.XVX        = loadArmaMat(model_dir + "/XVX.arma");
    std::cout << "  XVX:        " << data.XVX.n_rows << " x " << data.XVX.n_cols << std::endl;

    data.XVX_inv    = loadArmaMat(model_dir + "/XVX_inv.arma");
    std::cout << "  XVX_inv:    " << data.XVX_inv.n_rows << " x " << data.XVX_inv.n_cols << std::endl;

    data.XXVX_inv   = loadArmaMat(model_dir + "/XXVX_inv.arma");
    std::cout << "  XXVX_inv:   " << data.XXVX_inv.n_rows << " x " << data.XXVX_inv.n_cols << std::endl;

    data.XV         = loadArmaMat(model_dir + "/XV.arma");
    std::cout << "  XV:         " << data.XV.n_rows << " x " << data.XV.n_cols << std::endl;

    data.XVX_inv_XV = loadArmaMat(model_dir + "/XVX_inv_XV.arma");
    std::cout << "  XVX_inv_XV: " << data.XVX_inv_XV.n_rows << " x " << data.XVX_inv_XV.n_cols << std::endl;

    // Update p from matrix dimensions if not set from JSON
    if (data.p == 0) {
        data.p = static_cast<int>(data.XV.n_rows);
        std::cout << "  [inferred p = " << data.p << " from XV.arma]" << std::endl;
    }

    // --- Sigma_iXXSigma_iX: optional, default to 1x1 zero matrix ---
    if (!tryLoadArmaMat(model_dir + "/Sigma_iXXSigma_iX.arma", data.Sigma_iXXSigma_iX)) {
        data.Sigma_iXXSigma_iX = arma::mat(1, 1, arma::fill::zeros);
        std::cout << "  Sigma_iXXSigma_iX: not found, using 1x1 zero (dummy)" << std::endl;
    } else {
        std::cout << "  Sigma_iXXSigma_iX: " << data.Sigma_iXXSigma_iX.n_rows
                  << " x " << data.Sigma_iXXSigma_iX.n_cols << std::endl;
    }

    // =========================================================================
    // 4. Compute derived quantities
    // =========================================================================
    std::cout << "\n[4/5] Computing derived quantities..." << std::endl;

    // --- tauvec ---
    data.tauvec = arma::vec(2);
    data.tauvec(0) = data.tau0;
    data.tauvec(1) = data.tau1;
    std::cout << "  tauvec = [" << data.tauvec(0) << ", " << data.tauvec(1) << "]" << std::endl;

    // --- mu2: mu*(1-mu) for binary, 1/tau[0] for quantitative ---
    int n = data.n;
    if (data.traitType == "binary") {
        data.mu2 = data.mu % (1.0 - data.mu);
        std::cout << "  mu2 = mu*(1-mu) [binary], range: ["
                  << data.mu2.min() << ", " << data.mu2.max() << "]" << std::endl;
    } else {
        // quantitative: mu2 = 1/tau[0] for all elements
        // Note: In the test_scoretest.cpp, quantitative uses 1/tau[1] (the genetic tau).
        // But in the SAIGE R code (SAIGE_test.cpp constructor), for quantitative:
        //   mu2 is passed in from R, where it's set to mu*(1-mu) for binary
        //   and for quantitative the variance weight is 1/tau[0].
        // The CLAUDE.md says: "fill(1.0/tau0) for quantitative"
        data.mu2 = arma::vec(n, arma::fill::ones) * (1.0 / data.tau0);
        std::cout << "  mu2 = 1/tau[0] = " << (1.0 / data.tau0)
                  << " [quantitative, all " << n << " elements]" << std::endl;
    }

    // --- offset: try to load, else zeros ---
    if (!tryLoadArmaVec(model_dir + "/offset.arma", data.offset)) {
        data.offset = arma::vec(n, arma::fill::zeros);
        std::cout << "  offset: not found, using zeros(" << n << ")" << std::endl;
    } else {
        std::cout << "  offset: " << data.offset.n_elem << " elements" << std::endl;
    }

    // --- resout: try to load, else copy from res ---
    if (!tryLoadArmaVec(model_dir + "/resout.arma", data.resout)) {
        data.resout = data.res;
        std::cout << "  resout: not found, copying from res (" << data.resout.n_elem << " elements)" << std::endl;
    } else {
        std::cout << "  resout: " << data.resout.n_elem << " elements" << std::endl;
    }

    // =========================================================================
    // 5. Load variance ratios
    // =========================================================================
    std::cout << "\n[5/5] Loading variance ratios..." << std::endl;

    loadVarianceRatios(varianceRatio_file,
                       data.varRatio_null,
                       data.varRatio_sparse,
                       data.varRatio_null_noXadj,
                       data.cateVarRatioMinMACVecExclude,
                       data.cateVarRatioMaxMACVecInclude);

    // =========================================================================
    // Optional: Load sparse GRM
    // =========================================================================
    data.dimNum = 0;
    if (data.flagSparseGRM) {
        std::cout << "\n[Optional] Loading sparse GRM files..." << std::endl;

        bool locLoaded = tryLoadArmaUmat(model_dir + "/sparseGRM_locationMat.arma", data.locationMat);
        bool valLoaded = tryLoadArmaVec(model_dir + "/sparseGRM_valueVec.arma", data.valueVec);

        if (locLoaded && valLoaded) {
            data.dimNum = data.n;
            std::cout << "  sparseGRM_locationMat: " << data.locationMat.n_rows
                      << " x " << data.locationMat.n_cols << std::endl;
            std::cout << "  sparseGRM_valueVec: " << data.valueVec.n_elem << " elements" << std::endl;
            std::cout << "  dimNum = " << data.dimNum << std::endl;
        } else {
            std::cout << "  WARNING: flagSparseGRM=true but sparse GRM files not found." << std::endl;
            std::cout << "  Setting dimNum=0 (no sparse GRM)." << std::endl;
            data.locationMat = arma::umat(2, 1, arma::fill::zeros);
            data.valueVec = arma::vec(1, arma::fill::zeros);
            data.dimNum = 0;
        }
    } else {
        // No sparse GRM: set dummy values
        data.locationMat = arma::umat(2, 1, arma::fill::zeros);
        data.valueVec = arma::vec(1, arma::fill::zeros);
        std::cout << "\n  Sparse GRM: not requested (flagSparseGRM=false)" << std::endl;
    }

    // =========================================================================
    // Summary
    // =========================================================================
    std::cout << "\n========================================" << std::endl;
    std::cout << "Null model loaded successfully." << std::endl;
    std::cout << "  n=" << data.n << "  p=" << data.p
              << "  traitType=" << data.traitType << std::endl;
    std::cout << "  tau=[" << data.tau0 << ", " << data.tau1 << "]" << std::endl;
    std::cout << "  VR categories: " << data.varRatio_null.n_elem << std::endl;
    std::cout << "  Sparse GRM: " << (data.dimNum > 0 ? "YES" : "NO") << std::endl;
    std::cout << "========================================" << std::endl;

    // === CHECKPOINT OUTPUT BEGIN ===
    if (g_writeCheckpoints && !g_checkpointDir.empty()) {
        // Create checkpoint directory if it doesn't exist
        std::filesystem::create_directories(g_checkpointDir);

        std::cout << "\nWriting checkpoint files to: " << g_checkpointDir << std::endl;

        // --- ckpt_01_null_model.txt ---
        {
            std::ofstream ofs(g_checkpointDir + "/ckpt_01_null_model.txt");
            if (ofs.is_open()) {
                ofs << std::setprecision(15);
                ofs << "parameter\tvalue" << std::endl;
                ofs << "n\t" << data.n << std::endl;
                ofs << "p\t" << data.p << std::endl;
                ofs << "traitType\t" << data.traitType << std::endl;
                ofs << "tau0\t" << data.tau0 << std::endl;
                ofs << "tau1\t" << data.tau1 << std::endl;
                ofs << "SPA_Cutoff\t" << data.SPA_Cutoff << std::endl;
                ofs.close();
                std::cout << "  Written: ckpt_01_null_model.txt" << std::endl;
            } else {
                std::cerr << "  WARNING: Could not write ckpt_01_null_model.txt" << std::endl;
            }
        }

        // --- ckpt_02_mu.txt ---
        {
            std::ofstream ofs(g_checkpointDir + "/ckpt_02_mu.txt");
            if (ofs.is_open()) {
                ofs << std::setprecision(15);
                ofs << "index\tmu" << std::endl;
                int count = std::min((int)data.mu.n_elem, 10);
                for (int i = 0; i < count; i++) {
                    ofs << (i + 1) << "\t" << data.mu(i) << std::endl;
                }
                ofs.close();
                std::cout << "  Written: ckpt_02_mu.txt" << std::endl;
            } else {
                std::cerr << "  WARNING: Could not write ckpt_02_mu.txt" << std::endl;
            }
        }

        // --- ckpt_03_res.txt ---
        {
            std::ofstream ofs(g_checkpointDir + "/ckpt_03_res.txt");
            if (ofs.is_open()) {
                ofs << std::setprecision(15);
                ofs << "index\tres" << std::endl;
                int count = std::min((int)data.res.n_elem, 10);
                for (int i = 0; i < count; i++) {
                    ofs << (i + 1) << "\t" << data.res(i) << std::endl;
                }
                ofs.close();
                std::cout << "  Written: ckpt_03_res.txt" << std::endl;
            } else {
                std::cerr << "  WARNING: Could not write ckpt_03_res.txt" << std::endl;
            }
        }

        // --- ckpt_04_XVX.txt ---
        {
            std::ofstream ofs(g_checkpointDir + "/ckpt_04_XVX.txt");
            if (ofs.is_open()) {
                ofs << std::setprecision(15);
                ofs << "# XVX matrix [" << data.XVX.n_rows << " x " << data.XVX.n_cols << "]" << std::endl;
                for (arma::uword r = 0; r < data.XVX.n_rows; r++) {
                    for (arma::uword c = 0; c < data.XVX.n_cols; c++) {
                        if (c > 0) ofs << "\t";
                        ofs << data.XVX(r, c);
                    }
                    ofs << std::endl;
                }
                ofs.close();
                std::cout << "  Written: ckpt_04_XVX.txt" << std::endl;
            } else {
                std::cerr << "  WARNING: Could not write ckpt_04_XVX.txt" << std::endl;
            }
        }

        // --- ckpt_05_variance_ratios.txt ---
        {
            std::ofstream ofs(g_checkpointDir + "/ckpt_05_variance_ratios.txt");
            if (ofs.is_open()) {
                ofs << std::setprecision(15);
                ofs << "index\tvarRatio_null\tvarRatio_sparse\tMAC_min\tMAC_max\tvarRatio_noXadj" << std::endl;
                for (arma::uword i = 0; i < data.varRatio_null.n_elem; i++) {
                    ofs << (i + 1) << "\t"
                        << data.varRatio_null(i) << "\t"
                        << data.varRatio_sparse(i) << "\t"
                        << data.cateVarRatioMinMACVecExclude(i) << "\t"
                        << data.cateVarRatioMaxMACVecInclude(i) << "\t"
                        << data.varRatio_null_noXadj(i) << std::endl;
                }
                ofs.close();
                std::cout << "  Written: ckpt_05_variance_ratios.txt" << std::endl;
            } else {
                std::cerr << "  WARNING: Could not write ckpt_05_variance_ratios.txt" << std::endl;
            }
        }

        std::cout << "Checkpoint output complete." << std::endl;
    }
    // === CHECKPOINT OUTPUT END ===

    return data;
}
