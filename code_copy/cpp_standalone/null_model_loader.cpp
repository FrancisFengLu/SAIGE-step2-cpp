// Null model loader for SAIGE Step 2
// Loads Step 1 output: JSON manifest + .arma binary files + variance ratios
//
// This is a STUB implementation. The interface is complete but the
// actual loading logic needs to be implemented once the Step 1 output
// format is finalized.

#include "null_model_loader.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

// TODO: Add JSON parsing (either yaml-cpp for JSON subset, or nlohmann/json)
// #include <yaml-cpp/yaml.h>


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


void loadVarianceRatios(const std::string & filepath,
                        arma::vec & varRatio_null,
                        arma::vec & varRatio_sparse,
                        arma::vec & varRatio_null_noXadj,
                        arma::vec & cateVarRatioMinMACVecExclude,
                        arma::vec & cateVarRatioMaxMACVecInclude) {
    // TODO: Implement variance ratio loading from text file
    // The file format from SAIGE is typically:
    //   Line 1: header (optional)
    //   Subsequent lines: MAC_min MAC_max VR_null VR_sparse VR_null_noXadj
    // Or it may be a simpler single-column format.
    //
    // For now, throw to indicate this needs implementation.
    throw std::runtime_error("loadVarianceRatios: not yet implemented. File: " + filepath);
}


NullModelData loadNullModel(const std::string & model_dir,
                            const std::string & varianceRatio_file) {
    // TODO: Implement full null model loading
    //
    // Steps:
    // 1. Parse nullmodel.json from model_dir to get scalars (tau, traitType, n, p, etc.)
    // 2. Load each .arma file using loadArmaVec/loadArmaMat:
    //    - mu.arma, res.arma, y.arma, V.arma, S_a.arma (vectors)
    //    - X.arma, XVX.arma, XVX_inv.arma, XXVX_inv.arma, XV.arma, XVX_inv_XV.arma (matrices)
    // 3. Load variance ratios from varianceRatio_file
    // 4. Optionally load sparse GRM if flagSparseGRM is true
    // 5. Compute derived quantities:
    //    - mu2 = mu % (1 - mu) for binary; fill with 1/tau[0] for quantitative
    //    - tauvec = {tau0, tau1}
    //    - resout (if needed for ER)
    //
    // Example skeleton:
    //
    //   NullModelData data;
    //   YAML::Node config = YAML::LoadFile(model_dir + "/nullmodel.json");
    //   data.tau0 = config["tau"][0].as<double>();
    //   data.tau1 = config["tau"][1].as<double>();
    //   data.traitType = config["traitType"].as<std::string>();
    //   data.n = config["n"].as<int>();
    //   data.p = config["p"].as<int>();
    //   ...
    //   data.mu = loadArmaVec(model_dir + "/mu.arma");
    //   data.res = loadArmaVec(model_dir + "/res.arma");
    //   data.y = loadArmaVec(model_dir + "/y.arma");
    //   data.X = loadArmaMat(model_dir + "/X.arma");
    //   data.XVX = loadArmaMat(model_dir + "/XVX.arma");
    //   ...
    //   loadVarianceRatios(varianceRatio_file, data.varRatio_null, ...);
    //   return data;

    throw std::runtime_error("loadNullModel: not yet implemented. Model dir: " + model_dir);
}
