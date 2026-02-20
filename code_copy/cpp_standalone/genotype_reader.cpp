// Standalone port of SAIGE/src/PLINK.cpp
// PLINK genotype reader -- all Rcpp dependencies removed
//
// Conversions from original SAIGE PLINK.cpp:
//   1. #include <RcppArmadillo.h>  -->  #include <armadillo>
//   2. Rcpp::stop(...)             -->  throw std::runtime_error(...)
//   3. boost::split(...)           -->  splitLine() using std::istringstream
//   4. boost::replace_all(...)     -->  manual \r removal
//   5. Rcpp::match(...)            -->  std::unordered_map lookup
//   6. Rcpp::CharacterVector       -->  std::vector<std::string>
//   7. Rcpp::IntegerVector         -->  std::vector<uint32_t>

#include <armadillo>
#include "genotype_reader.hpp"
#include "UTIL.hpp"

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <stdexcept>
#include <cstdio>
#include <cstdint>
#include <sys/stat.h>
#include <filesystem>
#include <iomanip>

// Global checkpoint flag for debugging
bool g_writeCheckpoints = false;
std::string g_checkpointDir = "";

// Global pointer to PLINK object (mirrors SAIGE's ptr_gPLINKobj in Main.cpp)
PLINK::PlinkClass* ptr_gPLINKobj = nullptr;

// Static counter for checkpoint tracking
static uint64_t s_checkpointMarkerCount = 0;

// ============================================================
// Helper: split a line on whitespace (replaces boost::split)
// ============================================================
static std::vector<std::string> splitLine(const std::string& line) {
    std::vector<std::string> tokens;
    std::istringstream iss(line);
    std::string token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

// ============================================================
// Helper: remove trailing \r from string (replaces boost::replace_all for \r)
// ============================================================
static void removeTrailingCR(std::string& s) {
    if (!s.empty() && s.back() == '\r') {
        s.pop_back();
    }
}

namespace PLINK {

// ============================================================
// Constructor
// ============================================================
PlinkClass::PlinkClass(std::string t_bimFile,
                       std::string t_famFile,
                       std::string t_bedFile,
                       std::string t_AlleleOrder)
{
    setPlinkobj(t_bimFile, t_famFile, t_bedFile);
    m_AlleleOrder = t_AlleleOrder;
}

// ============================================================
// setPlinkobj: read PLINK files and open .bed
// ============================================================
void PlinkClass::setPlinkobj(std::string t_bimFile,
                             std::string t_famFile,
                             std::string t_bedFile)
{
    m_bimFile = t_bimFile;
    m_famFile = t_famFile;
    m_bedFile = t_bedFile;

    readBimFile();
    readFamFile();

    m_fin = fopen(t_bedFile.c_str(), "rb");
    if (!m_fin) {
        throw std::runtime_error("Cannot open PLINK .bed file: " + t_bedFile);
    }

    char magicNumber_2[65536];
    fread(magicNumber_2, 2, 1, m_fin);
    char magicNumber3;
    fread(&magicNumber3, 1, 1, m_fin);
    fseek(m_fin, 3, SEEK_SET);

    if (magicNumber3 != 1) {
        throw std::runtime_error(
            "The third magic number of the plink bed file is not 00000001. "
            "Please use SNP-major plink (plink version >= 1.9) files.");
    }
}

// ============================================================
// readBimFile: parse the .bim file
// ============================================================
void PlinkClass::readBimFile()
{
    std::cout << "Reading bim file...." << std::endl;
    std::ifstream bim(m_bimFile);
    if (!bim.is_open()) {
        throw std::runtime_error("Cannot open PLINK .bim file: " + m_bimFile);
    }

    m_M0 = 0;
    std::string line;

    while (getline(bim, line)) {
        m_M0++;
        std::vector<std::string> line_elements = splitLine(line);

        if (line_elements.size() < 6) {
            throw std::runtime_error(
                "PLINK .bim file has fewer than 6 columns at line " +
                std::to_string(m_M0));
        }

        // Remove trailing \r from the last element (Windows line endings)
        removeTrailingCR(line_elements[line_elements.size() - 1]);

        m_chr.push_back(line_elements[0]);
        m_MarkerInPlink.push_back(line_elements[1]);
        m_gd.push_back(std::stof(line_elements[2]));
        m_pd.push_back(std::stoi(line_elements[3]));

        // Convert alleles to uppercase (matching SAIGE behavior)
        std::transform(line_elements[4].begin(), line_elements[4].end(),
                       line_elements[4].begin(), ::toupper);
        std::transform(line_elements[5].begin(), line_elements[5].end(),
                       line_elements[5].begin(), ::toupper);

        m_alt.push_back(line_elements[4]);  // allele 1, usually minor allele, alt-first
        m_ref.push_back(line_elements[5]);  // allele 2, usually major allele
    }
    m_M = m_M0;

    std::cout << "Number of markers in bim file: " << m_M0 << std::endl;
}

// ============================================================
// readFamFile: parse the .fam file
// ============================================================
void PlinkClass::readFamFile()
{
    std::cout << "Reading fam file...." << std::endl;
    std::ifstream fam(m_famFile);
    if (!fam.is_open()) {
        throw std::runtime_error("Cannot open PLINK .fam file: " + m_famFile);
    }

    m_N0 = 0;
    std::string line;

    while (getline(fam, line)) {
        m_N0++;
        std::vector<std::string> line_elements = splitLine(line);

        if (line_elements.size() < 2) {
            throw std::runtime_error(
                "PLINK .fam file has fewer than 2 columns at line " +
                std::to_string(m_N0));
        }

        m_SampleInPlink.push_back(line_elements[1]);  // put IID to m_SampleInPlink
    }

    m_N = m_N0;
    m_numBytesofEachMarker0 = (m_N0 + 3) / 4;
    m_OneMarkerG4.reserve(m_numBytesofEachMarker0);
    m_OneMarkerG4.resize(m_numBytesofEachMarker0);

    std::cout << "Number of samples in fam file: " << m_N0 << std::endl;
}

// ============================================================
// setPosSampleInPlink: find positions of model samples in PLINK file
//
// Replaces Rcpp::match() with std::unordered_map lookup
// ============================================================
void PlinkClass::setPosSampleInPlink(std::vector<std::string>& t_SampleInModel)
{
    std::cout << "Setting position of samples in PLINK files...." << std::endl;

    // If no sample IDs provided, use all fam samples in order (identity mapping)
    if (t_SampleInModel.empty()) {
        std::cout << "  No sampleIDs in null model; using all " << m_N0
                  << " samples from fam file." << std::endl;
        m_N = m_N0;
        m_numBytesofEachMarker = (m_N + 3) / 4;
        m_posSampleInPlink.resize(m_N);
        for (uint32_t i = 0; i < m_N; i++) {
            m_posSampleInPlink[i] = i;
        }
        std::cout << "Number of samples in analysis: " << m_N << std::endl;
        return;
    }

    m_N = t_SampleInModel.size();
    m_numBytesofEachMarker = (m_N + 3) / 4;

    // Build map of plink sample IID -> index (0-based)
    std::unordered_map<std::string, uint32_t> plinkSampleMap;
    for (uint32_t i = 0; i < m_SampleInPlink.size(); i++) {
        plinkSampleMap[m_SampleInPlink[i]] = i;
    }

    m_posSampleInPlink.resize(m_N);
    for (uint32_t i = 0; i < m_N; i++) {
        auto it = plinkSampleMap.find(t_SampleInModel[i]);
        if (it == plinkSampleMap.end()) {
            throw std::runtime_error(
                "At least one subject requested is not in Plink file. "
                "Sample not found: " + t_SampleInModel[i]);
        }
        m_posSampleInPlink[i] = it->second;  // already 0-based
    }

    std::cout << "Number of samples in analysis: " << m_N << std::endl;
}

// ============================================================
// closegenofile: close the .bed file handle
// ============================================================
void PlinkClass::closegenofile()
{
    if (m_fin) {
        fclose(m_fin);
        m_fin = nullptr;
    }
}

// ============================================================
// getOneMarker: read one marker from .bed and decode genotypes
//
// This is a direct port of SAIGE/src/PLINK.cpp::getOneMarker
// with the EXACT same algorithm:
//   - Seek in .bed file based on previous/current index
//   - Read m_numBytesofEachMarker0 bytes
//   - Decode 2-bit genotypes
//   - Apply allele order mapping
//   - Calculate alt freq, alt counts, missing rate
// ============================================================
void PlinkClass::getOneMarker(uint64_t& t_gIndex_prev,
                              uint64_t& t_gIndex,
                              std::string& t_ref,
                              std::string& t_alt,
                              std::string& t_marker,
                              uint32_t& t_pd,
                              std::string& t_chr,
                              double& t_altFreq,
                              double& t_altCounts,
                              double& t_missingRate,
                              double& t_imputeInfo,
                              bool& t_isOutputIndexForMissing,
                              std::vector<uint>& t_indexForMissing,
                              bool& t_isOnlyOutputNonZero,
                              std::vector<uint>& t_indexForNonZero,
                              bool& t_isTrueGenotype,
                              arma::vec& OneMarkerG1)
{
    // t_isTrueGenotype = FALSE is used only when calculating GRM
    if (!t_isTrueGenotype) {
        if (t_isOutputIndexForMissing) {
            throw std::runtime_error(
                "Check PlinkClass::getOneMarker, if t_isTrueGenotype = FALSE, "
                "then t_isOutputIndexForMissing should be FALSE.");
        }
        if (t_isOnlyOutputNonZero) {
            throw std::runtime_error(
                "Check PlinkClass::getOneMarker, if t_isTrueGenotype = FALSE, "
                "then t_isOnlyOutputNonZero should be FALSE.");
        }
    }

    // Seek to the correct position in the .bed file
    uint64_t posSeek;
    if (t_gIndex > 0) {
        if (t_gIndex_prev == 0) {  // if it is the first element
            posSeek = 3 + m_numBytesofEachMarker0 * t_gIndex;
            fseek(m_fin, posSeek, SEEK_SET);
        } else {
            posSeek = m_numBytesofEachMarker0 * (t_gIndex - t_gIndex_prev - 1);
            if (posSeek > 0) {
                fseek(m_fin, posSeek, SEEK_CUR);
            }
        }
    }

    // Read the raw genotype bytes for this marker
    fread((char*)(&m_OneMarkerG4[0]), 1, m_numBytesofEachMarker0, m_fin);

    t_indexForMissing.clear();
    t_indexForNonZero.clear();

    t_marker = m_MarkerInPlink[t_gIndex];
    t_pd = m_pd[t_gIndex];
    t_chr = m_chr[t_gIndex];

    // Select genotype mapping based on allele order
    std::vector<int8_t> genoMaps;

    if (m_AlleleOrder == "alt-first") {
        t_ref = m_ref[t_gIndex];
        t_alt = m_alt[t_gIndex];
        genoMaps = m_genoMaps_alt_first;
    }

    if (m_AlleleOrder == "ref-first") {
        t_ref = m_alt[t_gIndex];
        t_alt = m_ref[t_gIndex];
        genoMaps = m_genoMaps_ref_first;
    }

    // Decode genotypes for each sample in the analysis
    uint j = 0;
    int counts[] = {0, 0, 0, 0};

    for (uint32_t i = 0; i < m_N; i++) {
        auto ind = m_posSampleInPlink[i];             // C++ start from 0
        unsigned char bufferG4 = m_OneMarkerG4[ind / 4];  // 1 byte for 4 genotypes
        size_t bufferG1;                                   // 1 genotype (1 sample)

        // https://www.cog-genomics.org/plink/1.9/formats#bed
        getGenotype(bufferG4, ind % 4, bufferG1);     // bufferG4 -> bufferG1

        counts[bufferG1]++;

        if (bufferG1 == MISSING && t_isOutputIndexForMissing) {
            t_indexForMissing.push_back(i);
        }

        if (t_isTrueGenotype) {
            bufferG1 = genoMaps[bufferG1];
        }

        if (bufferG1 > 0) {
            t_indexForNonZero.push_back(i);
        }

        if (t_isOnlyOutputNonZero) {
            if (bufferG1 > 0) {
                OneMarkerG1[j] = bufferG1;
                j = j + 1;
            }
        } else {
            OneMarkerG1[i] = bufferG1;
        }
    }

    int numMissing = counts[MISSING];

    int count = m_N - numMissing;
    t_missingRate = (double)numMissing / (double)m_N;
    t_imputeInfo = 1;
    t_altCounts = (double)(counts[HET] + 2 * counts[HOM_ALT]);

    if (count > 0) {
        t_altFreq = t_altCounts / (double)count / 2;
    } else {
        t_altFreq = 0;
    }

    // updated on 03/14/2021
    if (m_AlleleOrder == "ref-first") {
        t_altFreq = 1 - t_altFreq;
        t_altCounts = 2 * (double)count * t_altFreq;
    }

    if (t_isOnlyOutputNonZero) {
        OneMarkerG1.resize(j);
    }

    // === CHECKPOINT OUTPUT BEGIN ===
    if (g_writeCheckpoints && s_checkpointMarkerCount == 0 && !g_checkpointDir.empty()) {
        // Create checkpoint directory if it does not exist
        std::filesystem::create_directories(g_checkpointDir);

        // Write marker info checkpoint
        std::string infoFile = g_checkpointDir + "/ckpt_07_marker0_raw.txt";
        std::ofstream ofs_info(infoFile);
        if (ofs_info.is_open()) {
            ofs_info << std::setprecision(15);
            ofs_info << "field\tvalue" << std::endl;
            ofs_info << "chr\t" << t_chr << std::endl;
            ofs_info << "pos\t" << t_pd << std::endl;
            ofs_info << "ref\t" << t_ref << std::endl;
            ofs_info << "alt\t" << t_alt << std::endl;
            ofs_info << "marker\t" << t_marker << std::endl;
            ofs_info << "altFreq\t" << t_altFreq << std::endl;
            ofs_info << "altCounts\t" << t_altCounts << std::endl;
            ofs_info << "missingRate\t" << t_missingRate << std::endl;
            ofs_info << "imputeInfo\t" << t_imputeInfo << std::endl;
            ofs_info.close();
            std::cout << "[CHECKPOINT] Wrote ckpt_07_marker0_raw.txt to " << g_checkpointDir << std::endl;
        }

        // Write first 20 genotype values checkpoint
        std::string genoFile = g_checkpointDir + "/ckpt_08_marker0_geno.txt";
        std::ofstream ofs_geno(genoFile);
        if (ofs_geno.is_open()) {
            ofs_geno << std::setprecision(15);
            ofs_geno << "index\tgenotype" << std::endl;
            uint32_t nToWrite = std::min((uint32_t)20, (uint32_t)OneMarkerG1.n_elem);
            for (uint32_t k = 0; k < nToWrite; k++) {
                ofs_geno << (k + 1) << "\t" << OneMarkerG1[k] << std::endl;
            }
            ofs_geno.close();
            std::cout << "[CHECKPOINT] Wrote ckpt_08_marker0_geno.txt to " << g_checkpointDir << std::endl;
        }
    }
    s_checkpointMarkerCount++;
    // === CHECKPOINT OUTPUT END ===
}

} // namespace PLINK


// ============================================================
// whichCPP: C++ version of which(). Note: start from 0, not 1
// (Outside PLINK namespace, as declared in header)
// ============================================================
std::vector<unsigned int> whichCPP(std::vector<std::string>& strVec,
                                   std::string strValue)
{
    std::vector<unsigned int> indexVec;
    for (unsigned int i = 0; i < strVec.size(); i++) {
        if (strVec.at(i) == strValue)
            indexVec.push_back(i);
    }
    return indexVec;
}


// ============================================================
// setPLINKobjInCPP: create and configure the global PLINK object
// Mirrors SAIGE Main.cpp::setPLINKobjInCPP
// ============================================================
void setPLINKobjInCPP(std::string t_bimFile,
                      std::string t_famFile,
                      std::string t_bedFile,
                      std::vector<std::string>& t_SampleInModel,
                      std::string t_AlleleOrder)
{
    // Clean up any existing object
    if (ptr_gPLINKobj != nullptr) {
        ptr_gPLINKobj->closegenofile();
        delete ptr_gPLINKobj;
        ptr_gPLINKobj = nullptr;
    }

    ptr_gPLINKobj = new PLINK::PlinkClass(t_bimFile,
                                          t_famFile,
                                          t_bedFile,
                                          t_AlleleOrder);

    ptr_gPLINKobj->setPosSampleInPlink(t_SampleInModel);

    int n = ptr_gPLINKobj->getN();
    std::cout << "n:\t" << n << std::endl;
}


// ============================================================
// Unified_getOneMarker: unified dispatcher for genotype reading
//
// Ported from SAIGE Main.cpp. Currently only supports "plink".
// Throws for "bgen", "vcf", "pgen" (not yet implemented).
// ============================================================
bool Unified_getOneMarker(std::string& t_genoType,
                          uint64_t& t_gIndex_prev,
                          uint64_t& t_gIndex,
                          std::string& t_ref,
                          std::string& t_alt,
                          std::string& t_marker,
                          uint32_t& t_pd,
                          std::string& t_chr,
                          double& t_altFreq,
                          double& t_altCounts,
                          double& t_missingRate,
                          double& t_imputeInfo,
                          bool& t_isOutputIndexForMissing,
                          std::vector<uint>& t_indexForMissing,
                          bool& t_isOnlyOutputNonZero,
                          std::vector<uint>& t_indexForNonZero,
                          arma::vec& t_GVec,
                          bool t_isImputation)
{
    bool isBoolRead = true;

    if (t_genoType == "plink") {
        bool isTrueGenotype = true;
        // t_gIndex_prev is after reading the last marker

        if (ptr_gPLINKobj == nullptr) {
            throw std::runtime_error(
                "Unified_getOneMarker: PLINK object not initialized. "
                "Call setPLINKobjInCPP first.");
        }

        ptr_gPLINKobj->getOneMarker(t_gIndex_prev, t_gIndex,
                                    t_ref, t_alt, t_marker, t_pd, t_chr,
                                    t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
                                    t_isOutputIndexForMissing, t_indexForMissing,
                                    t_isOnlyOutputNonZero, t_indexForNonZero,
                                    isTrueGenotype, t_GVec);
    } else if (t_genoType == "bgen") {
        throw std::runtime_error(
            "Unified_getOneMarker: BGEN format not yet supported in standalone C++ port.");
    } else if (t_genoType == "vcf") {
        throw std::runtime_error(
            "Unified_getOneMarker: VCF format not yet supported in standalone C++ port.");
    } else if (t_genoType == "pgen") {
        throw std::runtime_error(
            "Unified_getOneMarker: PGEN format not yet supported in standalone C++ port.");
    } else {
        throw std::runtime_error(
            "Unified_getOneMarker: Unknown genotype type '" + t_genoType + "'. "
            "Supported types: plink");
    }

    return isBoolRead;
}


// ============================================================
// Unified_getSampleSizeinGeno: get total sample size in genotype file
// ============================================================
uint32_t Unified_getSampleSizeinGeno(std::string& t_genoType)
{
    uint32_t N0;
    if (t_genoType == "plink") {
        if (ptr_gPLINKobj == nullptr) {
            throw std::runtime_error(
                "Unified_getSampleSizeinGeno: PLINK object not initialized.");
        }
        N0 = ptr_gPLINKobj->getN0();
    } else {
        throw std::runtime_error(
            "Unified_getSampleSizeinGeno: Only 'plink' type supported. Got: " + t_genoType);
    }
    return N0;
}


// ============================================================
// Unified_getSampleSizeinAnalysis: get sample size in analysis (after subsetting)
// ============================================================
uint32_t Unified_getSampleSizeinAnalysis(std::string& t_genoType)
{
    uint32_t N;
    if (t_genoType == "plink") {
        if (ptr_gPLINKobj == nullptr) {
            throw std::runtime_error(
                "Unified_getSampleSizeinAnalysis: PLINK object not initialized.");
        }
        N = ptr_gPLINKobj->getN();
    } else {
        throw std::runtime_error(
            "Unified_getSampleSizeinAnalysis: Only 'plink' type supported. Got: " + t_genoType);
    }
    return N;
}


// ============================================================
// closeGenoFile: close the genotype file
// ============================================================
void closeGenoFile(std::string& t_genoType)
{
    if (t_genoType == "plink") {
        if (ptr_gPLINKobj != nullptr) {
            ptr_gPLINKobj->closegenofile();
        }
    } else {
        throw std::runtime_error(
            "closeGenoFile: Only 'plink' type supported. Got: " + t_genoType);
    }
}
