// Standalone port of SAIGE/src/PLINK.cpp, SAIGE/src/VCF.cpp, and SAIGE/src/BGEN.cpp
// PLINK + VCF + BGEN genotype readers -- all Rcpp dependencies removed
//
// Conversions from original SAIGE:
//   1. #include <RcppArmadillo.h>  -->  #include <armadillo>
//   2. Rcpp::stop(...)             -->  throw std::runtime_error(...)
//   3. boost::split(...)           -->  splitLine() using std::istringstream
//   4. boost::replace_all(...)     -->  manual \r removal
//   5. Rcpp::match(...)            -->  std::unordered_map lookup
//   6. Rcpp::CharacterVector       -->  std::vector<std::string>
//   7. Rcpp::IntegerVector         -->  std::vector<uint32_t>
//   8. savvy library (VCF)         -->  htslib (VCF/BCF/VCF.GZ)
//   9. Rcpp BGEN                   -->  standalone BGEN with zstd/zlib

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
#include <cmath>
#include <sys/stat.h>
#include <filesystem>
#include <iomanip>

// htslib headers for VCF reading
#include <htslib/hts.h>
#include <htslib/vcf.h>

// zstd and zlib headers for BGEN decompression
#include <zstd.h>
#include <zlib.h>

// Global checkpoint flag for debugging
bool g_writeCheckpoints = false;
std::string g_checkpointDir = "";

// Global pointer to PLINK object (mirrors SAIGE's ptr_gPLINKobj in Main.cpp)
PLINK::PlinkClass* ptr_gPLINKobj = nullptr;

// Global pointer to VCF object (mirrors SAIGE's ptr_gVCFobj in Main.cpp)
VCF::VcfClass* ptr_gVCFobj = nullptr;

// Global pointer to BGEN object (mirrors SAIGE's ptr_gBGENobj in Main.cpp)
BGEN::BgenClass* ptr_gBGENobj = nullptr;

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
// VCF namespace: VCF/BCF/VCF.GZ reader using htslib
// Ported from SAIGE/src/VCF.cpp with savvy -> htslib conversion
// ============================================================
namespace VCF {

// ============================================================
// Constructor
// ============================================================
VcfClass::VcfClass(const std::string& t_vcfFileName,
                   const std::string& t_vcfField,
                   std::vector<std::string>& t_SampleInModel)
    : m_htsFile(nullptr), m_hdr(nullptr), m_rec(nullptr),
      m_vcfFileName(t_vcfFileName), m_fmtField(t_vcfField),
      m_N0(0), m_N(0), m_M0(0), m_totalMarkers(0), m_isPreScanned(false)
{
    // Open VCF/BCF file
    m_htsFile = hts_open(t_vcfFileName.c_str(), "r");
    if (!m_htsFile) {
        throw std::runtime_error("Cannot open VCF file: " + t_vcfFileName);
    }

    // Read header
    m_hdr = bcf_hdr_read(m_htsFile);
    if (!m_hdr) {
        hts_close(m_htsFile);
        m_htsFile = nullptr;
        throw std::runtime_error("Cannot read VCF header from: " + t_vcfFileName);
    }

    // Allocate record
    m_rec = bcf_init();

    // Get sample IDs
    getSampleIDlist();
    std::cout << "Open VCF done" << std::endl;
    std::cout << "To read the field " << t_vcfField << std::endl;
    std::cout << "Number of samples in the vcf file: " << m_N0 << std::endl;

    // Verify the format field exists
    // Check if the requested format field is in the header
    int fmt_id = bcf_hdr_id2int(m_hdr, BCF_DT_ID, t_vcfField.c_str());
    if (fmt_id < 0 || !bcf_hdr_idinfo_exists(m_hdr, BCF_HL_FMT, fmt_id)) {
        // Try fallback to HDS if DS or GT was requested
        if ((t_vcfField == "DS" || t_vcfField == "GT")) {
            int hds_id = bcf_hdr_id2int(m_hdr, BCF_DT_ID, "HDS");
            if (hds_id >= 0 && bcf_hdr_idinfo_exists(m_hdr, BCF_HL_FMT, hds_id)) {
                m_fmtField = "HDS";
                std::cout << "Format field " << t_vcfField << " not found, using HDS instead" << std::endl;
            } else {
                std::cerr << "WARNING: vcfField (" << t_vcfField << ") not present in genotype file." << std::endl;
            }
        } else {
            std::cerr << "WARNING: vcfField (" << t_vcfField << ") not present in genotype file." << std::endl;
        }
    }

    // Set up sample mapping
    setPosSampleInVcf(t_SampleInModel);
}

// ============================================================
// Destructor
// ============================================================
VcfClass::~VcfClass()
{
    closegenofile();
}

// ============================================================
// getSampleIDlist: extract sample IDs from VCF header
// ============================================================
void VcfClass::getSampleIDlist()
{
    m_N0 = bcf_hdr_nsamples(m_hdr);
    m_SampleInVcf.resize(m_N0);
    for (uint32_t i = 0; i < m_N0; i++) {
        m_SampleInVcf[i] = m_hdr->samples[i];
    }
}

// ============================================================
// setPosSampleInVcf: find positions of model samples in VCF file
// Mirrors SAIGE VCF.cpp::setPosSampleInVcf but replaces Rcpp::match
// with std::unordered_map lookup
// ============================================================
void VcfClass::setPosSampleInVcf(std::vector<std::string>& t_SampleInModel)
{
    std::cout << "Setting position of samples in VCF files...." << std::endl;

    // If no sample IDs provided, use all VCF samples in order
    if (t_SampleInModel.empty()) {
        m_N = m_N0;
        std::cout << "  No sampleIDs in null model; using all " << m_N0
                  << " samples from VCF file." << std::endl;

        m_posSampleInModel.resize(m_N0);
        for (uint32_t i = 0; i < m_N0; i++) {
            m_posSampleInModel[i] = i;  // identity mapping
        }
    } else {
        m_N = t_SampleInModel.size();
        std::cout << "m_N " << m_N << std::endl;

        // Build map of model sample ID -> model index (0-based)
        std::unordered_map<std::string, int32_t> modelSampleMap;
        for (uint32_t i = 0; i < m_N; i++) {
            modelSampleMap[t_SampleInModel[i]] = (int32_t)i;
        }

        // Verify all model samples exist in VCF
        for (uint32_t i = 0; i < m_N; i++) {
            bool found = false;
            for (uint32_t j = 0; j < m_N0; j++) {
                if (m_SampleInVcf[j] == t_SampleInModel[i]) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                throw std::runtime_error(
                    "At least one subject requested is not in VCF file. "
                    "Sample not found: " + t_SampleInModel[i]);
            }
        }

        // Build mapping: for each VCF sample, what is its position in the model?
        // m_posSampleInModel[vcf_idx] = model_idx, or -1 if not in model
        m_posSampleInModel.resize(m_N0);
        for (uint32_t i = 0; i < m_N0; i++) {
            auto it = modelSampleMap.find(m_SampleInVcf[i]);
            if (it != modelSampleMap.end()) {
                m_posSampleInModel[i] = it->second;
            } else {
                m_posSampleInModel[i] = -1;
            }
        }
    }

    std::cout << "Number of samples in analysis: " << m_N << std::endl;
}

// ============================================================
// getOneMarker: read one marker from VCF and extract dosages
//
// This mirrors SAIGE VCF.cpp::getOneMarker but uses htslib
// instead of savvy. Reads sequentially (no random access).
//
// Supports GT (hard calls -> 0/1/2) and DS (dosage) fields.
// ============================================================
bool VcfClass::getOneMarker(
    std::string& t_ref,
    std::string& t_alt,
    std::string& t_marker,
    uint32_t& t_pd,
    std::string& t_chr,
    double& t_altFreq,
    double& t_altCounts,
    double& t_missingRate,
    double& t_imputeInfo,
    bool t_isOutputIndexForMissing,
    std::vector<uint>& t_indexForMissing,
    bool t_isOnlyOutputNonZero,
    std::vector<uint>& t_indexForNonZero,
    arma::vec& dosages,
    bool t_isImputation)
{
    t_indexForMissing.clear();
    t_indexForNonZero.clear();

    // Read next record
    int ret = bcf_read(m_htsFile, m_hdr, m_rec);
    if (ret < 0) {
        // End of file or error
        std::cout << "Reach the end of the vcf file" << std::endl;
        return false;
    }

    // Unpack the record (we need INFO, FORMAT, and shared fields)
    bcf_unpack(m_rec, BCF_UN_ALL);

    // Extract chromosome
    t_chr = bcf_hdr_id2name(m_hdr, m_rec->rid);

    // Extract position (VCF is 1-based, htslib stores 0-based)
    t_pd = (uint32_t)(m_rec->pos + 1);

    // Extract REF
    t_ref = m_rec->d.allele[0];

    // Extract ALT (first ALT allele only, skip multiallelic)
    if (m_rec->n_allele < 2) {
        t_alt = ".";
    } else {
        if (m_rec->n_allele > 2) {
            std::cerr << "Warning: skipping multiallelic variant at "
                      << t_chr << ":" << t_pd << std::endl;
        }
        t_alt = m_rec->d.allele[1];
    }

    // Extract marker ID
    if (m_rec->d.id && std::string(m_rec->d.id) != ".") {
        t_marker = m_rec->d.id;
    } else {
        t_marker = t_chr + ":" + std::to_string(t_pd) + ":" + t_ref + ":" + t_alt;
    }

    // Store marker info for later lookup
    m_chr.push_back(t_chr);
    m_pd.push_back(t_pd);
    m_ref.push_back(t_ref);
    m_alt.push_back(t_alt);
    m_MarkerInVcf.push_back(t_marker);
    m_M0++;

    // Initialize dosage vector
    dosages.set_size(m_N);
    dosages.fill(0.0);

    t_altCounts = 0;
    int missing_cnt = 0;

    // Extract imputation info (R2) if requested
    if (t_isImputation) {
        float* info_val = nullptr;
        int n_info = 0;
        int info_ret = bcf_get_info_float(m_hdr, m_rec, "R2", &info_val, &n_info);
        if (info_ret > 0 && n_info > 0) {
            t_imputeInfo = (double)info_val[0];
        } else {
            t_imputeInfo = 1.0;
        }
        if (info_val) free(info_val);
    } else {
        t_imputeInfo = 1.0;
    }

    // Read genotype data based on format field
    if (m_fmtField == "GT") {
        // Read GT field (hard-call genotypes)
        int32_t* gt_arr = nullptr;
        int n_gt = 0;
        int ngt_ret = bcf_get_genotypes(m_hdr, m_rec, &gt_arr, &n_gt);

        if (ngt_ret <= 0) {
            // No GT data available
            dosages.fill(-1.0);
            missing_cnt = m_N;
        } else {
            int ploidy = n_gt / m_N0;

            for (uint32_t i = 0; i < m_N0; i++) {
                int32_t model_idx = m_posSampleInModel[i];
                if (model_idx < 0) continue;  // sample not in model

                int dose = 0;
                bool is_missing = false;

                for (int p = 0; p < ploidy; p++) {
                    int32_t allele = gt_arr[i * ploidy + p];

                    if (allele == bcf_int32_vector_end) {
                        break;  // fewer ploidy for this sample
                    }
                    if (bcf_gt_is_missing(allele)) {
                        is_missing = true;
                        break;
                    }
                    // bcf_gt_allele extracts the allele index (0 = REF, 1 = ALT)
                    dose += bcf_gt_allele(allele);
                }

                if (is_missing) {
                    dosages[model_idx] = -1.0;
                    missing_cnt++;
                    if (t_isOutputIndexForMissing) {
                        t_indexForMissing.push_back((uint)model_idx);
                    }
                } else {
                    dosages[model_idx] = (double)dose;
                    t_altCounts += dose;
                    if (dose > 0) {
                        t_indexForNonZero.push_back((uint)model_idx);
                    }
                }
            }
        }

        if (gt_arr) free(gt_arr);

    } else {
        // Read DS or HDS field (dosage)
        float* ds_arr = nullptr;
        int n_ds = 0;
        int nds_ret = bcf_get_format_float(m_hdr, m_rec, m_fmtField.c_str(), &ds_arr, &n_ds);

        if (nds_ret <= 0) {
            // No dosage data available
            dosages.fill(-1.0);
            missing_cnt = m_N;
        } else {
            int stride = n_ds / m_N0;

            for (uint32_t i = 0; i < m_N0; i++) {
                int32_t model_idx = m_posSampleInModel[i];
                if (model_idx < 0) continue;  // sample not in model

                float dose_val;
                if (stride == 1) {
                    // DS field: single dosage value per sample
                    dose_val = ds_arr[i];
                } else {
                    // HDS field: sum haplotype dosages (like savvy stride_reduce)
                    dose_val = 0.0f;
                    for (int p = 0; p < stride; p++) {
                        float v = ds_arr[i * stride + p];
                        if (!bcf_float_is_missing(v) && !bcf_float_is_vector_end(v)) {
                            dose_val += v;
                        }
                    }
                }

                if (bcf_float_is_missing(dose_val) || std::isnan(dose_val)) {
                    dosages[model_idx] = -1.0;
                    missing_cnt++;
                    if (t_isOutputIndexForMissing) {
                        t_indexForMissing.push_back((uint)model_idx);
                    }
                } else {
                    dosages[model_idx] = (double)dose_val;
                    t_altCounts += dose_val;
                    if (dose_val > 0) {
                        t_indexForNonZero.push_back((uint)model_idx);
                    }
                }
            }
        }

        if (ds_arr) free(ds_arr);
    }

    // Calculate alt frequency and missing rate (same logic as SAIGE VCF.cpp)
    if (missing_cnt > 0) {
        if (missing_cnt == (int)m_N) {
            t_altFreq = 0;
        } else {
            t_altFreq = t_altCounts / 2.0 / (double)(m_N - missing_cnt);
        }
        t_missingRate = (double)missing_cnt / (double)m_N;
    } else {
        t_altFreq = t_altCounts / 2.0 / (double)m_N;
        t_missingRate = 0;
    }

    // Convert -1 (missing) values in dosages to the proper SAIGE convention
    // SAIGE VCF.cpp stores -1 for missing, but the downstream pipeline uses
    // the same missing convention as PLINK: impute or exclude.
    // We leave -1 in the dosage vector; the imputation step will handle it.
    // Actually, let's match PLINK convention: missing genotypes get value -1
    // (which the downstream imputeGenoAndFlip handles via indexForMissing).

    return true;
}

// ============================================================
// prescanMarkerCount: count total markers in VCF by scanning
// ============================================================
uint32_t VcfClass::prescanMarkerCount()
{
    if (m_isPreScanned) {
        return m_totalMarkers;
    }

    // Save current position and re-open file for counting
    htsFile* countFile = hts_open(m_vcfFileName.c_str(), "r");
    if (!countFile) {
        throw std::runtime_error("Cannot re-open VCF file for marker counting: " + m_vcfFileName);
    }
    bcf_hdr_t* countHdr = bcf_hdr_read(countFile);
    bcf1_t* countRec = bcf_init();

    m_totalMarkers = 0;
    while (bcf_read(countFile, countHdr, countRec) >= 0) {
        m_totalMarkers++;
    }

    bcf_destroy(countRec);
    bcf_hdr_destroy(countHdr);
    hts_close(countFile);

    m_isPreScanned = true;
    std::cout << "Number of markers in VCF file: " << m_totalMarkers << std::endl;
    return m_totalMarkers;
}

// ============================================================
// resetFile: close and re-open the VCF file to restart reading
// ============================================================
void VcfClass::resetFile()
{
    // Close current file
    if (m_rec) {
        bcf_destroy(m_rec);
        m_rec = nullptr;
    }
    if (m_hdr) {
        bcf_hdr_destroy(m_hdr);
        m_hdr = nullptr;
    }
    if (m_htsFile) {
        hts_close(m_htsFile);
        m_htsFile = nullptr;
    }

    // Clear accumulated marker info
    m_chr.clear();
    m_pd.clear();
    m_ref.clear();
    m_alt.clear();
    m_MarkerInVcf.clear();
    m_M0 = 0;

    // Re-open
    m_htsFile = hts_open(m_vcfFileName.c_str(), "r");
    if (!m_htsFile) {
        throw std::runtime_error("Cannot re-open VCF file: " + m_vcfFileName);
    }
    m_hdr = bcf_hdr_read(m_htsFile);
    if (!m_hdr) {
        throw std::runtime_error("Cannot re-read VCF header from: " + m_vcfFileName);
    }
    m_rec = bcf_init();
}

// ============================================================
// getMarkerIDToIndex: build chr:pos:ref:alt -> index map
// Requires markers to have been read (or pre-scanned)
// ============================================================
std::unordered_map<std::string, uint32_t> VcfClass::getMarkerIDToIndex()
{
    std::unordered_map<std::string, uint32_t> idMap;
    for (uint32_t i = 0; i < m_chr.size(); i++) {
        // Primary key: chr:pos:ref:alt
        std::string id1 = m_chr[i] + ":" + std::to_string(m_pd[i]) + ":"
                         + m_ref[i] + ":" + m_alt[i];
        idMap[id1] = i;
        // Secondary key: chr:pos:alt:ref (reversed allele order)
        std::string id2 = m_chr[i] + ":" + std::to_string(m_pd[i]) + ":"
                         + m_alt[i] + ":" + m_ref[i];
        if (idMap.find(id2) == idMap.end()) {
            idMap[id2] = i;
        }
    }
    return idMap;
}

// ============================================================
// getMarkerNameToIndex: build marker ID -> index map
// ============================================================
std::unordered_map<std::string, uint32_t> VcfClass::getMarkerNameToIndex()
{
    std::unordered_map<std::string, uint32_t> nameMap;
    for (uint32_t i = 0; i < m_MarkerInVcf.size(); i++) {
        nameMap[m_MarkerInVcf[i]] = i;
    }
    return nameMap;
}

// ============================================================
// closegenofile: clean up htslib resources
// ============================================================
void VcfClass::closegenofile()
{
    if (m_rec) {
        bcf_destroy(m_rec);
        m_rec = nullptr;
    }
    if (m_hdr) {
        bcf_hdr_destroy(m_hdr);
        m_hdr = nullptr;
    }
    if (m_htsFile) {
        hts_close(m_htsFile);
        m_htsFile = nullptr;
    }
}

} // namespace VCF


// ============================================================
// BGEN namespace: BGEN v1.2 reader (manual binary parsing)
// Ported from SAIGE/src/BGEN.cpp with all Rcpp dependencies removed
// ============================================================
namespace BGEN {

// ============================================================
// Constructor
// ============================================================
BgenClass::BgenClass(const std::string& t_bgenFileName,
                     const std::vector<std::string>& t_SampleInBgen,
                     std::vector<std::string>& t_SampleInModel,
                     const std::string& t_AlleleOrder)
    : m_fin(nullptr), m_zBufLens(0), m_bufLens(0),
      CompressedSNPBlocks(0), m_N0(0), m_N(0), m_M0(0)
{
    setBgenObj(t_bgenFileName, t_SampleInBgen);
    setPosSampleInBgen(t_SampleInModel);
    m_AlleleOrder = t_AlleleOrder;
    allele0.resize(65536);
    allele1.resize(65536);
}

// ============================================================
// Destructor
// ============================================================
BgenClass::~BgenClass()
{
    closegenofile();
}

// ============================================================
// setBgenObj: read BGEN header, validate format
// Ported from SAIGE/src/BGEN.cpp::setBgenObj
// ============================================================
void BgenClass::setBgenObj(const std::string& t_bgenFileName,
                            const std::vector<std::string>& t_SampleInBgen)
{
    m_bgenFileName = t_bgenFileName;
    m_SampleInBgen = t_SampleInBgen;

    m_fin = fopen(t_bgenFileName.c_str(), "rb");
    if (!m_fin) {
        throw std::runtime_error("Cannot open BGEN file: " + t_bgenFileName);
    }

    // Read BGEN v1.2 header
    uint32_t offset;
    fread(&offset, 4, 1, m_fin);

    uint32_t L_H;
    fread(&L_H, 4, 1, m_fin);

    fread(&m_M0, 4, 1, m_fin);
    std::cout << "snpBlocks (Mbgen): " << m_M0 << std::endl;
    if (m_M0 == 0) {
        throw std::runtime_error("No snpBlocks in BGEN file.");
    }

    fread(&m_N0, 4, 1, m_fin);
    std::cout << "samples (Nbgen): " << m_N0 << std::endl;

    unsigned int m_Nsample = t_SampleInBgen.size();
    if (m_N0 != m_Nsample) {
        throw std::runtime_error(
            "Number of samples in BGEN header (" + std::to_string(m_N0) +
            ") does not match sample file (" + std::to_string(m_Nsample) + ")");
    }

    char magic[5];
    fread(magic, 1, 4, m_fin);
    magic[4] = '\0';

    // Skip free data area
    fseek(m_fin, L_H - 20, SEEK_CUR);

    uint32_t flags;
    fread(&flags, 4, 1, m_fin);
    CompressedSNPBlocks = flags & 3;
    std::cout << "CompressedSNPBlocks: " << CompressedSNPBlocks << std::endl;

    if (CompressedSNPBlocks == COMPRESSION_ZLIB) {
        std::cout << "Warning: your bgen file uses zlib compression, which is slow." << std::endl;
    } else if (CompressedSNPBlocks != COMPRESSION_ZSTD) {
        throw std::runtime_error("Bgenreader only supports zlib or zstd compression.");
    }

    uint32_t Layout = (flags >> 2) & 0xf;
    std::cout << "Layout: " << Layout << std::endl;
    if (Layout != 1 && Layout != 2) {
        throw std::runtime_error("Bgenreader only supports Layout 1 or 2.");
    }

    // Seek past the header to the data area
    fseek(m_fin, offset + 4, SEEK_SET);
}

// ============================================================
// setPosSampleInBgen: find positions of model samples in BGEN
// Replaces Rcpp::match() with std::unordered_map lookup
// ============================================================
void BgenClass::setPosSampleInBgen(std::vector<std::string>& t_SampleInModel)
{
    std::cout << "Setting position of samples in Bgen files...." << std::endl;

    // If no sample IDs provided, use all BGEN samples in order
    if (t_SampleInModel.empty()) {
        m_N = m_N0;
        std::cout << "  No sampleIDs in null model; using all " << m_N0
                  << " samples from BGEN file." << std::endl;
        m_posSampleInModel.resize(m_N0);
        for (uint32_t i = 0; i < m_N0; i++) {
            m_posSampleInModel[i] = i;
        }
        std::cout << "Number of samples in analysis: " << m_N << std::endl;
        return;
    }

    m_N = t_SampleInModel.size();

    // Build map of BGEN sample -> index (0-based)
    std::unordered_map<std::string, uint32_t> bgenSampleMap;
    for (uint32_t i = 0; i < m_N0; i++) {
        bgenSampleMap[m_SampleInBgen[i]] = i;
    }

    // Verify all model samples exist in BGEN
    for (uint32_t i = 0; i < m_N; i++) {
        if (bgenSampleMap.find(t_SampleInModel[i]) == bgenSampleMap.end()) {
            throw std::runtime_error(
                "At least one subject requested is not in BGEN file. "
                "Sample not found: " + t_SampleInModel[i]);
        }
    }

    // Build map of model sample -> model index (0-based)
    std::unordered_map<std::string, int32_t> modelSampleMap;
    for (uint32_t i = 0; i < m_N; i++) {
        modelSampleMap[t_SampleInModel[i]] = (int32_t)i;
    }

    // For each BGEN sample, what is its position in the model?
    // m_posSampleInModel[bgen_idx] = model_idx, or -1 if not in model
    m_posSampleInModel.resize(m_N0);
    for (uint32_t i = 0; i < m_N0; i++) {
        auto it = modelSampleMap.find(m_SampleInBgen[i]);
        if (it != modelSampleMap.end()) {
            m_posSampleInModel[i] = it->second;
        } else {
            m_posSampleInModel[i] = -1;
        }
    }

    std::cout << "Number of samples in analysis: " << m_N << std::endl;
}

// ============================================================
// Parse2: decompress and parse genotype probabilities
// Direct port from SAIGE/src/BGEN.cpp::Parse2
// ============================================================
void BgenClass::Parse2(unsigned char* buf, uint32_t bufLen,
                        const unsigned char* zBuf, uint32_t zBufLen,
                        std::string& snpName,
                        arma::vec& dosages,
                        double& AC, double& AF,
                        std::vector<uint>& indexforMissing,
                        double& info,
                        std::vector<uint>& indexNonZero,
                        bool isImputation)
{
    // Decompress
    if (CompressedSNPBlocks == COMPRESSION_ZLIB) {
        z_stream strm = {};
        strm.next_in = const_cast<Bytef*>(zBuf);
        strm.avail_in = zBufLen;
        strm.next_out = buf;
        strm.avail_out = bufLen;

        if (inflateInit(&strm) != Z_OK) {
            std::cerr << "inflateInit failed" << std::endl;
            return;
        }

        int ret = inflate(&strm, Z_FINISH);
        if (ret != Z_STREAM_END) {
            std::cerr << "inflate failed with code " << ret << std::endl;
        }

        inflateEnd(&strm);
    }
    if (CompressedSNPBlocks == COMPRESSION_ZSTD) {
        ZSTD_DCtx* dctx = ZSTD_createDCtx();
        size_t actual_decompressed_size = ZSTD_decompressDCtx(dctx, buf, bufLen, zBuf, zBufLen);
        if (ZSTD_isError(actual_decompressed_size)) {
            std::cerr << "Decompression failed: " << ZSTD_getErrorName(actual_decompressed_size) << std::endl;
            ZSTD_freeDCtx(dctx);
            return;
        }
        ZSTD_freeDCtx(dctx);
    }

    unsigned char* bufAt = buf;
    uint32_t N = bufAt[0] | (bufAt[1] << 8) | (bufAt[2] << 16) | (bufAt[3] << 24);
    bufAt += 4;

    if (N != m_N0) {
        std::cerr << "ERROR: " << snpName << " has N = " << N
                  << " (mismatch with header block)" << std::endl;
        throw std::runtime_error("BGEN sample count mismatch in variant " + snpName);
    }

    uint32_t K = bufAt[0] | (bufAt[1] << 8);
    bufAt += 2;
    if (K != 2U) {
        std::cerr << "ERROR: " << snpName << " has K = " << K
                  << " (non-bi-allelic)" << std::endl;
        throw std::runtime_error("Non-bi-allelic variant in BGEN: " + snpName);
    }

    uint32_t Pmin = *bufAt; bufAt++;
    if (Pmin != 2U) {
        std::cerr << "ERROR: " << snpName << " has minimum ploidy = " << Pmin
                  << " (not 2)" << std::endl;
        throw std::runtime_error("Unsupported minimum ploidy in BGEN: " + snpName);
    }

    uint32_t Pmax = *bufAt; bufAt++;
    if (Pmax != 2U) {
        std::cerr << "ERROR: " << snpName << " has maximum ploidy = " << Pmax
                  << " (not 2)" << std::endl;
        throw std::runtime_error("Unsupported maximum ploidy in BGEN: " + snpName);
    }

    // Read per-sample ploidy/missingness bytes
    const unsigned char* ploidyMissBytes = bufAt;
    for (uint32_t i = 0; i < N; i++) {
        uint32_t ploidyMiss = *bufAt; bufAt++;
        if (ploidyMiss != 2U && ploidyMiss != 130U) {
            std::cerr << "ERROR: " << snpName << " has ploidy/missingness byte = "
                      << ploidyMiss << " (not 2 or 130)" << std::endl;
            throw std::runtime_error("Unsupported ploidy/missingness in BGEN: " + snpName);
        }
    }

    uint32_t Phased = *bufAt; bufAt++;
    if (Phased != 0U) {
        std::cerr << "ERROR: " << snpName << " has Phased = " << Phased
                  << " (not 0)" << std::endl;
        throw std::runtime_error("Phased data not supported in BGEN reader: " + snpName);
    }

    uint32_t B = *bufAt; bufAt++;
    if (B != 8U) {
        std::cerr << "ERROR: " << snpName << " has B = " << B
                  << " (not 8)" << std::endl;
        throw std::runtime_error("Unsupported bit depth in BGEN: " + snpName);
    }

    // Build lookup table for 8-bit probabilities
    double lut[256];
    for (int i = 0; i <= 255; i++) {
        lut[i] = i / 255.0;
    }

    // Parse genotype probabilities
    double sum_eij = 0, sum_fij_minus_eij2 = 0, sum_eij_sub = 0;
    double p11, p10, dosage, eij, fij;
    double dosage_new;

    dosages.set_size(m_N);
    dosages.fill(arma::datum::nan);

    std::size_t missing_cnt = 0;

    for (uint32_t i = 0; i < N; i++) {
        if (ploidyMissBytes[i] != 130U) {
            // Not missing
            p11 = lut[*bufAt]; bufAt++;
            p10 = lut[*bufAt]; bufAt++;

            if (m_posSampleInModel[i] >= 0) {
                dosage = 2 * p11 + p10;

                // SAIGE default is ref-first: dosage_new = 2 - dosage
                dosage_new = 2 - dosage;

                eij = dosage;
                fij = 4 * p11 + p10;
                sum_eij += eij;
                sum_fij_minus_eij2 += fij - eij * eij;

                dosages[m_posSampleInModel[i]] = dosage_new;
                if (dosage_new > 0) {
                    indexNonZero.push_back(m_posSampleInModel[i]);
                }
                sum_eij_sub += eij;
            }
        } else {
            // Missing (ploidy = 130)
            bufAt += 2;
            if (m_posSampleInModel[i] >= 0) {
                indexforMissing.push_back(m_posSampleInModel[i]);
                ++missing_cnt;
                dosages[m_posSampleInModel[i]] = -1;
            }
        }
    }

    // Compute AC and AF
    AC = 2 * ((double)(m_N - missing_cnt)) - sum_eij_sub;
    if (m_N == missing_cnt) {
        AF = 0;
    } else {
        AF = AC / 2 / ((double)(m_N - missing_cnt));
    }

    // Compute imputation info
    double thetaHat = sum_eij / (2 * (m_N - missing_cnt));
    if (isImputation) {
        info = (thetaHat == 0 || thetaHat == 1) ? 1 :
            1 - sum_fij_minus_eij2 / (2 * (m_N - missing_cnt) * thetaHat * (1 - thetaHat));
    } else {
        info = 1.0;
    }
}

// ============================================================
// getOneMarker: read one marker from BGEN
// Direct port from SAIGE/src/BGEN.cpp::getOneMarker
// t_gIndex / t_gIndex_prev are byte positions (BGEN uses byte seeking)
// ============================================================
void BgenClass::getOneMarker(uint64_t& t_gIndex_prev,
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
                              bool& t_isBoolRead,
                              arma::vec& dosages,
                              bool t_isImputation)
{
    // Seek to correct position (BGEN uses byte offsets)
    if (t_gIndex > 0) {
        if (t_gIndex_prev > 0) {
            uint64_t posSeek = t_gIndex - t_gIndex_prev;
            if (posSeek > 0) {
                fseek(m_fin, posSeek, SEEK_CUR);
            }
        } else {
            fseek(m_fin, t_gIndex, SEEK_SET);
        }
    }

    std::string SNPID, RSID, chromosome, first_allele, second_allele;
    uint32_t position;
    double AC, AF, info;

    t_indexForMissing.clear();
    t_indexForNonZero.clear();

    char snpID[65536], rsID[65536], chrStr[65536];
    uint16_t LS;
    size_t numBoolRead = fread(&LS, 2, 1, m_fin);

    if (numBoolRead > 0) {
        t_isBoolRead = true;

        fread(snpID, 1, LS, m_fin); snpID[LS] = '\0';

        uint16_t LR;
        fread(&LR, 2, 1, m_fin);
        fread(rsID, 1, LR, m_fin); rsID[LR] = '\0';
        RSID = std::string(rsID) == "." ? snpID : rsID;

        uint16_t LC;
        fread(&LC, 2, 1, m_fin);
        fread(chrStr, 1, LC, m_fin); chrStr[LC] = '\0';
        chromosome = std::string(chrStr);

        uint32_t physpos;
        fread(&physpos, 4, 1, m_fin);
        position = physpos;

        uint16_t K;
        fread(&K, 2, 1, m_fin);

        uint32_t LA;
        fread(&LA, 4, 1, m_fin);
        if (LA >= allele1.size()) {
            allele1.resize(2 * LA);
        }
        fread(allele1.data(), 1, LA, m_fin); allele1[LA] = '\0';
        first_allele = std::string(allele1.data());

        uint32_t LB;
        fread(&LB, 4, 1, m_fin);
        if (LB >= allele0.size()) {
            allele0.resize(2 * LB);
        }
        fread(allele0.data(), 1, LB, m_fin); allele0[LB] = '\0';
        second_allele = std::string(allele0.data());

        uint32_t C;
        fread(&C, 4, 1, m_fin);
        if (C > m_zBuf.size()) m_zBuf.resize(C - 4);

        uint32_t D;
        fread(&D, 4, 1, m_fin);
        m_zBufLens = C - 4;
        m_bufLens = D;
        fread(&m_zBuf[0], 1, C - 4, m_fin);

        AC = 0;
        AF = 0;
        info = 0;

        if (m_bufLens > m_buf.size()) m_buf.resize(m_bufLens);

        Parse2(&m_buf[0], m_bufLens, &m_zBuf[0], m_zBufLens,
               RSID, dosages, AC, AF, t_indexForMissing, info,
               t_indexForNonZero, t_isImputation);

        // Default setting is "ref-first" (matching SAIGE 03-14-2021 update)
        t_alt = second_allele;
        t_ref = first_allele;
        t_marker = RSID;
        t_pd = position;
        t_chr = chromosome;
        t_altFreq = AF;
        t_altCounts = AC;
        t_imputeInfo = info;
        t_missingRate = (double)t_indexForMissing.size() / (double)m_N;

        if (m_AlleleOrder == "alt-first") {
            t_alt = first_allele;
            t_ref = second_allele;
            t_altFreq = 1 - t_altFreq;
            t_altCounts = t_altFreq * 2 * ((double)m_N - (double)t_indexForMissing.size());
            for (unsigned int i = 0; i < dosages.n_elem; i++) {
                if (dosages[i] >= 0) {  // don't flip missing values (-1)
                    dosages[i] = 2 - dosages[i];
                }
            }
        }

        // Store marker info for index maps
        m_chr.push_back(t_chr);
        m_pd.push_back(t_pd);
        m_ref.push_back(t_ref);
        m_alt.push_back(t_alt);
        m_MarkerInBgen.push_back(t_marker);

    } else {
        t_isBoolRead = false;
    }
}

// ============================================================
// getMarkerIDToIndex: build chr:pos:ref:alt -> index map
// ============================================================
std::unordered_map<std::string, uint32_t> BgenClass::getMarkerIDToIndex()
{
    std::unordered_map<std::string, uint32_t> idMap;
    for (uint32_t i = 0; i < m_chr.size(); i++) {
        std::string id1 = m_chr[i] + ":" + std::to_string(m_pd[i]) + ":"
                         + m_ref[i] + ":" + m_alt[i];
        idMap[id1] = i;
        std::string id2 = m_chr[i] + ":" + std::to_string(m_pd[i]) + ":"
                         + m_alt[i] + ":" + m_ref[i];
        if (idMap.find(id2) == idMap.end()) {
            idMap[id2] = i;
        }
    }
    return idMap;
}

// ============================================================
// getMarkerNameToIndex: build marker ID -> index map
// ============================================================
std::unordered_map<std::string, uint32_t> BgenClass::getMarkerNameToIndex()
{
    std::unordered_map<std::string, uint32_t> nameMap;
    for (uint32_t i = 0; i < m_MarkerInBgen.size(); i++) {
        nameMap[m_MarkerInBgen[i]] = i;
    }
    return nameMap;
}

// ============================================================
// closegenofile: close the BGEN file handle
// ============================================================
void BgenClass::closegenofile()
{
    if (m_fin) {
        fclose(m_fin);
        m_fin = nullptr;
    }
}

} // namespace BGEN


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
// setVCFobjInCPP: create and configure the global VCF object
// Mirrors SAIGE Main.cpp::setVCFobjInCPP
// ============================================================
void setVCFobjInCPP(const std::string& t_vcfFileName,
                     const std::string& t_vcfField,
                     std::vector<std::string>& t_SampleInModel)
{
    // Clean up any existing object
    if (ptr_gVCFobj != nullptr) {
        ptr_gVCFobj->closegenofile();
        delete ptr_gVCFobj;
        ptr_gVCFobj = nullptr;
    }

    ptr_gVCFobj = new VCF::VcfClass(t_vcfFileName, t_vcfField, t_SampleInModel);

    int n = ptr_gVCFobj->getN();
    std::cout << "n:\t" << n << std::endl;
}


// ============================================================
// setBGENobjInCPP: create and configure the global BGEN object
// Mirrors SAIGE Main.cpp::setBGENobjInCPP
// ============================================================
void setBGENobjInCPP(const std::string& t_bgenFileName,
                      const std::vector<std::string>& t_SampleInBgen,
                      std::vector<std::string>& t_SampleInModel,
                      const std::string& t_AlleleOrder)
{
    // Clean up any existing object
    if (ptr_gBGENobj != nullptr) {
        ptr_gBGENobj->closegenofile();
        delete ptr_gBGENobj;
        ptr_gBGENobj = nullptr;
    }

    std::cout << "t_SampleInBgen " << t_SampleInBgen.size() << std::endl;
    ptr_gBGENobj = new BGEN::BgenClass(t_bgenFileName,
                                        t_SampleInBgen,
                                        t_SampleInModel,
                                        t_AlleleOrder);

    int n = ptr_gBGENobj->getN();
    std::cout << "n:\t" << n << std::endl;
}


// ============================================================
// Unified_getOneMarker: unified dispatcher for genotype reading
//
// Ported from SAIGE Main.cpp. Supports "plink", "vcf", and "bgen".
// Throws for "pgen" (not yet implemented).
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
    } else if (t_genoType == "vcf") {
        if (ptr_gVCFobj == nullptr) {
            throw std::runtime_error(
                "Unified_getOneMarker: VCF object not initialized. "
                "Call setVCFobjInCPP first.");
        }

        isBoolRead = ptr_gVCFobj->getOneMarker(
            t_ref, t_alt, t_marker, t_pd, t_chr,
            t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
            t_isOutputIndexForMissing, t_indexForMissing,
            t_isOnlyOutputNonZero, t_indexForNonZero,
            t_GVec, t_isImputation);
    } else if (t_genoType == "bgen") {
        if (ptr_gBGENobj == nullptr) {
            throw std::runtime_error(
                "Unified_getOneMarker: BGEN object not initialized. "
                "Call setBGENobjInCPP first.");
        }

        ptr_gBGENobj->getOneMarker(t_gIndex_prev, t_gIndex,
                                    t_ref, t_alt, t_marker, t_pd, t_chr,
                                    t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
                                    t_isOutputIndexForMissing, t_indexForMissing,
                                    t_isOnlyOutputNonZero, t_indexForNonZero,
                                    isBoolRead, t_GVec, t_isImputation);
    } else if (t_genoType == "pgen") {
        throw std::runtime_error(
            "Unified_getOneMarker: PGEN format not yet supported in standalone C++ port.");
    } else {
        throw std::runtime_error(
            "Unified_getOneMarker: Unknown genotype type '" + t_genoType + "'. "
            "Supported types: plink, vcf, bgen");
    }

    return isBoolRead;
}


// ============================================================
// Unified_getMarkerCount: get total marker count from genotype file
// ============================================================
uint32_t Unified_getMarkerCount(std::string& t_genoType)
{
    if (t_genoType == "plink") {
        if (ptr_gPLINKobj == nullptr) {
            throw std::runtime_error("Unified_getMarkerCount: PLINK object not initialized.");
        }
        return ptr_gPLINKobj->getM();
    } else if (t_genoType == "vcf") {
        if (ptr_gVCFobj == nullptr) {
            throw std::runtime_error("Unified_getMarkerCount: VCF object not initialized.");
        }
        return ptr_gVCFobj->getM0();
    } else if (t_genoType == "bgen") {
        if (ptr_gBGENobj == nullptr) {
            throw std::runtime_error("Unified_getMarkerCount: BGEN object not initialized.");
        }
        return ptr_gBGENobj->getM0();
    } else {
        throw std::runtime_error(
            "Unified_getMarkerCount: Unsupported type: " + t_genoType);
    }
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
    } else if (t_genoType == "vcf") {
        if (ptr_gVCFobj == nullptr) {
            throw std::runtime_error(
                "Unified_getSampleSizeinGeno: VCF object not initialized.");
        }
        N0 = ptr_gVCFobj->getN0();
    } else if (t_genoType == "bgen") {
        if (ptr_gBGENobj == nullptr) {
            throw std::runtime_error(
                "Unified_getSampleSizeinGeno: BGEN object not initialized.");
        }
        N0 = ptr_gBGENobj->getN0();
    } else {
        throw std::runtime_error(
            "Unified_getSampleSizeinGeno: Unsupported type: " + t_genoType);
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
    } else if (t_genoType == "vcf") {
        if (ptr_gVCFobj == nullptr) {
            throw std::runtime_error(
                "Unified_getSampleSizeinAnalysis: VCF object not initialized.");
        }
        N = ptr_gVCFobj->getN();
    } else if (t_genoType == "bgen") {
        if (ptr_gBGENobj == nullptr) {
            throw std::runtime_error(
                "Unified_getSampleSizeinAnalysis: BGEN object not initialized.");
        }
        N = ptr_gBGENobj->getN();
    } else {
        throw std::runtime_error(
            "Unified_getSampleSizeinAnalysis: Unsupported type: " + t_genoType);
    }
    return N;
}


// ============================================================
// Unified_getMarkerIDToIndex: get marker ID to index map
// ============================================================
std::unordered_map<std::string, uint32_t> Unified_getMarkerIDToIndex(std::string& t_genoType)
{
    if (t_genoType == "plink") {
        if (ptr_gPLINKobj == nullptr) {
            throw std::runtime_error("Unified_getMarkerIDToIndex: PLINK object not initialized.");
        }
        return ptr_gPLINKobj->getMarkerIDToIndex();
    } else if (t_genoType == "vcf") {
        if (ptr_gVCFobj == nullptr) {
            throw std::runtime_error("Unified_getMarkerIDToIndex: VCF object not initialized.");
        }
        return ptr_gVCFobj->getMarkerIDToIndex();
    } else if (t_genoType == "bgen") {
        if (ptr_gBGENobj == nullptr) {
            throw std::runtime_error("Unified_getMarkerIDToIndex: BGEN object not initialized.");
        }
        return ptr_gBGENobj->getMarkerIDToIndex();
    } else {
        throw std::runtime_error(
            "Unified_getMarkerIDToIndex: Unsupported type: " + t_genoType);
    }
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
    } else if (t_genoType == "vcf") {
        if (ptr_gVCFobj != nullptr) {
            ptr_gVCFobj->closegenofile();
        }
    } else if (t_genoType == "bgen") {
        if (ptr_gBGENobj != nullptr) {
            ptr_gBGENobj->closegenofile();
        }
    } else {
        throw std::runtime_error(
            "closeGenoFile: Unsupported type: " + t_genoType);
    }
}
