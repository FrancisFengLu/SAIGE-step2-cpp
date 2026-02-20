// Standalone port of SAIGE/src/PLINK.hpp
// PLINK genotype reader -- all Rcpp dependencies removed
// Ported from: /SAIGE/src/PLINK.hpp and /SAIGE/src/PLINK.cpp

#ifndef GENOTYPE_READER_HPP
#define GENOTYPE_READER_HPP

#include <armadillo>
#include <string>
#include <vector>
#include <cstdio>
#include <unordered_map>

namespace PLINK {

class PlinkClass {
private:

    // added on 03/14/2021
    std::string m_AlleleOrder;  // "alt-first" or "ref-first"

    // information from bim file
    uint32_t m_M0, m_M;  // total markers, markers in analysis
    std::vector<std::string> m_chr;               // Chromosome code
    std::vector<std::string> m_MarkerInPlink;     // Variant identifier
    std::vector<float> m_gd;                      // Position in morgans or centimorgans
    std::vector<uint32_t> m_pd;                   // Base-pair coordinate (1-based)
    std::vector<std::string> m_alt;               // Allele 1 (clear bits in .bed; usually minor)
    std::vector<std::string> m_ref;               // Allele 2 (set bits in .bed; usually major)

    // information from fam file
    std::vector<std::string> m_SampleInPlink;
    uint32_t m_N0, m_N;
    unsigned long long int m_numBytesofEachMarker0, m_numBytesofEachMarker;

    // input file stream of .bed file
    FILE* m_fin;

    // PLINK files
    std::string m_bimFile, m_famFile, m_bedFile;
    std::vector<uint32_t> m_posSampleInPlink;

    // https://www.cog-genomics.org/plink/1.9/formats#bed
    // PLINK format
    // The two-bit genotype codes have the following meanings:
    // 00  Homozygous for first allele in .bim file
    // 01  Missing genotype
    // 10  Heterozygous
    // 11  Homozygous for second allele in .bim file
    static const unsigned char HOM_REF = 0x3;   // 0b11
    static const unsigned char HET = 0x2;       // 0b10
    static const unsigned char HOM_ALT = 0x0;   // 0b00
    static const unsigned char MISSING = 0x1;   // 0b01

    // Genotype mapping vectors
    std::vector<int8_t> m_genoMaps_alt_first = {2, -1, 1, 0};
    std::vector<int8_t> m_genoMaps_ref_first = {0, -1, 1, 2};

    // pipeline: OneMarkerG4 --> bufferG4 --> bufferG1 --> OneMarkerG1
    std::vector<unsigned char> m_OneMarkerG4;

    void readBimFile();
    void readFamFile();

    // extract geno (0,1,2,3) at specific pos (0,1,2,3) of address c (1 byte)
    inline void getGenotype(const unsigned char c, const uint32_t pos, size_t& geno) {
        geno = (c >> (pos << 1)) & 0x3;  // 0b11 = 0x3
    }

public:

    PlinkClass(std::string t_bimFile,
               std::string t_famFile,
               std::string t_bedFile,
               std::string t_AlleleOrder);

    // setup PlinkClass
    void setPlinkobj(std::string t_bimFile,
                     std::string t_famFile,
                     std::string t_bedFile);

    void setPosSampleInPlink(std::vector<std::string>& t_SampleInModel);

    void getOneMarker(uint64_t& t_gIndex_prev,
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
                      arma::vec& OneMarkerG1);

    // Convenience overload: simplified getOneMarker (like SAIGE PLINK.hpp inline overloads)
    void getOneMarker(uint64_t t_gIndex_prev,
                      uint64_t t_gIndex,
                      double& t_altFreq,
                      double& t_missingRate,
                      std::string& t_chr,
                      arma::vec& OneMarkerG1)
    {
        std::string ref, alt, marker;
        uint32_t pd;
        double altCounts, imputeInfo;
        std::vector<uint> indexForMissing, indexForNonZero;
        bool isOutputIndexForMissing = false;
        bool isOnlyOutputNonZero = false;
        bool isTrueGenotype = false;
        getOneMarker(t_gIndex_prev, t_gIndex, ref, alt, marker, pd, t_chr,
                     t_altFreq, altCounts, t_missingRate, imputeInfo,
                     isOutputIndexForMissing, indexForMissing,
                     isOnlyOutputNonZero, indexForNonZero,
                     isTrueGenotype, OneMarkerG1);
    }

    // Convenience overload: with indexForMissing output
    void getOneMarker(uint64_t t_gIndex_prev,
                      uint64_t t_gIndex,
                      double& t_altFreq,
                      double& t_missingRate,
                      std::vector<uint>& t_indexForMissing,
                      arma::vec& OneMarkerG1)
    {
        std::string ref, alt, marker, chr;
        uint32_t pd;
        double altCounts, imputeInfo;
        std::vector<uint> indexForNonZero;
        bool isOutputIndexForMissing = false;
        bool isOnlyOutputNonZero = false;
        bool isTrueGenotype = true;
        getOneMarker(t_gIndex_prev, t_gIndex, ref, alt, marker, pd, chr,
                     t_altFreq, altCounts, t_missingRate, imputeInfo,
                     isOutputIndexForMissing, t_indexForMissing,
                     isOnlyOutputNonZero, indexForNonZero,
                     isTrueGenotype, OneMarkerG1);
    }

    uint32_t getN0() { return m_N0; }
    uint32_t getN() { return m_N; }
    uint32_t getM0() { return m_M0; }
    uint32_t getM() { return m_M; }
    uint32_t getnumBytesofEachMarker0() { return m_numBytesofEachMarker0; }
    uint32_t getnumBytesofEachMarker() { return m_numBytesofEachMarker; }

    std::vector<std::string> getChrVec() { return m_chr; }

    // Build a map from "chr:pos:ref:alt" (and "chr:pos:alt:ref") to marker index
    // Used by group file parser to look up variant genotype indices
    std::unordered_map<std::string, uint32_t> getMarkerIDToIndex() {
        std::unordered_map<std::string, uint32_t> idMap;
        for (uint32_t i = 0; i < m_M0; i++) {
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

    void closegenofile();
};

} // namespace PLINK

// Global checkpoint flag for debugging
extern bool g_writeCheckpoints;
extern std::string g_checkpointDir;

// Global pointer to PLINK object (mirrors SAIGE's ptr_gPLINKobj in Main.cpp)
extern PLINK::PlinkClass* ptr_gPLINKobj;

// Unified dispatcher (like SAIGE Main.cpp)
// For now, only handles "plink" type. Throws for bgen/vcf/pgen.
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
                          bool t_isImputation);

// Helper: set up PLINK object (mirrors SAIGE Main.cpp::setPLINKobjInCPP)
void setPLINKobjInCPP(std::string t_bimFile,
                      std::string t_famFile,
                      std::string t_bedFile,
                      std::vector<std::string>& t_SampleInModel,
                      std::string t_AlleleOrder);

// Helper: get sample size from genotype file
uint32_t Unified_getSampleSizeinGeno(std::string& t_genoType);

// Helper: get sample size in analysis (after subsetting)
uint32_t Unified_getSampleSizeinAnalysis(std::string& t_genoType);

// Helper: close genotype file
void closeGenoFile(std::string& t_genoType);

// C++ version of which(). Note: start from 0, not 1
std::vector<unsigned int> whichCPP(std::vector<std::string>& strVec,
                                   std::string strValue);

#endif
