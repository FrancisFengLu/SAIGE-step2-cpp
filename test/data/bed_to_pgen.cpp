// bed_to_pgen.cpp
// Convert PLINK 1.9 .bed/.bim/.fam to PGEN 2.0 .pgen/.pvar/.psam
// Output PGEN mode 0x02 (basic variant-major, hard-calls only)
//
// PLINK .bed 2-bit encoding per sample (within each byte, LSB first):
//   00 = homozygous for A1 (alt allele in .bim col 5)
//   01 = missing
//   10 = heterozygous
//   11 = homozygous for A2 (ref allele in .bim col 6)
//
// PGEN mode 0x02 encoding:
//   00 = homozygous ref (REF in .pvar)
//   01 = heterozygous
//   10 = homozygous alt (ALT in .pvar)
//   11 = missing
//
// Since PLINK .bim lists A1 (alt) first and A2 (ref) second,
// and PGEN .pvar lists REF then ALT:
//   BED 00 (hom A1/alt) -> PGEN 10 (hom alt)
//   BED 01 (missing)    -> PGEN 11 (missing)
//   BED 10 (het)        -> PGEN 01 (het)
//   BED 11 (hom A2/ref) -> PGEN 00 (hom ref)
//
// Mapping per 2-bit pair: bed_geno -> pgen_geno
//   0 -> 2, 1 -> 3, 2 -> 1, 3 -> 0
//
// Usage: bed_to_pgen <plink_prefix> <output_prefix>

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <plink_prefix> <output_prefix>" << std::endl;
        return 1;
    }

    std::string plinkPrefix = argv[1];
    std::string outPrefix = argv[2];

    std::string bimFile = plinkPrefix + ".bim";
    std::string famFile = plinkPrefix + ".fam";
    std::string bedFile = plinkPrefix + ".bed";

    std::string pvarFile = outPrefix + ".pvar";
    std::string psamFile = outPrefix + ".psam";
    std::string pgenFile = outPrefix + ".pgen";

    // ---- Read .fam to get sample count ----
    std::ifstream fam(famFile);
    if (!fam.is_open()) {
        std::cerr << "Cannot open .fam file: " << famFile << std::endl;
        return 1;
    }
    std::vector<std::pair<std::string, std::string>> samples; // FID, IID
    std::string line;
    while (std::getline(fam, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string fid, iid;
        iss >> fid >> iid;
        samples.push_back({fid, iid});
    }
    fam.close();
    uint32_t N = samples.size();
    std::cout << "Samples: " << N << std::endl;

    // ---- Write .psam ----
    std::ofstream psam(psamFile);
    psam << "#FID\tIID" << std::endl;
    for (auto& s : samples) {
        psam << s.first << "\t" << s.second << std::endl;
    }
    psam.close();
    std::cout << "Wrote " << psamFile << std::endl;

    // ---- Read .bim and write .pvar ----
    std::ifstream bim(bimFile);
    if (!bim.is_open()) {
        std::cerr << "Cannot open .bim file: " << bimFile << std::endl;
        return 1;
    }
    std::ofstream pvar(pvarFile);
    pvar << "#CHROM\tPOS\tID\tREF\tALT" << std::endl;

    uint32_t M = 0;
    while (std::getline(bim, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string chr, id, cm, pos, a1, a2;
        iss >> chr >> id >> cm >> pos >> a1 >> a2;

        // .bim: A1 = alt (minor), A2 = ref (major)
        // .pvar: REF then ALT
        // Convert to uppercase
        std::transform(a1.begin(), a1.end(), a1.begin(), ::toupper);
        std::transform(a2.begin(), a2.end(), a2.begin(), ::toupper);

        pvar << chr << "\t" << pos << "\t" << id << "\t" << a2 << "\t" << a1 << std::endl;
        M++;
    }
    bim.close();
    pvar.close();
    std::cout << "Variants: " << M << std::endl;
    std::cout << "Wrote " << pvarFile << std::endl;

    // ---- Build byte-level translation table ----
    // For each byte in .bed, translate all 4 genotypes at once
    // BED byte encodes 4 samples: bits [1:0]=s0, [3:2]=s1, [5:4]=s2, [7:6]=s3
    // Mapping per 2-bit pair: 0->2, 1->3, 2->1, 3->0
    uint8_t xlate[256];
    for (int b = 0; b < 256; b++) {
        uint8_t out = 0;
        for (int s = 0; s < 4; s++) {
            uint8_t bed_geno = (b >> (s * 2)) & 0x03;
            uint8_t pgen_geno;
            switch (bed_geno) {
                case 0: pgen_geno = 2; break;  // hom A1(alt) -> hom alt
                case 1: pgen_geno = 3; break;  // missing -> missing
                case 2: pgen_geno = 1; break;  // het -> het
                case 3: pgen_geno = 0; break;  // hom A2(ref) -> hom ref
                default: pgen_geno = 3; break;
            }
            out |= (pgen_geno << (s * 2));
        }
        xlate[b] = out;
    }

    // ---- Read .bed and write .pgen ----
    FILE* bed = fopen(bedFile.c_str(), "rb");
    if (!bed) {
        std::cerr << "Cannot open .bed file: " << bedFile << std::endl;
        return 1;
    }

    // Validate .bed magic number
    uint8_t magic[3];
    fread(magic, 1, 3, bed);
    if (magic[0] != 0x6c || magic[1] != 0x1b || magic[2] != 0x01) {
        std::cerr << "Invalid .bed magic number or not SNP-major" << std::endl;
        fclose(bed);
        return 1;
    }

    FILE* pgen = fopen(pgenFile.c_str(), "wb");
    if (!pgen) {
        std::cerr << "Cannot open output .pgen file: " << pgenFile << std::endl;
        fclose(bed);
        return 1;
    }

    // Write PGEN header (mode 0x02)
    // Bytes 0-1: magic 0x6c 0x1b
    uint8_t pgen_magic[2] = {0x6c, 0x1b};
    fwrite(pgen_magic, 1, 2, pgen);

    // Byte 2: mode 0x02
    uint8_t mode = 0x02;
    fwrite(&mode, 1, 1, pgen);

    // Bytes 3-6: variant count
    fwrite(&M, 4, 1, pgen);

    // Bytes 7-10: sample count
    fwrite(&N, 4, 1, pgen);

    // Byte 11: header control (0x00 for basic mode, no nonref flags)
    uint8_t headerCtrl = 0x00;
    fwrite(&headerCtrl, 1, 1, pgen);

    // Now write variant records
    uint64_t bytesPerVariant = (N + 3) / 4;
    std::vector<uint8_t> bedBuf(bytesPerVariant);
    std::vector<uint8_t> pgenBuf(bytesPerVariant);

    for (uint32_t v = 0; v < M; v++) {
        fread(bedBuf.data(), 1, bytesPerVariant, bed);

        // Translate each byte
        for (uint64_t i = 0; i < bytesPerVariant; i++) {
            pgenBuf[i] = xlate[bedBuf[i]];
        }

        fwrite(pgenBuf.data(), 1, bytesPerVariant, pgen);
    }

    fclose(bed);
    fclose(pgen);

    std::cout << "Wrote " << pgenFile << " (mode 0x02, " << M << " variants, "
              << N << " samples, " << (12 + bytesPerVariant * M) << " bytes)" << std::endl;
    std::cout << "Done." << std::endl;

    return 0;
}
