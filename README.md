# SAIGE-GENE Step 2: C++ Standalone

Standalone C++17 implementation of SAIGE-GENE Step 2 (gene-level and single-variant association testing), converted from the original R/Rcpp codebase to pure C++ with no R dependency.

Takes the null model output from SAIGE Step 1 and runs association tests against a genotype file.

## Features

- **Single-variant testing**: score test with variance ratio adjustment
- **Region/gene-based testing**: BURDEN, SKAT, SKAT-O via Davies method and Liu approximation
- **SPA correction**: saddlepoint approximation for binary traits (Newton-Raphson root finding, Lugannani-Rice formula)
- **Firth correction**: penalized logistic regression fallback for binary traits
- **Efficient resampling**: exact test for ultra-rare binary variants (MAC <= 4)
- **Cauchy combination test (CCT)**: combine p-values across annotation/MAF strata
- **Conditional analysis**: condition on specified markers
- **Sparse GRM support**: PCG solver for re-evaluation with sparse genetic relationship matrix
- **LD matrix output**: optional pairwise LD computation for regions

## Supported Genotype Formats

| Format | Library | Config Key |
|--------|---------|------------|
| PLINK (.bed/.bim/.fam) | Built-in | `plinkFile` |
| VCF/BCF/VCF.GZ | htslib | `vcfFile`, `vcfField` |
| BGEN v1.2 | zstd + zlib | `bgenFile`, `bgenSampleFile` |
| PGEN (.pgen/.pvar/.psam) | Built-in (modes 0x01/0x02) | `pgenFile`, `pvarFile`, `psamFile` |

All formats produce identical output for the same genotype data.

## Build

### Dependencies

- C++17 compiler (clang++ or g++)
- Armadillo + OpenBLAS + LAPACK
- yaml-cpp
- Boost (math: normal, chi-squared, cauchy, beta distributions)
- SuperLU (sparse solver for PCG)
- htslib (VCF/BCF)
- zstd + zlib (BGEN decompression)

### Compile

```bash
cd code_copy/cpp_standalone
make -j4
```

Produces the `saige-step2` binary.

## Usage

```bash
./saige-step2 config.yaml
```

### Config File

```yaml
# Required
modelFile: /path/to/nullmodel_dir          # Directory with nullmodel.json + .arma files
varianceRatioFile: /path/to/varianceRatio.txt
plinkFile: /path/to/plink_prefix           # Or vcfFile/bgenFile/pgenFile
outputFile: /path/to/output.txt

# Single-variant parameters
minMAF: 0
minMAC: 0.5
maxMissRate: 0.15
AlleleOrder: alt-first                     # or ref-first

# Region/gene-based testing (optional — omit for single-variant only)
groupFile: /path/to/group_file.txt
annotationList:
  - "lof"
  - "missense;lof"
maxMAFList:
  - 0.0001
  - 0.001
  - 0.01
MACCutoff_to_CollapseUltraRare: 10
weights_beta: [1, 25]
isSingleInGroupTest: true

# Conditional analysis (optional)
condition: ["rs1", "rs2"]

# Sparse GRM (optional)
sparseGRMFile: /path/to/sparseGRM.mtx

# LD matrix output (optional)
isLDMatrix: true
```

### Input

- **Null model directory**: Step 1 output containing `nullmodel.json` and 11 `.arma` binary files (mu, res, y, V, X, XVX, XVX_inv, XXVX_inv, XV, XVX_inv_XV, S_a)
- **Variance ratio file**: Tab-separated file from Step 1 (`value type nMarkers` format or legacy numeric format)
- **Genotype file**: PLINK/VCF/BGEN/PGEN
- **Group file** (for region tests): 3-line-per-gene format (variants, annotations, weights)

### Output

- `{prefix}.txt` — main results (single-variant or region-level)
- `{prefix}.singleAssoc.txt` — per-variant results within each region (when `isSingleInGroupTest: true`)

## Validation

Validated against R SAIGE 1.5.1 across 12 test configurations:

| Test | Description | Result |
|------|-------------|--------|
| Single-variant (quantitative) | 128,868 markers, 5 columns | EXACT (644,340/644,340) |
| Single-variant (binary) | 70 markers, SPA + Firth + ER | EXACT (350/350) |
| Region BURDEN/SKAT (quantitative) | 2 regions, 6 strata | EXACT |
| Region BURDEN/SKAT (binary) | 2 regions, SPA Phi adjustment | EXACT |
| Region SKAT-O | Liu moment-matching approximation | ~1% max relative error |
| Sparse GRM (scoreTestFast + PCG) | 128,868 markers | EXACT (644,340/644,340) |
| No-covariate-adjustment path | scoreTestFast_noadjCov | EXACT (644,340/644,340) |
| Conditional analysis | 3 conditioning markers, 10 columns | EXACT (644,325/644,325) |
| Multiple variance ratio categories | 2 MAC bins | EXACT (644,340/644,340) |
| VCF format | htslib, GT field | EXACT match with PLINK |
| BGEN format | zstd/zlib decompression | EXACT match with PLINK |
| PGEN format | Modes 0x01/0x02 | EXACT match with PLINK |

## Project Structure

```
code_copy/cpp_standalone/
  main.cpp              CLI entry point, mainMarkerInCPP, mainRegionInCPP
  saige_test.cpp/hpp    SAIGEClass: score tests, SPA dispatch, Firth, conditional
  spa.cpp/hpp           SPA dispatcher
  spa_binary.cpp/hpp    Binary trait SPA (Newton-Raphson, saddlepoint probability)
  cct.cpp/hpp           Cauchy combination test
  skat.cpp/hpp          SKAT/BURDEN/SKAT-O, Davies method
  group_file.cpp/hpp    Group file parser
  null_model_loader.cpp/hpp  Load Step 1 null model
  genotype_reader.cpp/hpp    Unified genotype reader (PLINK/VCF/BGEN/PGEN)
  er_binary.cpp/hpp     Efficient resampling for rare binary variants
  UTIL.cpp/hpp          Math utilities
  getMem.cpp/hpp        Memory reporting
  Makefile
test/
  R/                    R scripts for validation (convert .rda, run R SAIGE, compare)
  data/                 Test data (null models, group files)
  output/               C++ test output
```
