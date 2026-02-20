# SAIGE Step 2: Gene-Level Association Testing — C++ Standalone Conversion

## Project Goal

Convert SAIGE-GENE Step 2 (region/gene-based association tests) from Rcpp C++ to **standalone C++** (no R dependency). The original C++ code lives in the SAIGE R package at `/SAIGE/src/`. Our job is to extract it, remove Rcpp wrappers, and build a standalone CLI.

## What This Agent CAN Do

- **WRITE** files only to `/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Step_2_Feb_11/`
- **READ** any file in the project (especially `/SAIGE/src/` for reference C++ code)
- **READ** Step 1 C++ code at `/Jan_30_comparison/code_copy/cpp_standalone/` for patterns
- **RUN** build commands (`make`, compile) within the Step_2 directory
- **RUN** the Step 2 binary for testing

## What This Agent CANNOT Do

- **DO NOT** modify files in `/SAIGE/` (original R package — read-only reference)
- **DO NOT** modify files in `/Jan_30_comparison/` (Step 1 — finished, read-only)
- **DO NOT** change the algorithm — port the exact logic from SAIGE's Rcpp code
- **DO NOT** add features not in the original SAIGE Step 2
- **DO NOT** use R or Rcpp in the standalone code

## Directory Structure

```
Step_2_Feb_11/
├── CLAUDE.md                  # This file
├── AI_INSTRUCTIONS.txt        # Detailed agent rules (create if needed)
├── SESSION_NOTES.txt          # Chronological session log
├── code_copy/
│   └── cpp_standalone/        # C++ standalone Step 2 implementation
│       ├── main.cpp           # CLI entry point (YAML config)
│       ├── saige_test.cpp/hpp # SAIGEClass: score test, getMarkerPval
│       ├── spa.cpp/hpp        # SPA dispatcher
│       ├── spa_binary.cpp/hpp # Binary trait SPA (Newton-Raphson, saddle prob)
│       ├── cct.cpp/hpp        # Cauchy combination test
│       ├── group_file.cpp/hpp # Group file parser (gene definitions)
│       ├── null_model_loader.cpp/hpp  # Load Step 1 output (JSON + .arma files)
│       ├── genotype_reader.cpp/hpp    # Unified genotype reading (PLINK first)
│       ├── skat.cpp/hpp       # SKAT/BURDEN/SKAT-O + Davies method
│       ├── er_binary.cpp/hpp  # Efficient resampling for rare variants
│       ├── UTIL.cpp/hpp       # Math utilities (shared with Step 1)
│       ├── getMem.cpp/hpp     # Memory reporting
│       ├── Makefile
│       ├── config_test.yaml   # Test configuration (Step 1 .arma model)
│       └── config_compare.yaml # Comparison config (converted .rda model)
├── test/
│   ├── R/
│   │   ├── convert_rda_to_arma.R  # .rda → .arma + JSON converter
│   │   ├── run_r_step2.R          # R SAIGE Step 2 reference runner
│   │   ├── run_r_step2.sh         # Shell wrapper for pixi invocation
│   │   └── r_step2_results.txt    # R output (128,858 markers)
│   ├── compare_step2.R            # Column-by-column R vs C++ comparison
│   ├── run_comparison.sh          # End-to-end orchestrator (4 steps)
│   ├── RESULTS_step2.txt          # Comparison verdict
│   ├── data/
│   │   └── nullmodel_from_rda/    # Converted null model (12 .arma + JSON)
│   ├── output/
│   │   └── cpp_compare_results.txt # C++ output for comparison
│   └── Cpp_compare/               # C++ checkpoint files
├── reference/
│   ├── CPP_STANDALONE_CALL_GRAPH.html
│   ├── STEP2_CALL_GRAPH_v2.html
│   └── STEP2_MATH_EXPLAINED.html  # Math visualization
└── function_test/                  # Per-function unit tests (CCT, SPA, etc.)
```

## Reference Source Files (READ-ONLY)

These are the SAIGE Rcpp C++ files we are porting:

| Original File | Lines | Standalone Target | Conversion Difficulty |
|---|---|---|---|
| `/SAIGE/src/SAIGE_test.cpp` | 1,257 | `saige_test.cpp` | HIGH — R::pchisq/pnorm, Rcpp::Environment |
| `/SAIGE/src/SAIGE_test.hpp` | 265 | `saige_test.hpp` | LOW — just header change |
| `/SAIGE/src/SPA.cpp` | 300 | `spa.cpp` | HIGH — Rcpp::List returns |
| `/SAIGE/src/SPA_binary.cpp` | 529 | `spa_binary.cpp` | HIGH — 12× Rcpp::List, R::pnorm |
| `/SAIGE/src/SPA_binary.hpp` | 16 | `spa_binary.hpp` | MEDIUM — change return types |
| `/SAIGE/src/Main.cpp` | 3,044 | `main.cpp` (partial) | MEDIUM — Rcpp::export, Rcpp::stop |
| `/SAIGE/src/Main.hpp` | 382 | — | LOW — header only |
| `/SAIGE/src/CCT.cpp` | 82 | `cct.cpp` | LOW — just header change |
| `/SAIGE/src/UTIL.cpp` | ~200 | `UTIL.cpp` | LOW — mostly pure C++ |
| `/SAIGE/src/PLINK.cpp` | ~500 | `genotype_reader.cpp` | MEDIUM — may reuse Step 1 code |
| `/SAIGE/src/ER_binary_func.cpp` | ~300 | `er_binary.cpp` | MEDIUM |

## Rcpp-to-Standalone Conversion Rules

### Header Change (ALL files)
```cpp
// REMOVE:
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// ADD:
#include <armadillo>
```

### Replace Rcpp::List with C++ Structs
```cpp
// REMOVE patterns like:
Rcpp::List result;
result["root"] = root;
result["Isconverge"] = converge;
return result;

// REPLACE with:
struct RootResult { double root; int niter; bool converged; };
return RootResult{root, niter, converge};
```

### Replace R Statistical Functions
| R Function | C++ Replacement |
|---|---|
| `R::pchisq(x, df, lower, log)` | `boost::math::cdf(complement(boost::math::chi_squared(df), x))` |
| `R::pnorm(x, 0, 1, lower, log)` | `boost::math::cdf(boost::math::normal(0,1), x)` or complement |
| `R::qnorm(p, 0, 1, lower, log)` | `boost::math::quantile(boost::math::normal(0,1), p)` |

### Replace Rcpp::stop with Exceptions
```cpp
// REMOVE: Rcpp::stop("error message");
// REPLACE: throw std::runtime_error("error message");
```

### Replace Rcpp::Environment
```cpp
// REMOVE:
Rcpp::Environment base_env("package:base");
Rcpp::Function set_seed_r = base_env["set.seed"];
set_seed_r(seed);

// REPLACE:
static std::mt19937 rng_engine(seed);
```

## Step 2 Logic Structure

```
┌─────────────────────────────────────────────────────────────────┐
│                        SAIGE Step 2 Flow                        │
└─────────────────────────────────────────────────────────────────┘

1. INITIALIZATION
   ├── Load null model from Step 1 (JSON manifest + .arma binary files)
   │   ├── tau (variance components)
   │   ├── mu, res, y, V (N-vectors)
   │   ├── X, XVX, XVX_inv, XXVX_inv, XV, XVX_inv_XV (matrices)
   │   └── S_a (score vector)
   ├── Load variance ratios (varianceRatio.txt)
   ├── Load sparse GRM (if used in Step 1)
   ├── Initialize SAIGEClass with all above
   └── Open genotype file (PLINK/VCF/BGEN)

2. SINGLE-VARIANT TESTING (mainMarkerInCPP)
   ├── For each marker in genotype file:
   │   ├── Read genotype vector G [N×1]
   │   ├── Filter: MAC ≥ threshold, missing rate ≤ cutoff
   │   ├── Impute missing genotypes (mean imputation)
   │   ├── Assign variance ratio based on MAC category
   │   ├── Score test:
   │   │   ├── S = G' * res                    (score statistic)
   │   │   ├── var_S = varRatio * G' * P * G   (variance, P = projection)
   │   │   ├── Z = S / sqrt(var_S)             (Z-statistic)
   │   │   └── p = 2 * pnorm(-|Z|)             (p-value, normal approx)
   │   ├── If binary trait AND |Z| > SPA_cutoff (default 2):
   │   │   └── SPA correction:
   │   │       ├── Find saddlepoint ζ: K'(ζ) = q  (Newton-Raphson)
   │   │       ├── K(t) = Σ log(1 - mu_i + mu_i * exp(g_i * t))
   │   │       ├── p = saddle_prob(ζ, q)
   │   │       └── Two-sided: p = p_upper + p_lower
   │   └── Output: CHR, POS, REF, ALT, AF, MAC, Beta, SE, p-value
   └── Write results to output file

3. GENE/REGION-BASED TESTING (mainRegionInCPP)
   ├── Parse group file (gene → variant list + annotations + weights)
   ├── For each gene/region:
   │   ├── Read all markers for this region
   │   ├── For each annotation mask (e.g., lof, missense:lof):
   │   │   ├── For each MAF threshold (e.g., 0.0001, 0.001, 0.01):
   │   │   │   ├── Select variants matching annotation + MAF filter
   │   │   │   ├── Collapse ultra-rare variants (MAC ≤ 10)
   │   │   │   ├── Build score vector S [m×1] and covariance G_adj [N×m]
   │   │   │   ├── Compute variance-covariance matrix:
   │   │   │   │   ├── P1 = sqrt(VR) * G_adj'          [m × N]
   │   │   │   │   ├── P2 = sqrt(VR) * Sigma_inv * G   [N × m]
   │   │   │   │   └── Sigma = P1 * P2                  [m × m]
   │   │   │   ├── Apply weights: Beta(MAF, 1, 25) by default
   │   │   │   ├── Run tests:
   │   │   │   │   ├── BURDEN: T = (Σ w_i * S_i)^2 / (w' * Sigma * w)
   │   │   │   │   ├── SKAT: Q = S' * W * S, p via Davies method
   │   │   │   │   └── SKAT-O: optimal ρ over BURDEN↔SKAT spectrum
   │   │   │   ├── SPA adjustment (binary traits):
   │   │   │   │   └── Resampling-based p-value correction
   │   │   │   └── Output: Region, Annotation, MAF, p_SKATO, p_Burden, p_SKAT
   │   │   └── End MAF loop
   │   └── End annotation loop
   ├── Cauchy combination test (CCT) across masks
   └── Write results to output file

4. OUTPUT FILES
   ├── {prefix}.txt                    # Main results (gene-level)
   ├── {prefix}.singleAssoc.txt        # Single-variant within groups
   └── {prefix}.markerList.txt         # Markers used per test (optional)
```

## Step 1 → Step 2 Interface

Step 2 loads Step 1 output via `null_model_loader.cpp`:

```
Step 1 Output Directory:
├── nullmodel.json           # Scalars: tau, alpha, traitType, n, p, sampleIDs
├── mu.arma                  # arma::vec [N] — fitted values
├── res.arma                 # arma::vec [N] — residuals
├── y.arma                   # arma::vec [N] — phenotype
├── V.arma                   # arma::vec [N] — variance weights
├── X.arma                   # arma::mat [N×p] — design matrix
├── XVX.arma                 # arma::mat [p×p]
├── XVX_inv.arma             # arma::mat [p×p]
├── XXVX_inv.arma            # arma::mat [N×p]
├── XV.arma                  # arma::mat [p×N]
├── XVX_inv_XV.arma          # arma::mat [N×p]
├── S_a.arma                 # arma::vec [p] — score vector
├── varianceRatio.txt        # Text: VR values by MAC
└── sparseGRM.mtx            # Optional: sparse GRM (MatrixMarket)
```

## Build Dependencies

Same as Step 1:
- C++17 compiler (clang++ or g++)
- Armadillo (linear algebra)
- OpenBLAS + LAPACK
- yaml-cpp (config parsing)
- Boost (math distributions: normal, chi-squared, cauchy, beta)
- SuperLU (sparse solver, for sparse GRM PCG in Step 2)

## Implementation Order

1. **null_model_loader** — Load Step 1 output (JSON + .arma files)
2. **saige_test** — SAIGEClass + score test (single variant)
3. **spa_binary** — SPA for binary traits
4. **spa** — SPA dispatcher
5. **cct** — Cauchy combination test
6. **genotype_reader** — PLINK reading (reuse Step 1 code)
7. **main** — CLI with mainMarkerInCPP (single-variant testing)
8. **group_file** — Group file parser
9. **main** — Add mainRegionInCPP (gene-level testing)
10. **er_binary** — Efficient resampling (rare variants)

## Pixi Environment (Running SAIGE R)

SAIGE R package is NOT installed in system R. It runs via pixi:

```bash
# Pixi manifest location
PIXI_MANIFEST="/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/code_copy/SAIGE_isolated/pixi.toml"

# Command pattern
~/.pixi/bin/pixi run --manifest-path=$PIXI_MANIFEST Rscript <script.R>

# R version in pixi: 4.4.3 (osx-64)
# SAIGE version: 1.5.1

# To reinstall SAIGE after R code changes:
~/.pixi/bin/pixi run --manifest-path=$PIXI_MANIFEST R CMD INSTALL --preclean /path/to/SAIGE
```

**Important**: System R (4.5.1) does NOT have SAIGE. All SAIGE R operations must use pixi.
The R test scripts in `test/R/` do NOT modify SAIGE code — they only call SAIGE functions.

## Armadillo Binary Format (.arma files)

Both vectors and matrices use `ARMA_MAT_BIN_FN008` header:
```
ARMA_MAT_BIN_FN008
<rows> <cols>
<raw little-endian double-precision data, column-major>
```
- Vectors: `rows=N, cols=1`
- Matrices: `rows=nrow, cols=ncol`, column-major storage (same as R)
- Do NOT use `ARMA_COL_BIN_FN008` — Armadillo C++ `vec::load()` expects MAT format

## Validation Status

### Single-Variant Testing: VALIDATED (Feb 18, 2026)

**EXACT MATCH**: 644,340/644,340 values identical across all 5 columns
(BETA, SE, Tstat, var, p.value) for 128,868 markers.

Comparison pipeline: `test/run_comparison.sh`
1. `convert_rda_to_arma.R` — Extract .rda null model → .arma + JSON
2. `run_r_step2.R` — Run R SAIGE Step 2 (LOCO=TRUE, chr 1)
3. `./saige-step2 config_compare.yaml` — Run C++ Step 2
4. `compare_step2.R` — Column-by-column comparison → RESULTS_step2.txt

Key requirement: Both R and C++ must use the **same null model**. The .rda has
LOCO=TRUE (per-chromosome models). The converter extracts chr 1 LOCO values.
R must also run with LOCO=TRUE so it uses the same chr 1 model.

### Region/Gene-Based Testing: NOT YET VALIDATED

Code is implemented (mainRegionInCPP, SKAT/BURDEN/SKAT-O, Davies method)
but no R reference output has been generated for comparison yet.

## Key Variables

| Variable | Type | Description |
|---|---|---|
| `m_mu` | vec [N] | Fitted values from null model |
| `m_res` | vec [N] | Residuals (y - mu) |
| `m_mu2` | vec [N] | Variance weights: mu*(1-mu) for binary, 1/tau[0] for quant |
| `m_varRatioVal` | double | Current variance ratio for test marker |
| `m_XVX` | mat [p×p] | Weighted covariate information matrix |
| `m_XXVX_inv` | mat [N×p] | Projection: X * (X'VX)^{-1} |
| `m_XV` | mat [p×N] | Weighted covariates transposed |
| `m_S_a` | vec [p] | Score vector from null model |
| `P1Mat` | mat [m×N] | Gene marker score contribution (built per region) |
| `P2Mat` | mat [N×m] | Gene marker variance contribution (built per region) |

## Lessons from Step 1 and Step 2

1. **RNG bypass**: R and C++ RNGs are incompatible. If Step 2 uses random numbers, use bypass files.
2. **QR sign conventions**: Minor differences are OK (< 0.01%).
3. **Inner loop convergence**: Must replicate R's exact iteration logic.
4. **Solver selection**: Match R's default (direct sparse solve, not PCG) unless configured otherwise.
5. **Float precision**: Step 1 uses `arma::fvec`/`arma::fmat` (32-bit) — Step 2 uses double precision throughout.
6. **LOCO model consistency**: When comparing R vs C++, both must use the same null model (same LOCO chromosome). The .rda contains per-chromosome LOCO models; our converter extracts chr 1. R must also run with LOCO=TRUE.
7. **Armadillo binary format**: Always use `ARMA_MAT_BIN_FN008` header, even for vectors. `ARMA_COL_BIN_FN008` is not recognized by `vec::load(arma_binary)`.
8. **Empty sampleIDs fallback**: When null model has no sampleIDs, the genotype reader uses all fam samples in order (identity mapping).
9. **Pixi required**: System R does not have SAIGE. All R SAIGE operations must go through pixi.
