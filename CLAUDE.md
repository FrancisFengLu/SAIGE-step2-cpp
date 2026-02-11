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
│       ├── er_binary.cpp/hpp  # Efficient resampling for rare variants
│       ├── UTIL.cpp/hpp       # Math utilities (shared with Step 1)
│       ├── getMem.cpp/hpp     # Memory reporting
│       ├── Makefile
│       └── config_test.yaml   # Test configuration
└── output/
    ├── checkpoints/           # R vs C++ comparison files
    └── test_results/          # Step 2 output (p-values, etc.)
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

## Validation Strategy

Same as Step 1: checkpoint comparison.
- R Step 2 outputs reference p-values
- C++ Step 2 outputs same format
- Compare: differences < 0.01% are acceptable (floating-point precision)
- Differences > 1% are bugs

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

## Lessons from Step 1

1. **RNG bypass**: R and C++ RNGs are incompatible. If Step 2 uses random numbers, use bypass files.
2. **QR sign conventions**: Minor differences are OK (< 0.01%).
3. **Inner loop convergence**: Must replicate R's exact iteration logic.
4. **Solver selection**: Match R's default (direct sparse solve, not PCG) unless configured otherwise.
5. **Float precision**: Step 1 uses `arma::fvec`/`arma::fmat` (32-bit) — Step 2 may need double precision for p-values.
