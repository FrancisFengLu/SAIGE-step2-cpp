#!/usr/bin/env Rscript
# ==============================================================================
# SAIGE Step 2 — Checkpoint Extraction Script
# ==============================================================================
#
# Purpose:
#   Run SAIGE Step 2 on test data and output intermediate/checkpoint values
#   for comparison with the standalone C++ port.
#
# Usage:
#   Rscript run_step2_checkpoints.R
#
# Output:
#   Writes checkpoint files to the same directory as this script (test/R/).
#   All numeric values use 15 significant digits for full-precision comparison.
#
# Requires:
#   - SAIGE R package installed
#   - SAIGE test data in extdata/
#
# ==============================================================================

# --- Configuration -----------------------------------------------------------

# Resolve paths relative to this script's location
script_dir <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) {
  # Fallback: use current working directory
  getwd()
})

# Base paths
saige_root    <- normalizePath(file.path(script_dir, "..", "..", "..", "SAIGE"))
extdata_dir   <- file.path(saige_root, "extdata")
input_dir     <- file.path(extdata_dir, "input")
output_dir    <- file.path(extdata_dir, "output")
checkpoint_dir <- normalizePath(script_dir, mustWork = FALSE)

# Ensure checkpoint output directory exists
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

# Helper: write a value with full precision to a checkpoint file
write_checkpoint <- function(filename, ..., header = NULL) {
  filepath <- file.path(checkpoint_dir, filename)
  cat("  Writing checkpoint:", filepath, "\n")

  lines <- c()
  if (!is.null(header)) {
    lines <- c(lines, header)
  }

  args <- list(...)
  for (item in args) {
    if (is.character(item)) {
      lines <- c(lines, item)
    } else if (is.matrix(item) || is.data.frame(item)) {
      # Write matrix row by row, tab-delimited, 15 significant digits
      for (i in seq_len(nrow(item))) {
        row_vals <- sapply(item[i, ], function(v) {
          if (is.numeric(v)) formatC(v, digits = 15, format = "g") else as.character(v)
        })
        lines <- c(lines, paste(row_vals, collapse = "\t"))
      }
    } else if (is.numeric(item)) {
      lines <- c(lines, paste(formatC(item, digits = 15, format = "g"), collapse = "\t"))
    } else {
      lines <- c(lines, paste(as.character(item), collapse = "\t"))
    }
  }

  writeLines(lines, filepath)
}

# --- Preflight checks --------------------------------------------------------

cat("=== SAIGE Step 2 Checkpoint Extraction ===\n\n")

cat("Script directory:", checkpoint_dir, "\n")
cat("SAIGE root:     ", saige_root, "\n")
cat("extdata:        ", extdata_dir, "\n\n")

# Check SAIGE is installed
if (!requireNamespace("SAIGE", quietly = TRUE)) {
  cat("ERROR: SAIGE R package is not installed.\n\n")
  cat("To install SAIGE:\n")
  cat("  # Option 1: From GitHub\n")
  cat("  remotes::install_github('saigegit/SAIGE')\n\n")
  cat("  # Option 2: From local source\n")
  cat("  R CMD INSTALL /path/to/SAIGE\n\n")
  cat("  # Option 3: Using conda\n")
  cat("  conda install -c bioconda r-saige\n\n")
  stop("SAIGE package not found")
}

library(SAIGE)
cat("SAIGE version:", as.character(packageVersion("SAIGE")), "\n\n")

# --- Discover test data -------------------------------------------------------

cat("--- Discovering test data ---\n")

# List input files
cat("Input files:\n")
input_files <- list.files(input_dir, full.names = TRUE)
for (f in input_files) {
  cat("  ", basename(f), "\n")
}

# List output files (null model from Step 1)
cat("\nOutput (null model) files:\n")
output_files <- list.files(output_dir, full.names = TRUE, recursive = TRUE)
for (f in output_files) {
  cat("  ", f, "\n")
}

# --- Identify test data paths -------------------------------------------------
#
# SAIGE extdata typically contains:
#   input/nfam_100_nindep_0_step1_includeMoreCovars_catOffsetVar.bed/bim/fam
#   input/nfam_100_nindep_0_step2_*.bed/bim/fam (genotypes for association)
#   input/group_new_chrposa1a2.txt (group file for gene-based tests)
#   output/extdata_binary_*.rda (Step 1 null model, binary trait)
#   output/extdata_quantitative_*.rda (Step 1 null model, quantitative trait)

# Find PLINK files for Step 2 genotypes
step2_bfiles <- unique(gsub("\\.(bed|bim|fam)$", "",
  list.files(input_dir, pattern = "step2.*\\.(bed|bim|fam)$")))

if (length(step2_bfiles) == 0) {
  # Fallback: use any PLINK set
  step2_bfiles <- unique(gsub("\\.(bed|bim|fam)$", "",
    list.files(input_dir, pattern = "\\.(bed|bim|fam)$")))
}

cat("\nPLINK file stems for Step 2:\n")
for (b in step2_bfiles) {
  cat("  ", b, "\n")
}

# Find null model RDA files
rda_files <- list.files(output_dir, pattern = "\\.rda$", full.names = TRUE, recursive = TRUE)
cat("\nNull model .rda files:\n")
for (f in rda_files) {
  cat("  ", f, "\n")
}

# Find group file
group_files <- list.files(input_dir, pattern = "group", full.names = TRUE)
cat("\nGroup files:\n")
for (f in group_files) {
  cat("  ", f, "\n")
}

# Find variance ratio files
vr_files <- list.files(output_dir, pattern = "varianceRatio|varRatio", full.names = TRUE, recursive = TRUE)
cat("\nVariance ratio files:\n")
for (f in vr_files) {
  cat("  ", f, "\n")
}

# Find sparse GRM files
grm_files <- list.files(output_dir, pattern = "sparseGRM|sparse.*GRM|GRM.*sparse",
                        full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
if (length(grm_files) > 0) {
  cat("\nSparse GRM files:\n")
  for (f in grm_files) {
    cat("  ", f, "\n")
  }
}

cat("\n")

# ==============================================================================
# MAIN: Run Step 2 and extract checkpoints
# ==============================================================================

# We'll try binary trait first, then quantitative

run_step2_with_checkpoints <- function(trait_type = "binary") {

  cat("===================================================================\n")
  cat("Running Step 2 checkpoints for trait type:", trait_type, "\n")
  cat("===================================================================\n\n")

  prefix <- trait_type

  # --- Step A: Find the right null model and test data ----------------------

  # Find the null model .rda for this trait type
  model_rda <- grep(trait_type, rda_files, value = TRUE)
  if (length(model_rda) == 0) {
    cat("WARNING: No .rda file found for trait type:", trait_type, "\n")
    cat("Available .rda files:", paste(rda_files, collapse = ", "), "\n")
    return(invisible(NULL))
  }
  # Use the first match that looks like the GMMATmodel
  model_rda_gmmat <- grep("GMMATmodel", model_rda, value = TRUE)
  if (length(model_rda_gmmat) > 0) {
    model_rda <- model_rda_gmmat[1]
  } else {
    model_rda <- model_rda[1]
  }
  cat("Using null model:", model_rda, "\n")

  # Find the variance ratio file for this trait type
  vr_file <- grep(trait_type, vr_files, value = TRUE)
  if (length(vr_file) == 0) {
    vr_file <- vr_files
  }
  vr_file_match <- grep("varianceRatio\\.txt$", vr_file, value = TRUE)
  if (length(vr_file_match) > 0) {
    vr_file <- vr_file_match[1]
  } else if (length(vr_file) > 0) {
    vr_file <- vr_file[1]
  } else {
    cat("WARNING: No variance ratio file found for trait type:", trait_type, "\n")
    return(invisible(NULL))
  }
  cat("Using variance ratio file:", vr_file, "\n")

  # PLINK bfile for Step 2 genotypes
  step2_bfile <- file.path(input_dir, step2_bfiles[1])
  cat("Using PLINK bfile:", step2_bfile, "\n")

  # Group file (for region-based tests)
  group_file <- if (length(group_files) > 0) group_files[1] else NULL
  if (!is.null(group_file)) {
    cat("Using group file:", group_file, "\n")
  }

  # Sparse GRM
  sparse_grm_file <- NULL
  sparse_grm_sample_file <- NULL
  if (length(grm_files) > 0) {
    sparse_grm_file <- grep("\\.mtx$|\\.txt$", grm_files, value = TRUE)
    if (length(sparse_grm_file) > 0) {
      sparse_grm_file <- sparse_grm_file[1]
      # Look for corresponding sample file
      sparse_grm_sample_file <- gsub("\\.[^.]+$", ".sampleIDs.txt", sparse_grm_file)
      if (!file.exists(sparse_grm_sample_file)) {
        sparse_grm_sample_file <- gsub("\\.[^.]+$", "SampleID.txt", sparse_grm_file)
      }
      if (!file.exists(sparse_grm_sample_file)) {
        sparse_grm_sample_file <- NULL
        sparse_grm_file <- NULL
      }
    } else {
      sparse_grm_file <- NULL
    }
  }

  cat("\n")

  # --- Step B: Load the null model and extract internals --------------------

  cat("--- Checkpoint B: Loading null model internals ---\n")

  tryCatch({
    # Load the .rda file to inspect its contents
    env <- new.env()
    load(model_rda, envir = env)
    model_objects <- ls(env)
    cat("Objects in .rda file:", paste(model_objects, collapse = ", "), "\n")

    # The main model object is usually called 'modglmm' or similar
    # Let's find it
    model_obj <- NULL
    for (obj_name in model_objects) {
      obj <- get(obj_name, envir = env)
      if (is.list(obj)) {
        cat("  Object '", obj_name, "' is a list with names: ",
            paste(head(names(obj), 20), collapse = ", "), "\n")
        if ("theta" %in% names(obj) || "tau" %in% names(obj) ||
            "mu" %in% names(obj) || "residuals" %in% names(obj) ||
            "fitted.values" %in% names(obj)) {
          model_obj <- obj
          cat("  -> Using this as the null model object\n")
        }
      }
    }

    if (!is.null(model_obj)) {
      # Extract key dimensions and parameters
      n_samples <- length(model_obj$y)
      trait <- if (!is.null(model_obj$traitType)) model_obj$traitType else "unknown"
      tau <- model_obj$theta  # variance components (tau[0], tau[1])
      if (is.null(tau)) tau <- model_obj$tau
      mu <- model_obj$fitted.values
      if (is.null(mu)) mu <- model_obj$mu
      res <- model_obj$residuals
      if (is.null(res)) res <- model_obj$res

      # Number of covariates
      X <- model_obj$X
      p_covars <- if (!is.null(X)) ncol(X) else NA

      cat("\n  n =", n_samples, "\n")
      cat("  p (covariates) =", p_covars, "\n")
      cat("  traitType =", trait, "\n")
      cat("  tau =", paste(formatC(tau, digits = 15, format = "g"), collapse = ", "), "\n")

      # Checkpoint a: null_model_dimensions.txt
      write_checkpoint(
        paste0(prefix, "_null_model_dimensions.txt"),
        header = "parameter\tvalue",
        paste0("n\t", n_samples),
        paste0("p\t", p_covars),
        paste0("traitType\t", trait),
        paste0("tau0\t", formatC(tau[1], digits = 15, format = "g")),
        paste0("tau1\t", if (length(tau) > 1) formatC(tau[2], digits = 15, format = "g") else "NA"),
        paste0("SPA_Cutoff\t", 2)
      )

      # Checkpoint b: mu_first10.txt
      if (!is.null(mu) && length(mu) >= 10) {
        write_checkpoint(
          paste0(prefix, "_mu_first10.txt"),
          header = "index\tmu",
          paste(seq_len(min(10, length(mu))), "\t",
                formatC(mu[1:min(10, length(mu))], digits = 15, format = "g"),
                sep = "")
        )
      }

      # Checkpoint c: res_first10.txt
      if (!is.null(res) && length(res) >= 10) {
        write_checkpoint(
          paste0(prefix, "_res_first10.txt"),
          header = "index\tres",
          paste(seq_len(min(10, length(res))), "\t",
                formatC(res[1:min(10, length(res))], digits = 15, format = "g"),
                sep = "")
        )
      }

      # Checkpoint d: XVX_matrix.txt
      XVX <- model_obj$XVX
      if (is.null(XVX) && !is.null(model_obj$obj.noK)) {
        XVX <- model_obj$obj.noK$XVX
      }
      if (!is.null(XVX)) {
        write_checkpoint(paste0(prefix, "_XVX_matrix.txt"),
                         header = paste0("# XVX matrix [", nrow(XVX), " x ", ncol(XVX), "]"),
                         XVX)
      } else {
        cat("  WARNING: XVX matrix not found in model object\n")
        cat("  Available names:", paste(names(model_obj), collapse = ", "), "\n")
        # Try to compute it
        if (!is.null(X) && !is.null(mu)) {
          V <- mu * (1 - mu)  # binary trait
          XVX_computed <- t(X) %*% diag(V) %*% X
          write_checkpoint(paste0(prefix, "_XVX_matrix.txt"),
                           header = paste0("# XVX matrix (computed) [", nrow(XVX_computed), " x ", ncol(XVX_computed), "]"),
                           XVX_computed)
        }
      }
    }
  }, error = function(e) {
    cat("ERROR loading null model internals:", conditionMessage(e), "\n")
    cat(traceback(), "\n")
  })

  # --- Step C: Load variance ratios ------------------------------------------

  cat("\n--- Checkpoint C: Variance ratios ---\n")

  tryCatch({
    vr_data <- read.table(vr_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    cat("  Variance ratio file has", nrow(vr_data), "rows and", ncol(vr_data), "columns\n")
    cat("  Columns:", paste(names(vr_data), collapse = ", "), "\n")
    print(head(vr_data, 20))

    # Checkpoint e: variance_ratios.txt
    write_checkpoint(
      paste0(prefix, "_variance_ratios.txt"),
      header = paste(names(vr_data), collapse = "\t")
    )
    # Append data rows
    filepath <- file.path(checkpoint_dir, paste0(prefix, "_variance_ratios.txt"))
    write.table(vr_data, file = filepath, append = TRUE,
                sep = "\t", row.names = FALSE, col.names = FALSE,
                quote = FALSE)
  }, error = function(e) {
    cat("ERROR reading variance ratio file:", conditionMessage(e), "\n")
    # Try comma-separated
    tryCatch({
      vr_data <- read.csv(vr_file, stringsAsFactors = FALSE)
      cat("  (Read as CSV) Variance ratio file has", nrow(vr_data), "rows\n")
      print(head(vr_data, 20))
    }, error = function(e2) {
      cat("ERROR reading as CSV too:", conditionMessage(e2), "\n")
    })
  })

  # --- Step D: Run SAIGE Step 2 (single-variant) ----------------------------

  cat("\n--- Checkpoint D: Running Step 2 single-variant test ---\n")

  step2_output <- file.path(checkpoint_dir, paste0(prefix, "_step2_output.txt"))

  tryCatch({
    # Use SAIGE's SPAGMMATtest function
    # The exact function name depends on SAIGE version
    # Try the newer API first (SAIGE >= 1.0)

    cat("  Attempting SPATests_Score_SPA_MAIN...\n")

    # First, we need to set up the model via SAIGE's internal functions
    # SAIGE Step 2 typically uses:
    #   1. SAIGE.Pheno.readGMMATmodelFromR() or similar to load null model
    #   2. SPAGMMATtest() to run tests

    # Check what functions are available
    saige_exports <- getNamespaceExports("SAIGE")
    cat("  Available SAIGE exports (first 50):\n")
    cat("  ", paste(head(sort(saige_exports), 50), collapse = ", "), "\n\n")

    # Look for the Step 2 main function
    step2_funcs <- grep("SPAGMMATtest|ScoreTest_wSPA|Step2|step2|mainMarker|mainRegion|SPA_G",
                        saige_exports, value = TRUE, ignore.case = TRUE)
    cat("  Step 2-related functions:", paste(step2_funcs, collapse = ", "), "\n\n")

    # Also look for model loading functions
    load_funcs <- grep("readGMMATmodel|loadNull|Load|GLMM|setgeno|Pheno",
                       saige_exports, value = TRUE, ignore.case = TRUE)
    cat("  Model loading functions:", paste(load_funcs, collapse = ", "), "\n\n")

    # Try to use SAIGE's documented Step 2 entry point
    # The standard call is:
    # SPAGMMATtest(
    #   vcfFile = "", vcfFileIndex = "", vcfField = "GT",
    #   savFile = "", savFileIndex = "",
    #   bgenFile = "", bgenFileIndex = "",
    #   sampleFile = "",
    #   bedFile = step2_bfile, bimFile = ..., famFile = ...,
    #   GMMATmodelFile = model_rda,
    #   varianceRatioFile = vr_file,
    #   SAIGEOutputFile = step2_output,
    #   ...
    # )

    # Build argument list for SPAGMMATtest
    bedFile <- paste0(step2_bfile, ".bed")
    bimFile <- paste0(step2_bfile, ".bim")
    famFile <- paste0(step2_bfile, ".fam")

    cat("  bedFile:", bedFile, "\n")
    cat("  bimFile:", bimFile, "\n")
    cat("  famFile:", famFile, "\n")
    cat("  Exist? bed:", file.exists(bedFile),
        "bim:", file.exists(bimFile),
        "fam:", file.exists(famFile), "\n\n")

    # Check if SPAGMMATtest exists
    if ("SPAGMMATtest" %in% saige_exports) {
      cat("  Calling SPAGMMATtest()...\n\n")

      SPAGMMATtest(
        bedFile = bedFile,
        bimFile = bimFile,
        famFile = famFile,
        AlleleOrder = "alt-first",
        GMMATmodelFile = model_rda,
        varianceRatioFile = vr_file,
        SAIGEOutputFile = step2_output,
        IsOutputAFinCaseCtrl = TRUE,
        LOCO = FALSE,
        is_output_moreDetails = TRUE
      )

    } else if ("SPAGMMATtest_v2" %in% saige_exports) {
      cat("  Calling SPAGMMATtest_v2()...\n\n")

      SPAGMMATtest_v2(
        bedFile = bedFile,
        bimFile = bimFile,
        famFile = famFile,
        AlleleOrder = "alt-first",
        GMMATmodelFile = model_rda,
        varianceRatioFile = vr_file,
        SAIGEOutputFile = step2_output,
        IsOutputAFinCaseCtrl = TRUE,
        LOCO = FALSE
      )

    } else {
      cat("  WARNING: SPAGMMATtest not found in SAIGE exports.\n")
      cat("  Trying alternative approach via internal functions...\n")
    }

    # Read and display results
    if (file.exists(step2_output)) {
      results <- read.table(step2_output, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE)
      cat("\n  Step 2 results: ", nrow(results), " markers\n")
      cat("  Columns:", paste(names(results), collapse = ", "), "\n\n")
      print(head(results, 5))

      # Checkpoint f: marker_0_info.txt
      if (nrow(results) >= 1) {
        m0 <- results[1, ]
        info_cols <- intersect(names(m0),
          c("CHR", "POS", "SNPID", "Allele1", "Allele2",
            "AC_Allele2", "AF_Allele2", "MissingRate",
            "chr", "pos", "ref", "alt", "altFreq", "altCounts", "missingRate",
            "MarkerID", "rsid"))

        write_checkpoint(
          paste0(prefix, "_marker_0_info.txt"),
          header = paste(names(m0), collapse = "\t"),
          paste(sapply(m0, function(v) {
            if (is.numeric(v)) formatC(v, digits = 15, format = "g") else as.character(v)
          }), collapse = "\t")
        )
      }

      # Checkpoint h: marker_0_pval.txt
      if (nrow(results) >= 1) {
        m0 <- results[1, ]
        pval_cols <- intersect(names(m0),
          c("Tstat", "var", "p.value", "p.value.NA", "Is.SPA",
            "BETA", "SE", "Tstat_unadj", "p.value_unadj",
            "Beta", "SE", "pval", "pval.norm"))

        write_checkpoint(
          paste0(prefix, "_marker_0_pval.txt"),
          header = "field\tvalue",
          sapply(names(m0), function(col) {
            val <- m0[[col]]
            if (is.numeric(val)) {
              paste0(col, "\t", formatC(val, digits = 15, format = "g"))
            } else {
              paste0(col, "\t", as.character(val))
            }
          })
        )
      }
    }

  }, error = function(e) {
    cat("ERROR running Step 2:", conditionMessage(e), "\n")
    cat("Traceback:\n")
    traceback()
  })

  # --- Step E: Extract genotype data for first marker -----------------------

  cat("\n--- Checkpoint E: First marker genotype values ---\n")

  tryCatch({
    # Read BIM file to get marker info
    bimFile <- paste0(step2_bfile, ".bim")
    bim <- read.table(bimFile, stringsAsFactors = FALSE,
                      col.names = c("CHR", "SNP", "CM", "POS", "A1", "A2"))
    cat("  BIM has", nrow(bim), "markers\n")
    cat("  First marker:", bim$SNP[1], "at", bim$CHR[1], ":", bim$POS[1], "\n")

    # Read FAM file to get sample count
    famFile <- paste0(step2_bfile, ".fam")
    fam <- read.table(famFile, stringsAsFactors = FALSE)
    n_samples <- nrow(fam)
    cat("  FAM has", n_samples, "samples\n")

    # Read first marker genotype from BED file using SAIGE or PLINK-reading code
    # We can use the BEDMatrix package or read raw bytes

    bedFile <- paste0(step2_bfile, ".bed")

    # Try BEDMatrix if available
    if (requireNamespace("BEDMatrix", quietly = TRUE)) {
      cat("  Using BEDMatrix to read genotypes...\n")
      bed <- BEDMatrix::BEDMatrix(step2_bfile)
      geno_first <- bed[, 1]  # First marker
      geno_first20 <- geno_first[1:min(20, length(geno_first))]

      write_checkpoint(
        paste0(prefix, "_marker_0_geno_first20.txt"),
        header = "index\tgenotype",
        paste(seq_along(geno_first20), "\t",
              formatC(geno_first20, digits = 15, format = "g"),
              sep = "")
      )
    } else {
      cat("  BEDMatrix not available. Reading BED file manually...\n")

      # Read BED file manually
      # BED format: 3-byte magic number, then packed genotypes
      # Each byte stores 4 genotypes (2 bits each)
      # Encoding: 00=homA (0), 01=missing (NA), 10=het (1), 11=homB (2)

      con <- file(bedFile, "rb")
      magic <- readBin(con, "raw", 3)
      cat("  BED magic bytes:", paste(as.integer(magic), collapse = " "), "\n")

      # Number of bytes per marker (SNP-major mode)
      bytes_per_marker <- ceiling(n_samples / 4)

      # Read first marker
      raw_bytes <- readBin(con, "raw", bytes_per_marker)
      close(con)

      # Decode genotypes
      geno <- integer(n_samples)
      idx <- 1
      for (byte_i in seq_along(raw_bytes)) {
        byte_val <- as.integer(raw_bytes[byte_i])
        for (bit_j in 0:3) {
          if (idx > n_samples) break
          bits <- bitwAnd(bitwShiftR(byte_val, 2 * bit_j), 3L)
          geno[idx] <- switch(as.character(bits),
                              "0" = 0L,    # homozygous A1
                              "1" = NA,    # missing
                              "2" = 1L,    # heterozygous
                              "3" = 2L)    # homozygous A2
          idx <- idx + 1
        }
      }

      geno_first20 <- geno[1:min(20, length(geno))]

      write_checkpoint(
        paste0(prefix, "_marker_0_geno_first20.txt"),
        header = "index\tgenotype",
        paste(seq_along(geno_first20), "\t",
              ifelse(is.na(geno_first20), "NA",
                     formatC(geno_first20, digits = 15, format = "g")),
              sep = "")
      )
    }

  }, error = function(e) {
    cat("ERROR reading genotypes:", conditionMessage(e), "\n")
  })

  # --- Step F: Run region/gene-based test if group file exists ---------------

  if (!is.null(group_file)) {
    cat("\n--- Checkpoint F: Running Step 2 region/gene-based test ---\n")

    step2_region_output <- file.path(checkpoint_dir,
                                     paste0(prefix, "_step2_region_output.txt"))

    tryCatch({
      # Check for region-based testing function
      region_funcs <- grep("Region|region|SKAT|Burden|gene",
                           saige_exports, value = TRUE, ignore.case = TRUE)
      cat("  Region-related functions:", paste(region_funcs, collapse = ", "), "\n\n")

      if ("SPAGMMATtest" %in% saige_exports) {
        cat("  Calling SPAGMMATtest() with gene-based options...\n\n")

        SPAGMMATtest(
          bedFile = bedFile,
          bimFile = bimFile,
          famFile = famFile,
          AlleleOrder = "alt-first",
          GMMATmodelFile = model_rda,
          varianceRatioFile = vr_file,
          SAIGEOutputFile = step2_region_output,
          groupFile = group_file,
          annotation_in_groupTest = "lof:missense",
          maxMAF_in_groupTest = 0.01,
          is_output_markerList_in_groupTest = TRUE,
          LOCO = FALSE,
          IsOutputAFinCaseCtrl = TRUE
        )

        # Read and display region results
        if (file.exists(step2_region_output)) {
          region_results <- read.table(step2_region_output, header = TRUE,
                                       sep = "\t", stringsAsFactors = FALSE)
          cat("\n  Region results:", nrow(region_results), "rows\n")
          cat("  Columns:", paste(names(region_results), collapse = ", "), "\n\n")
          print(head(region_results, 5))

          # Save as checkpoint
          write_checkpoint(
            paste0(prefix, "_region_results.txt"),
            header = paste(names(region_results), collapse = "\t")
          )
          filepath <- file.path(checkpoint_dir, paste0(prefix, "_region_results.txt"))
          write.table(region_results, file = filepath, append = TRUE,
                      sep = "\t", row.names = FALSE, col.names = FALSE,
                      quote = FALSE)
        }
      }
    }, error = function(e) {
      cat("ERROR in region-based test:", conditionMessage(e), "\n")
    })
  }

  cat("\n--- Done with", trait_type, "---\n\n")
}

# --- Step G: Deep internal checkpoint extraction via C++ hooks ----------------
#
# SAIGE's R code calls C++ functions via Rcpp. We can intercept these by
# tracing the relevant R functions or by calling the internal C++ functions
# directly through the SAIGE namespace.

extract_deep_checkpoints <- function(model_rda, vr_file, step2_bfile, prefix) {

  cat("===================================================================\n")
  cat("Extracting deep internal checkpoints for:", prefix, "\n")
  cat("===================================================================\n\n")

  tryCatch({
    # Access SAIGE's internal namespace
    ns <- asNamespace("SAIGE")

    # List internal C++ functions
    cpp_funcs <- grep("^[A-Z].*_cpp$|^set_|^get_|^initialize|^Score|^SPA",
                      ls(ns), value = TRUE)
    cat("Internal C++ functions:\n")
    cat("  ", paste(head(sort(cpp_funcs), 40), collapse = ", "), "\n\n")

    # Try to initialize the null model the way SAIGE does internally
    init_funcs <- grep("init|Init|load|Load|set|readGMMAT",
                       ls(ns), value = TRUE, ignore.case = TRUE)
    cat("Init/load functions:\n")
    cat("  ", paste(sort(init_funcs), collapse = ", "), "\n\n")

    # Look for the function that loads the model and sets up internal state
    # In newer SAIGE, it's typically:
    #   Prepare_Score_GLMM_NULL_Model()
    #   or readGMMATmodel()

    glmm_funcs <- grep("GLMM|glmm|Null|null|Model|model",
                        ls(ns), value = TRUE)
    cat("GLMM/null model functions:\n")
    cat("  ", paste(sort(glmm_funcs), collapse = ", "), "\n\n")

    # Try to load the model using SAIGE's internal pathway
    if (exists("Prepare_Score_GLMM_NULL_Model", envir = ns)) {
      cat("  Found Prepare_Score_GLMM_NULL_Model, calling...\n")
      # This function typically sets up internal C++ state
    }

    # Try to access internal state after SPAGMMATtest has run
    # (it should have been called above)
    state_funcs <- grep("^get_|^Get", ls(ns), value = TRUE)
    cat("Getter functions:\n")
    cat("  ", paste(sort(state_funcs), collapse = ", "), "\n\n")

    # Extract whatever internal values we can
    for (func_name in state_funcs) {
      tryCatch({
        func <- get(func_name, envir = ns)
        if (is.function(func) && length(formals(func)) == 0) {
          result <- func()
          if (is.numeric(result) && length(result) <= 100) {
            cat("  ", func_name, "():", paste(head(formatC(result, digits = 6, format = "g"), 10), collapse = ", "), "\n")
          } else if (is.numeric(result)) {
            cat("  ", func_name, "(): length =", length(result), ", first few:", paste(head(formatC(result, digits = 6, format = "g"), 5), collapse = ", "), "\n")
          } else {
            cat("  ", func_name, "():", str(result), "\n")
          }
        }
      }, error = function(e) {
        # Skip functions that fail
      })
    }

  }, error = function(e) {
    cat("ERROR in deep checkpoint extraction:", conditionMessage(e), "\n")
  })
}

# ==============================================================================
# Execute
# ==============================================================================

# Wrap everything in tryCatch so partial results are still saved
tryCatch({
  # Run binary trait
  run_step2_with_checkpoints("binary")
}, error = function(e) {
  cat("\n*** FATAL ERROR in binary trait run:", conditionMessage(e), "***\n\n")
})

tryCatch({
  # Run quantitative trait
  run_step2_with_checkpoints("quantitative")
}, error = function(e) {
  cat("\n*** FATAL ERROR in quantitative trait run:", conditionMessage(e), "***\n\n")
})

# Deep checkpoint extraction
tryCatch({
  # Use the first available model for deep extraction
  if (length(rda_files) > 0) {
    binary_rda <- grep("binary", rda_files, value = TRUE)
    binary_vr <- grep("binary", vr_files, value = TRUE)
    if (length(binary_rda) > 0 && length(binary_vr) > 0) {
      extract_deep_checkpoints(
        binary_rda[1], binary_vr[1],
        file.path(input_dir, step2_bfiles[1]),
        "binary_deep"
      )
    }
  }
}, error = function(e) {
  cat("\n*** FATAL ERROR in deep extraction:", conditionMessage(e), "***\n\n")
})

# ==============================================================================
# Summary
# ==============================================================================

cat("\n===================================================================\n")
cat("Checkpoint files written to:", checkpoint_dir, "\n")
cat("===================================================================\n\n")

checkpoint_files <- list.files(checkpoint_dir, pattern = "\\.(txt|tsv)$")
for (f in checkpoint_files) {
  fpath <- file.path(checkpoint_dir, f)
  fsize <- file.info(fpath)$size
  cat(sprintf("  %-50s  %s bytes\n", f, fsize))
}

cat("\nDone.\n")
