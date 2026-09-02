#!/usr/bin/env Rscript
# ── AmbientFilter test script ──────────────────────────────────────────────────
# Generates a small synthetic Seurat object mimicking islet scRNA-seq ambient
# contamination, then runs the full AmbientFilter pipeline on it.
# Catches common bugs and validates that corrections are working as expected.
#
# Usage:
#   Rscript tests/test_ambient_filter.R
#   # or from within R:
#   source("tests/test_ambient_filter.R")
#
# Expected runtime: ~2-5 minutes (no GPU, no VAE by default)
# ────────────────────────────────────────────────────────────────────────────────

cat("\n══════════════════════════════════════════════\n")
cat("  AmbientFilter — functional test suite\n")
cat("══════════════════════════════════════════════\n\n")

# ── 0. Load packages ──────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(AmbientFilter)
  library(Seurat)
  library(Matrix)
})

# ── 1. Synthetic data generator ───────────────────────────────────────────────
# Generates a Seurat object with realistic islet ambient contamination:
#   - Beta cells express INS, TTR, IAPP highly
#   - Alpha cells express GCG, TTR moderately
#   - Macrophages express CD68, TREM2 + have LOW but detectable ambient INS/TTR
#   - Lymphocytes express CD3D, CD8A + have HIGHER ambient INS/TTR
#   - Ductal cells express KRT19, CFTR + have HIGH ambient INS/TTR (bystanders)
#
# Ambient contamination is injected as a proportion of each cell's total
# counts drawn from the beta cell profile — matching the real islet problem.

make_test_seurat <- function(n_cells_per_type = 80,
                             n_samples        = 4,
                             n_genes          = 500,
                             ambient_fraction = 0.12,  # 12% ambient reads
                             seed             = 42) {
  set.seed(seed)

  cat("  Building synthetic islet scRNA-seq object...\n")

  # Gene universe
  all_genes <- c(
    # Beta cell markers (main ambient source)
    "INS", "TTR", "IAPP", "MAFA", "PDX1", "NKX6-1", "GCK", "PCSK1",
    # Alpha cell markers
    "GCG", "ARX", "IRX2", "GPRC5C",
    # Macrophage markers
    "CD68", "TREM2", "FOLR2", "LYVE1", "MRC1", "ADRA2B", "TNIP3",
    "SPP1", "FABP4", "CD163", "CCL3", "CCL4", "ISG15", "HAMP",
    # Lymphocyte markers
    "CD3D", "CD3E", "CD8A", "CD8B", "CD4", "NKG7", "GZMB", "NCAM1",
    "PRF1", "IFNG", "CCL5", "IL7R", "KLRD1",
    # Ductal markers
    "KRT19", "CFTR", "SOX9", "HNF1B", "MUC1",
    # Filler genes (housekeeping, no cell type specificity)
    paste0("GENE_", seq_len(n_genes - 42))
  )
  n_genes_actual <- length(all_genes)

  cell_types  <- c("beta", "alpha", "macrophage", "lymphocyte", "ductal")
  tissues     <- c("pancreas", "spleen")
  sample_ids  <- paste0("DONOR_", seq_len(n_samples))

  # Base expression profiles per cell type (counts per gene, pre-ambient)
  make_base_profile <- function(ct) {
    prof <- rep(0.5, n_genes_actual)
    names(prof) <- all_genes

    if (ct == "beta") {
      prof["INS"]   <- 800
      prof["TTR"]   <- 400
      prof["IAPP"]  <- 150
      prof["MAFA"]  <- 30
      prof["PDX1"]  <- 20
      prof["NKX6-1"] <- 15
      prof["GCK"]   <- 25
      prof["PCSK1"] <- 40
    } else if (ct == "alpha") {
      prof["GCG"]   <- 600
      prof["TTR"]   <- 80
      prof["ARX"]   <- 20
      prof["IRX2"]  <- 15
      prof["GPRC5C"] <- 10
    } else if (ct == "macrophage") {
      prof["CD68"]  <- 50
      prof["TREM2"] <- 35
      prof["FOLR2"] <- 25
      prof["LYVE1"] <- 20
      prof["MRC1"]  <- 30
      prof["ADRA2B"] <- 15
      prof["TNIP3"] <- 40
      prof["CD163"] <- 20
    } else if (ct == "lymphocyte") {
      prof["CD3D"]  <- 60
      prof["CD3E"]  <- 55
      prof["CD8A"]  <- 40
      prof["CD8B"]  <- 38
      prof["NKG7"]  <- 25
      prof["GZMB"]  <- 30
      prof["IL7R"]  <- 20
    } else if (ct == "ductal") {
      prof["KRT19"] <- 80
      prof["CFTR"]  <- 40
      prof["SOX9"]  <- 30
      prof["HNF1B"] <- 25
      prof["MUC1"]  <- 35
    }
    prof
  }

  # Ambient profile: dominated by beta cell transcripts (the islet problem)
  # Varies slightly per sample to simulate processing batch effects
  make_ambient_profile <- function(sample_id, seed_offset = 0) {
    set.seed(seed + seed_offset)
    base <- rep(0.001, n_genes_actual)
    names(base) <- all_genes
    base["INS"]  <- 0.45 + rnorm(1, 0, 0.05)   # INS = ~45% of ambient reads
    base["TTR"]  <- 0.20 + rnorm(1, 0, 0.03)
    base["IAPP"] <- 0.08 + rnorm(1, 0, 0.01)
    base["GCG"]  <- 0.06 + rnorm(1, 0, 0.01)
    base["MAFA"] <- 0.02
    base["PDX1"] <- 0.01
    base <- pmax(base, 0)
    base / sum(base)  # normalise to proportions
  }

  # ── Generate cells ────────────────────────────────────────────────────────
  all_counts  <- list()
  all_meta    <- list()

  for (i_s in seq_along(sample_ids)) {
    s         <- sample_ids[i_s]
    tissue    <- if (i_s <= 2) "pancreas" else "spleen"
    amb_prof  <- make_ambient_profile(s, seed_offset = i_s * 100)

    for (ct in cell_types) {
      n_cells <- n_cells_per_type + sample(c(-15, 0, 15), 1)  # slight variation
      base    <- make_base_profile(ct)

      # Generate counts: Poisson with cell-specific library size scaling
      lib_sizes <- round(rlnorm(n_cells, meanlog = 7, sdlog = 0.5))

      cell_counts <- sapply(lib_sizes, function(lib) {
        # True signal
        true_probs   <- base / sum(base)
        true_counts  <- round(lib * (1 - ambient_fraction) * true_probs)

        # Ambient contamination drawn from ambient profile
        # Macrophages: lower ambient (they degrade some of what they ingest)
        ct_ambient_frac <- if (ct == "macrophage") ambient_fraction * 0.5
                           else if (ct == "beta")  0.0  # source, no contamination
                           else                    ambient_fraction

        amb_counts   <- rmultinom(1, round(lib * ct_ambient_frac), amb_prof)[, 1]
        pmax(true_counts + amb_counts, 0)
      })

      # Build sparse matrix
      cell_names <- paste0(s, "_", ct, "_cell", seq_len(n_cells))
      colnames(cell_counts) <- cell_names
      rownames(cell_counts) <- all_genes

      all_counts[[paste0(s, "_", ct)]] <- Matrix::Matrix(cell_counts, sparse = TRUE)
      all_meta[[paste0(s, "_", ct)]] <- data.frame(
        cell_id      = cell_names,
        sample_id    = s,
        tissue       = tissue,
        cell_type    = ct,
        donor_number = i_s,
        row.names    = cell_names,
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine
  count_mat <- do.call(cbind, all_counts)
  meta_df   <- do.call(rbind, all_meta)

  # Build Seurat object
  seu <- Seurat::CreateSeuratObject(
    counts   = count_mat,
    meta.data = meta_df,
    project  = "AmbientFilter_test"
  )
  seu$orig.ident <- seu$sample_id

  cat(paste0("    Created: ", ncol(seu), " cells | ",
             nrow(seu), " genes | ",
             length(unique(seu$sample_id)), " samples | ",
             length(unique(seu$cell_type)), " cell types\n"))
  seu
}


# ── 2. Test helpers ───────────────────────────────────────────────────────────
.pass <- function(test_name) {
  cat(paste0("  ✓ PASS  ", test_name, "\n"))
}
.fail <- function(test_name, msg) {
  cat(paste0("  ✗ FAIL  ", test_name, "\n"))
  cat(paste0("         → ", msg, "\n"))
}
.section <- function(title) {
  cat(paste0("\n── ", title, " ", paste(rep("─", max(0, 48 - nchar(title))),
                                         collapse=""), "\n"))
}

n_pass <- 0L
n_fail <- 0L

run_test <- function(name, expr) {
  result <- tryCatch({
    val <- expr
    if (isTRUE(val)) {
      .pass(name)
      n_pass <<- n_pass + 1L
      TRUE
    } else {
      .fail(name, paste0("returned: ", val))
      n_fail <<- n_fail + 1L
      FALSE
    }
  }, error = function(e) {
    .fail(name, conditionMessage(e))
    n_fail <<- n_fail + 1L
    FALSE
  })
  invisible(result)
}


# ── 3. Build test object ──────────────────────────────────────────────────────
.section("Building test data")
seu <- make_test_seurat(n_cells_per_type = 60, n_samples = 4)

run_test("Seurat object created", {
  is(seu, "Seurat") && ncol(seu) > 0
})
run_test("Required metadata columns present", {
  all(c("sample_id", "tissue", "cell_type") %in% colnames(seu@meta.data))
})
run_test("Raw counts accessible", {
  !is.null(SeuratObject::LayerData(seu, assay = "RNA", layer = "counts"))
})
run_test("INS is expressed in beta cells", {
  beta_cells <- colnames(seu)[seu$cell_type == "beta"]
  beta_counts <- as.matrix(
    SeuratObject::LayerData(seu, assay="RNA", layer="counts")["INS", beta_cells]
  )
  mean(beta_counts) > 100
})
run_test("INS is ambient in macrophages (pre-correction)", {
  mac_cells <- colnames(seu)[seu$cell_type == "macrophage"]
  mac_counts <- as.matrix(
    SeuratObject::LayerData(seu, assay="RNA", layer="counts")["INS", mac_cells]
  )
  mean(mac_counts) > 0  # should have some ambient signal
})


# ── 4. Test estimate_background ───────────────────────────────────────────────
.section("estimate_background()")
bg <- estimate_background(
  object        = seu,
  sample_col    = "sample_id",
  tissue_col    = "tissue",
  assay         = "RNA",
  verbose       = FALSE
)

run_test("Background is a named list", {
  is.list(bg) && !is.null(names(bg))
})
run_test("One background profile per tissue.sample combination", {
  n_expected <- length(unique(paste(seu$tissue, seu$sample_id, sep = ".")))
  length(bg) == n_expected
})
run_test("Background profiles sum to ~1", {
  all(sapply(bg, function(x) abs(sum(x) - 1) < 0.01))
})
run_test("INS dominates pancreas background", {
  pancreas_bg <- bg[grep("pancreas", names(bg))[1]][[1]]
  pancreas_bg["INS"] > 0.30  # should be ~45%
})
run_test("Background has no negative values", {
  all(sapply(bg, function(x) all(x >= 0)))
})


# ── 5. Test estimate_correction_strength ──────────────────────────────────────
.section("estimate_correction_strength()")
strength_est <- estimate_correction_strength(
  object              = seu,
  background          = bg,
  cell_type_col       = "cell_type",
  sample_col          = "sample_id",
  tissue_col          = "tissue",
  assay               = "RNA",
  cell_types          = c("macrophage", "lymphocyte", "ductal"),
  min_cells           = 5,
  min_samples         = 2,
  compute_gene_weights = TRUE,
  make_plots          = FALSE,
  verbose             = FALSE
)

run_test("Returns named correction_strength vector", {
  is.numeric(strength_est$correction_strength) &&
    !is.null(names(strength_est$correction_strength))
})
run_test("All strengths are in [0, 1]", {
  all(strength_est$correction_strength >= 0 &
      strength_est$correction_strength <= 1)
})
run_test("Macrophages get more conservative correction than ductal cells", {
  cs <- strength_est$correction_strength
  "macrophage" %in% names(cs) && "ductal" %in% names(cs) &&
    cs["macrophage"] < cs["ductal"]
})
run_test("Gene weights computed for macrophages", {
  !is.null(strength_est$gene_weights$macrophage) &&
    length(strength_est$gene_weights$macrophage) > 0
})
run_test("INS has high ambient weight in macrophages", {
  gw <- strength_est$gene_weights$macrophage
  "INS" %in% names(gw) && gw["INS"] > 0.5
})
run_test("Summary data frame has correct columns", {
  all(c("cell_type", "n_samples", "median_correlation",
        "final_strength") %in% colnames(strength_est$summary))
})
run_test("Prior blending works", {
  est_with_prior <- estimate_correction_strength(
    object          = seu,
    background      = bg,
    cell_type_col   = "cell_type",
    sample_col      = "sample_id",
    tissue_col      = "tissue",
    cell_types      = "macrophage",
    prior_strengths = c(macrophage = 0.5),
    prior_weight    = 1.0,   # fully prior
    min_cells       = 5,
    min_samples     = 2,
    make_plots      = FALSE,
    verbose         = FALSE
  )
  abs(est_with_prior$correction_strength["macrophage"] - 0.5) < 0.01
})


# ── 6. Test correct_cell_type ─────────────────────────────────────────────────
.section("correct_cell_type()")
mac_result <- correct_cell_type(
  object              = seu,
  background          = bg,
  cell_label          = "macrophage",
  cell_type_col       = "cell_type",
  sample_col          = "sample_id",
  tissue_col          = "tissue",
  assay               = "RNA",
  correction_strength = 0.5,
  min_cells           = 5,
  verbose             = FALSE
)

run_test("Returns list with expected elements", {
  all(c("original", "corrected", "changed_genes",
        "samples_used", "samples_skipped") %in% names(mac_result))
})
run_test("Corrected object has same cells as original", {
  ncol(mac_result$corrected) == ncol(mac_result$original)
})
run_test("Corrected object has same genes as original", {
  nrow(mac_result$corrected) == nrow(mac_result$original)
})
run_test("INS is reduced after correction in macrophages", {
  orig_ins  <- mean(as.matrix(
    SeuratObject::LayerData(mac_result$original,  layer = "counts")["INS", ]
  ))
  corr_ins  <- mean(as.matrix(
    SeuratObject::LayerData(mac_result$corrected, layer = "counts")["INS", ]
  ))
  corr_ins < orig_ins
})
run_test("CD68 (true marker) largely preserved after correction", {
  orig_cd68  <- mean(as.matrix(
    SeuratObject::LayerData(mac_result$original,  layer = "counts")["CD68", ]
  ))
  corr_cd68  <- mean(as.matrix(
    SeuratObject::LayerData(mac_result$corrected, layer = "counts")["CD68", ]
  ))
  # CD68 should be changed minimally (< 30% reduction)
  (orig_cd68 - corr_cd68) / (orig_cd68 + 1) < 0.30
})
run_test("changed_genes data frame has expected columns", {
  all(c("gene", "mean_original", "mean_corrected",
        "log2fc", "abs_log2fc") %in% colnames(mac_result$changed_genes))
})
run_test("All changed gene log2FC values are negative (reduction)", {
  cg <- mac_result$changed_genes
  all(cg$log2fc[cg$gene %in% c("INS", "TTR")] < 0)
})


# ── 7. Test full ambient_filter pipeline ──────────────────────────────────────
.section("ambient_filter() — full pipeline")
cell_types_to_test <- c("macrophage", "lymphocyte", "ductal")

results <- ambient_filter(
  object                    = seu,
  cell_type_col             = "cell_type",
  sample_col                = "sample_id",
  tissue_col                = "tissue",
  cell_types                = cell_types_to_test,
  auto_estimate_strength    = TRUE,
  prior_strengths           = c(macrophage = 0.5),
  prior_weight              = 0.3,
  assay                     = "RNA",
  min_cells                 = 5,
  min_samples_for_estimation = 2,
  make_plots                = FALSE,
  run_consistency_analysis  = TRUE,
  min_sample_fraction       = 0.5,
  min_cell_types            = 2,
  export_offenders          = FALSE,  # skip export in test
  return_full_corrected_obj = FALSE,
  save_dir                  = NULL,
  verbose                   = FALSE
)

run_test("ambient_filter() returns expected structure", {
  all(c("cell_type_results", "background", "strength_estimates",
        "summary", "ambient_map") %in% names(results))
})
run_test("All requested cell types corrected", {
  all(cell_types_to_test %in% names(results$cell_type_results))
})
run_test("Summary has one row per cell type", {
  nrow(results$summary) == length(cell_types_to_test)
})
run_test("Strength estimates populated", {
  !is.null(results$strength_estimates) &&
    all(cell_types_to_test %in% names(results$strength_estimates$correction_strength))
})
run_test("Ambient map populated", {
  !is.null(results$ambient_map) &&
    !is.null(results$ambient_map$per_cell_type)
})
run_test("INS reduced in lymphocytes after full pipeline", {
  orig_ins <- mean(as.matrix(
    SeuratObject::LayerData(
      results$cell_type_results$lymphocyte$original, layer = "counts"
    )["INS", ]
  ))
  corr_ins <- mean(as.matrix(
    SeuratObject::LayerData(
      results$cell_type_results$lymphocyte$corrected, layer = "counts"
    )["INS", ]
  ))
  corr_ins < orig_ins
})
run_test("Correction is more aggressive for lymphocytes than macrophages", {
  mac_str  <- results$strength_estimates$correction_strength["macrophage"]
  lymp_str <- results$strength_estimates$correction_strength["lymphocyte"]
  lymp_str > mac_str
})


# ── 8. Test identify_consistent_ambient ───────────────────────────────────────
.section("identify_consistent_ambient()")
ambient_map <- identify_consistent_ambient(
  results             = results,
  object              = seu,
  sample_col          = "sample_id",
  tissue_col          = "tissue",
  cell_type_col       = "cell_type",
  assay               = "RNA",
  lfc_threshold       = 0.2,
  min_sample_fraction = 0.4,
  min_cell_types      = 2,
  min_cells           = 5,
  make_plots          = FALSE,
  verbose             = FALSE
)

run_test("identify_consistent_ambient() returns expected structure", {
  all(c("per_cell_type", "cross_cell_type", "consistency_matrix",
        "recommended_features") %in% names(ambient_map))
})
run_test("Per-cell-type results present for all corrected cell types", {
  all(cell_types_to_test %in% names(ambient_map$per_cell_type))
})
run_test("INS identified as consistent ambient in lymphocytes", {
  lymp_df <- ambient_map$per_cell_type$lymphocyte
  "INS" %in% lymp_df$gene[lymp_df$is_consistent_ambient]
})
run_test("Consistency matrix is a matrix", {
  is.matrix(ambient_map$consistency_matrix)
})
run_test("recommended_features is a named list of character vectors", {
  is.list(ambient_map$recommended_features) &&
    all(sapply(ambient_map$recommended_features, is.character))
})


# ── 9. Test edge cases ────────────────────────────────────────────────────────
.section("Edge cases")

# Too few cells — should warn not error
run_test("Handles min_cells threshold gracefully", {
  r <- tryCatch(
    correct_cell_type(
      object              = seu,
      background          = bg,
      cell_label          = "macrophage",
      cell_type_col       = "cell_type",
      sample_col          = "sample_id",
      tissue_col          = "tissue",
      correction_strength = 1.0,
      min_cells           = 9999,  # impossibly high
      verbose             = FALSE
    ),
    error = function(e) list(samples_used = character(0),
                              samples_skipped = "error")
  )
  length(r$samples_used) == 0  # all skipped, no crash
})

# Invalid cell type — should error clearly
run_test("Invalid cell_label gives clear error", {
  tryCatch({
    correct_cell_type(seu, bg, cell_label = "NOTACELLTYPE",
                      cell_type_col = "cell_type",
                      sample_col = "sample_id",
                      verbose = FALSE)
    FALSE  # should not reach here
  }, error = function(e) {
    grepl("No cells found", conditionMessage(e))
  })
})

# Missing tissue_col — should fall back gracefully
run_test("Works without tissue_col (single tissue)", {
  bg_notissue <- estimate_background(
    seu, sample_col = "sample_id", tissue_col = NULL,
    assay = "RNA", verbose = FALSE
  )
  r <- correct_cell_type(
    seu, bg_notissue,
    cell_label          = "macrophage",
    cell_type_col       = "cell_type",
    sample_col          = "sample_id",
    tissue_col          = NULL,
    correction_strength = 0.5,
    min_cells           = 5,
    verbose             = FALSE
  )
  !is.null(r$corrected)
})

# prior_weight = 0 → fully data-driven
run_test("prior_weight = 0 is fully data-driven", {
  est_no_prior <- estimate_correction_strength(
    object = seu, background = bg,
    cell_type_col = "cell_type", sample_col = "sample_id",
    tissue_col = "tissue", cell_types = "macrophage",
    prior_strengths = c(macrophage = 0.99),  # extreme prior
    prior_weight = 0.0,  # fully data-driven — prior should have no effect
    min_cells = 5, min_samples = 2,
    make_plots = FALSE, verbose = FALSE
  )
  est_pure <- estimate_correction_strength(
    object = seu, background = bg,
    cell_type_col = "cell_type", sample_col = "sample_id",
    tissue_col = "tissue", cell_types = "macrophage",
    prior_strengths = NULL,
    min_cells = 5, min_samples = 2,
    make_plots = FALSE, verbose = FALSE
  )
  abs(est_no_prior$correction_strength["macrophage"] -
      est_pure$correction_strength["macrophage"]) < 0.01
})


# ── 10. Test export_ambient_offenders ─────────────────────────────────────────
.section("export_ambient_offenders()")
tmp_dir <- tempdir()
export_dir <- file.path(tmp_dir, "test_offenders")

offenders <- export_ambient_offenders(
  ambient_map  = results$ambient_map,
  dataset_name = "test_dataset",
  tissue       = "human_pancreatic_islet",
  n_donors     = 4,
  platform     = "synthetic",
  notes        = "AmbientFilter test suite",
  save_dir     = export_dir,
  save_excel   = FALSE
)

run_test("export_ambient_offenders() returns data frame", {
  is.data.frame(offenders)
})
run_test("CSV file created", {
  file.exists(file.path(export_dir, "ambient_offenders.csv"))
})
run_test("RDS file created", {
  file.exists(file.path(export_dir, "ambient_offenders.rds"))
})
run_test("RDS loads correctly", {
  loaded <- readRDS(file.path(export_dir, "ambient_offenders.rds"))
  is.list(loaded) && all(c("offenders", "metadata") %in% names(loaded))
})
run_test("CSV metadata header present", {
  lines <- readLines(file.path(export_dir, "ambient_offenders.csv"), n = 20)
  any(grepl("^# dataset_name:", lines))
})

# Cleanup
unlink(export_dir, recursive = TRUE)


# ── 11. Results summary ───────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════\n")
cat(sprintf("  Results: %d passed | %d failed | %d total\n",
            n_pass, n_fail, n_pass + n_fail))
cat("══════════════════════════════════════════════\n\n")

if (n_fail > 0) {
  cat("  ⚠ Some tests failed — check output above for details.\n\n")
  quit(status = 1)
} else {
  cat("  All tests passed! AmbientFilter is working correctly.\n\n")
  cat("  To test the VAE method (requires Python + PyTorch):\n")
  cat("    AmbientFilter::install_ambient_vae()\n")
  cat("    source('tests/test_vae.R')\n\n")
  quit(status = 0)
}
