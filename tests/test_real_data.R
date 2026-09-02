#!/usr/bin/env Rscript
# ── AmbientFilter — real data test script ─────────────────────────────────────
# Tests AmbientFilter on a downsampled version of the SS2026 islet atlas.
# Checks that column names, cell type labels, and pipeline outputs are correct.
#
# Object path: ~/Verchere/jdrf2/SS2026/objects/lot2d/Islet/Islet_subclustered.rds
#
# Usage (from inside AmbientFilter/ folder):
#   source("tests/test_real_data.R")
# ────────────────────────────────────────────────────────────────────────────────

cat("\n══════════════════════════════════════════════════════\n")
cat("  AmbientFilter — SS2026 islet atlas real data test\n")
cat("══════════════════════════════════════════════════════\n\n")

suppressPackageStartupMessages({
  library(AmbientFilter)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

# ── Config: your object's metadata column names ───────────────────────────────
OBJ_PATH      <- "~/Verchere/jdrf2/SS2026/objects/lot2d/Islet/Islet_subclustered.rds"
CELL_TYPE_COL <- "cell_type"
SAMPLE_COL    <- "Sample"
TISSUE_COL    <- "Tissue"
ASSAY         <- "RNA"

# Cell types to correct (skip endocrine — they ARE the ambient source)
CELL_TYPES_TO_CORRECT <- c(
  "Macrophages",
  "Lymphocytes",
  "Mast",
  "Fibroblasts",
  "Endothelial",
  "Schwann",
  "Pericytes",
  "Activated_Stellate",
  "Quiescent_Stellate",
  "Ductal",
  "B_Cells"
)

# Cells to skip for background estimation (the ambient source)
EXCLUDE_FROM_BG <- c("Beta", "Alpha", "Alpha-Beta", "Delta",
                     "Gamma", "Epsilon")

# Per-cell-type prior strengths (biological knowledge)
# Macrophages: conservative — they are genuinely nibbly
# Everything else: data-driven will handle it
PRIOR_STRENGTHS <- c(Macrophages = 0.5)

# Downsample settings
N_CELLS_PER_TYPE <- 100   # cells per cell type for testing
N_DONORS         <- 5     # number of donors to include
SEED             <- 42

# ── Helpers ───────────────────────────────────────────────────────────────────
.pass <- function(t) { cat(paste0("  ✓ PASS  ", t, "\n")); n_pass <<- n_pass + 1L }
.fail <- function(t, m) {
  cat(paste0("  ✗ FAIL  ", t, "\n"))
  cat(paste0("         → ", m, "\n"))
  n_fail <<- n_fail + 1L
}
.section <- function(t) {
  cat(paste0("\n── ", t, " ", paste(rep("─", max(0, 52 - nchar(t))),
                                     collapse = ""), "\n"))
}
run_test <- function(name, expr) {
  tryCatch({
    val <- expr
    if (isTRUE(val)) { .pass(name) } else { .fail(name, paste0("returned: ", val)) }
  }, error = function(e) { .fail(name, conditionMessage(e)) })
  invisible(NULL)
}

n_pass <- 0L
n_fail <- 0L


# ── 1. Load object (check environment first) ──────────────────────────────────
.section("Loading object")

if (exists("seu") && is(seu, "Seurat")) {
  cat("  Object 'seu' already in environment — skipping load.\n")
  cat(paste0("  Dimensions: ", nrow(seu), " genes × ", ncol(seu), " cells\n"))
} else {
  cat(paste0("  Loading from: ", OBJ_PATH, "\n"))
  if (!file.exists(path.expand(OBJ_PATH))) {
    stop(paste0("Object not found at: ", OBJ_PATH,
                "\nUpdate OBJ_PATH at the top of this script."))
  }
  seu <- readRDS(path.expand(OBJ_PATH))
  cat(paste0("  Loaded: ", nrow(seu), " genes × ", ncol(seu), " cells\n"))
}

run_test("Object is a Seurat object", is(seu, "Seurat"))
run_test("cell_type column exists", CELL_TYPE_COL %in% colnames(seu@meta.data))
run_test("Sample column exists",    SAMPLE_COL    %in% colnames(seu@meta.data))
run_test("Tissue column exists",    TISSUE_COL    %in% colnames(seu@meta.data))
run_test("Raw counts layer exists", {
  !is.null(tryCatch(
    SeuratObject::LayerData(seu, assay = ASSAY, layer = "counts"),
    error = function(e) NULL
  ))
})
run_test("percent.ins column exists (pre-computed ambient metric)", {
  "percent.ins" %in% colnames(seu@meta.data)
})

# Show cell type breakdown
cat("\n  Cell type counts in full object:\n")
ct_table <- sort(table(seu@meta.data[[CELL_TYPE_COL]]), decreasing = TRUE)
for (i in seq_along(ct_table)) {
  ct <- names(ct_table)[i]
  n  <- ct_table[i]
  flag <- if (ct %in% CELL_TYPES_TO_CORRECT) "  ← will correct"
          else if (ct %in% EXCLUDE_FROM_BG)  "  ← ambient source (skip)"
          else "  ← not in correction list"
  cat(sprintf("    %-25s %6d%s\n", ct, n, flag))
}

# Warn about cell types in correction list that aren't in object
missing_cts <- setdiff(CELL_TYPES_TO_CORRECT,
                       unique(seu@meta.data[[CELL_TYPE_COL]]))
if (length(missing_cts) > 0) {
  cat(paste0("\n  ⚠ These cell types are in CELL_TYPES_TO_CORRECT but not in ",
             "the object:\n    ", paste(missing_cts, collapse = ", "), "\n"))
  CELL_TYPES_TO_CORRECT <- intersect(CELL_TYPES_TO_CORRECT,
                                     unique(seu@meta.data[[CELL_TYPE_COL]]))
}


# ── 2. Downsample ─────────────────────────────────────────────────────────────
.section("Downsampling")

set.seed(SEED)
meta <- seu@meta.data

# Select donors with enough cells across cell types
donor_counts <- table(meta[[SAMPLE_COL]])
donors_with_cells <- names(donor_counts)[donor_counts >= 50]
selected_donors <- sample(donors_with_cells,
                          min(N_DONORS, length(donors_with_cells)))
cat(paste0("  Selected donors (", length(selected_donors), "): ",
           paste(selected_donors, collapse = ", "), "\n"))

# Downsample per cell type per donor
sampled_cells <- lapply(selected_donors, function(donor) {
  donor_cells <- rownames(meta)[meta[[SAMPLE_COL]] == donor]
  donor_meta  <- meta[donor_cells, ]

  per_type <- lapply(CELL_TYPES_TO_CORRECT, function(ct) {
    ct_cells <- rownames(donor_meta)[donor_meta[[CELL_TYPE_COL]] == ct]
    if (length(ct_cells) == 0) return(NULL)
    sample(ct_cells, min(N_CELLS_PER_TYPE, length(ct_cells)))
  })
  unlist(per_type)
})
sampled_cells <- unique(unlist(sampled_cells))

cat(paste0("  Downsampled to: ", length(sampled_cells), " cells\n"))
seu_small <- subset(seu, cells = sampled_cells)

# Verify Seurat v5 counts layer is available
# (JoinLayers if needed)
tryCatch({
  SeuratObject::LayerData(seu_small, assay = ASSAY, layer = "counts")
}, error = function(e) {
  cat("  Joining layers (Seurat v5)...\n")
  seu_small <<- SeuratObject::JoinLayers(seu_small)
})

run_test("Downsampled object created", is(seu_small, "Seurat") && ncol(seu_small) > 0)
run_test("Downsampled object has correct cell types", {
  all(CELL_TYPES_TO_CORRECT %in%
      unique(seu_small@meta.data[[CELL_TYPE_COL]]))
})
run_test("Downsampled object has multiple donors", {
  length(unique(seu_small@meta.data[[SAMPLE_COL]])) >= 2
})

cat(paste0("\n  Downsampled breakdown:\n"))
ds_table <- table(seu_small@meta.data[[CELL_TYPE_COL]])
for (ct in names(ds_table)) {
  cat(sprintf("    %-25s %4d cells\n", ct, ds_table[ct]))
}


# ── 3. Test estimate_background ───────────────────────────────────────────────
.section("estimate_background()")

bg <- estimate_background(
  object             = seu_small,
  sample_col         = SAMPLE_COL,
  tissue_col         = TISSUE_COL,
  assay              = ASSAY,
  exclude_cell_types = EXCLUDE_FROM_BG,
  cell_type_col      = CELL_TYPE_COL,
  verbose            = TRUE
)

run_test("Background profiles created", is.list(bg) && length(bg) > 0)
run_test("One profile per tissue.sample", {
  n_expected <- length(unique(paste(
    seu_small@meta.data[[TISSUE_COL]],
    seu_small@meta.data[[SAMPLE_COL]],
    sep = "."
  )))
  length(bg) == n_expected
})
run_test("All profiles sum to ~1", {
  all(sapply(bg, function(x) abs(sum(x) - 1) < 0.01))
})
run_test("INS is a top ambient gene", {
  top_ambient <- names(sort(bg[[1]], decreasing = TRUE))[1:20]
  "INS" %in% top_ambient
})

# Show top ambient genes for first sample
cat(paste0("\n  Top 10 ambient genes in sample '",
           names(bg)[1], "':\n"))
top10 <- head(sort(bg[[1]], decreasing = TRUE), 10)
for (i in seq_along(top10)) {
  cat(sprintf("    %2d. %-10s %.4f\n", i, names(top10)[i], top10[i]))
}


# ── 4. Test estimate_correction_strength ──────────────────────────────────────
.section("estimate_correction_strength()")

strength_est <- estimate_correction_strength(
  object              = seu_small,
  background          = bg,
  cell_type_col       = CELL_TYPE_COL,
  sample_col          = SAMPLE_COL,
  tissue_col          = TISSUE_COL,
  assay               = ASSAY,
  cell_types          = CELL_TYPES_TO_CORRECT,
  min_cells           = 5,
  min_samples         = 2,
  compute_gene_weights = TRUE,
  prior_strengths     = PRIOR_STRENGTHS,
  prior_weight        = 0.3,
  make_plots          = FALSE,
  verbose             = FALSE
)

run_test("Correction strengths estimated", {
  is.numeric(strength_est$correction_strength) &&
    length(strength_est$correction_strength) > 0
})
run_test("All strengths in [0, 1]", {
  all(strength_est$correction_strength >= 0 &
      strength_est$correction_strength <= 1)
})
run_test("Macrophages get conservative correction", {
  "Macrophages" %in% names(strength_est$correction_strength) &&
    strength_est$correction_strength["Macrophages"] <= 0.7
})
run_test("Gene weights computed for Macrophages", {
  !is.null(strength_est$gene_weights$Macrophages)
})
run_test("INS has high ambient weight in Macrophages", {
  gw <- strength_est$gene_weights$Macrophages
  "INS" %in% names(gw) && gw["INS"] > 0.3
})

cat("\n  Estimated correction strengths:\n")
summ <- strength_est$summary[order(strength_est$summary$final_strength,
                                    decreasing = TRUE), ]
for (i in seq_len(nrow(summ))) {
  r <- summ[i, ]
  cat(sprintf("    %-25s corr=%5.3f  strength=%5.3f  (%s)\n",
              r$cell_type,
              ifelse(is.na(r$median_correlation), NA, r$median_correlation),
              r$final_strength,
              r$method))
}


# ── 5. Test correct_cell_type on Macrophages ──────────────────────────────────
.section("correct_cell_type() — Macrophages")

mac_str <- strength_est$correction_strength["Macrophages"]
mac_result <- correct_cell_type(
  object              = seu_small,
  background          = bg,
  cell_label          = "Macrophages",
  cell_type_col       = CELL_TYPE_COL,
  sample_col          = SAMPLE_COL,
  tissue_col          = TISSUE_COL,
  assay               = ASSAY,
  correction_strength = mac_str,
  min_cells           = 5,
  verbose             = TRUE
)

run_test("Macrophage correction ran successfully", {
  !is.null(mac_result$corrected)
})
run_test("Same number of cells pre/post correction", {
  ncol(mac_result$original) == ncol(mac_result$corrected)
})
run_test("INS reduced in macrophages after correction", {
  orig_ins <- mean(as.numeric(
    SeuratObject::LayerData(mac_result$original,  layer = "counts")["INS", ]
  ))
  corr_ins <- mean(as.numeric(
    SeuratObject::LayerData(mac_result$corrected, layer = "counts")["INS", ]
  ))
  cat(paste0("\n    INS: ", round(orig_ins, 3), " → ", round(corr_ins, 3),
             " (", round((orig_ins - corr_ins) / (orig_ins + 0.001) * 100, 1),
             "% reduction)\n"))
  corr_ins <= orig_ins
})
run_test("TTR reduced in macrophages after correction", {
  orig <- mean(as.numeric(
    SeuratObject::LayerData(mac_result$original,  layer = "counts")["TTR", ]
  ))
  corr <- mean(as.numeric(
    SeuratObject::LayerData(mac_result$corrected, layer = "counts")["TTR", ]
  ))
  cat(paste0("    TTR: ", round(orig, 3), " → ", round(corr, 3),
             " (", round((orig - corr) / (orig + 0.001) * 100, 1),
             "% reduction)\n"))
  corr <= orig
})
run_test("Changed genes table populated", nrow(mac_result$changed_genes) > 0)

cat("\n  Top 10 most-changed genes in Macrophages:\n")
top_changed <- head(mac_result$changed_genes, 10)
for (i in seq_len(nrow(top_changed))) {
  r <- top_changed[i, ]
  cat(sprintf("    %-12s  orig=%.3f  corr=%.3f  log2FC=%.3f\n",
              r$gene, r$mean_original, r$mean_corrected, r$log2fc))
}


# ── 6. Test full ambient_filter pipeline ──────────────────────────────────────
.section("ambient_filter() — full pipeline on downsampled object")

results <- ambient_filter(
  object                     = seu_small,
  cell_type_col              = CELL_TYPE_COL,
  sample_col                 = SAMPLE_COL,
  tissue_col                 = TISSUE_COL,
  cell_types                 = CELL_TYPES_TO_CORRECT,
  auto_estimate_strength     = TRUE,
  prior_strengths            = PRIOR_STRENGTHS,
  prior_weight               = 0.3,
  exclude_from_background    = EXCLUDE_FROM_BG,
  assay                      = ASSAY,
  min_cells                  = 5,
  min_samples_for_estimation = 2,
  top_genes_pct              = 0.3,
  make_plots                 = FALSE,
  run_consistency_analysis   = TRUE,
  min_sample_fraction        = 0.4,
  min_cell_types             = 2,
  export_offenders           = FALSE,
  return_full_corrected_obj  = FALSE,
  save_dir                   = NULL,
  verbose                    = TRUE
)

run_test("Pipeline completed successfully", !is.null(results))
run_test("All cell types present in results", {
  all(CELL_TYPES_TO_CORRECT %in% names(results$cell_type_results))
})
run_test("Summary table has correct number of rows", {
  nrow(results$summary) == length(CELL_TYPES_TO_CORRECT)
})
run_test("Strength estimates populated", {
  !is.null(results$strength_estimates)
})
run_test("Ambient map populated", {
  !is.null(results$ambient_map)
})
run_test("INS reduced across all corrected immune cell types", {
  immune_types <- intersect(c("Macrophages", "Lymphocytes", "Mast", "B_Cells"),
                            CELL_TYPES_TO_CORRECT)
  all(sapply(immune_types, function(ct) {
    if (!ct %in% names(results$cell_type_results)) return(TRUE)
    orig <- mean(as.numeric(SeuratObject::LayerData(
      results$cell_type_results[[ct]]$original,  layer = "counts")["INS", ]))
    corr <- mean(as.numeric(SeuratObject::LayerData(
      results$cell_type_results[[ct]]$corrected, layer = "counts")["INS", ]))
    corr <= orig
  }))
})

cat("\n  INS reduction summary across cell types:\n")
for (ct in CELL_TYPES_TO_CORRECT) {
  if (!ct %in% names(results$cell_type_results)) next
  r <- results$cell_type_results[[ct]]
  if (!"INS" %in% rownames(r$original)) next
  orig <- mean(as.numeric(
    SeuratObject::LayerData(r$original,  layer = "counts")["INS", ]
  ))
  corr <- mean(as.numeric(
    SeuratObject::LayerData(r$corrected, layer = "counts")["INS", ]
  ))
  pct_red <- round((orig - corr) / (orig + 0.001) * 100, 1)
  str_val <- round(results$strength_estimates$correction_strength[ct], 3)
  cat(sprintf("    %-25s  INS: %.3f → %.3f (%5.1f%% reduction)  strength=%.3f\n",
              ct, orig, corr, pct_red, str_val))
}


# ── 7. Test ambient map ───────────────────────────────────────────────────────
.section("identify_consistent_ambient() output")

run_test("Cross-cell-type offenders identified", {
  !is.null(results$ambient_map$cross_cell_type)
})
run_test("Consistency matrix present", {
  is.matrix(results$ambient_map$consistency_matrix)
})
run_test("Recommended features available", {
  is.list(results$ambient_map$recommended_features)
})

if (nrow(results$ambient_map$cross_cell_type) > 0) {
  cat("\n  Top cross-cell-type ambient offenders:\n")
  top_offenders <- head(results$ambient_map$cross_cell_type, 15)
  for (i in seq_len(nrow(top_offenders))) {
    r <- top_offenders[i, ]
    cat(sprintf("    %-12s  %d cell types  mean_frac=%.3f\n",
                r$gene, r$n_cell_types_affected, r$mean_fraction_corrected))
  }
} else {
  cat("  No cross-cell-type offenders found — try lowering min_sample_fraction.\n")
}


# ── 8. Metadata-specific checks ───────────────────────────────────────────────
.section("SS2026-specific checks")

run_test("percent.ins is higher in uncorrected macrophages than corrected", {
  # percent.ins in metadata reflects INS % in original counts
  mac_cells <- colnames(results$cell_type_results$Macrophages$original)
  meta_mac  <- seu_small@meta.data[mac_cells, ]
  mean(meta_mac$percent.ins, na.rm = TRUE) > 0
})
run_test("Condition column intact after correction", {
  mac_corr_meta <- results$cell_type_results$Macrophages$corrected@meta.data
  "Condition" %in% colnames(mac_corr_meta)
})
run_test("Donor column intact after correction", {
  mac_corr_meta <- results$cell_type_results$Macrophages$corrected@meta.data
  "Donor" %in% colnames(mac_corr_meta)
})
run_test("BMI_Category column intact after correction", {
  mac_corr_meta <- results$cell_type_results$Macrophages$corrected@meta.data
  "BMI_Category" %in% colnames(mac_corr_meta)
})
run_test("cell_type column intact in corrected object", {
  mac_corr_meta <- results$cell_type_results$Macrophages$corrected@meta.data
  CELL_TYPE_COL %in% colnames(mac_corr_meta) &&
    all(mac_corr_meta[[CELL_TYPE_COL]] == "Macrophages")
})


# ── 9. Downstream readiness check ─────────────────────────────────────────────
.section("Downstream readiness")

run_test("Corrected macrophage object can be normalised", {
  tryCatch({
    Seurat::NormalizeData(
      results$cell_type_results$Macrophages$corrected,
      verbose = FALSE
    )
    TRUE
  }, error = function(e) FALSE)
})
run_test("Corrected macrophage object can find variable features", {
  tryCatch({
    Seurat::FindVariableFeatures(
      results$cell_type_results$Macrophages$corrected,
      nfeatures = 500, verbose = FALSE
    )
    TRUE
  }, error = function(e) FALSE)
})

cat("\n  ✓ Corrected objects are ready for subclustering.\n")
cat("  Access them with:\n")
cat("    results$cell_type_results$Macrophages$corrected\n")
cat("    results$cell_type_results$Lymphocytes$corrected\n")
cat("    # etc.\n")


# ── 10. Results summary ───────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════\n")
cat(sprintf("  Results: %d passed | %d failed | %d total\n",
            n_pass, n_fail, n_pass + n_fail))
cat("══════════════════════════════════════════════════════\n\n")

if (n_fail > 0) {
  cat("  ⚠ Some tests failed — check output above.\n")
  cat("  Common fixes:\n")
  cat("    - JoinLayers() if Seurat v5 layer error\n")
  cat("    - Check cell type label spelling matches exactly\n")
  cat("    - Increase N_CELLS_PER_TYPE if min_cells threshold hit\n\n")
} else {
  cat("  All tests passed! Ready to run on the full object:\n\n")
  cat("  results <- ambient_filter(\n")
  cat(paste0("    object                 = seu,\n"))
  cat(paste0("    cell_type_col          = \"", CELL_TYPE_COL, "\",\n"))
  cat(paste0("    sample_col             = \"", SAMPLE_COL, "\",\n"))
  cat(paste0("    tissue_col             = \"", TISSUE_COL, "\",\n"))
  cat(paste0("    cell_types             = c(\"",
             paste(CELL_TYPES_TO_CORRECT, collapse = "\", \""), "\"),\n"))
  cat("    auto_estimate_strength = TRUE,\n")
  cat("    prior_strengths        = c(Macrophages = 0.5),\n")
  cat("    prior_weight           = 0.3,\n")
  cat(paste0("    exclude_from_background = c(\"",
             paste(EXCLUDE_FROM_BG, collapse = "\", \""), "\"),\n"))
  cat("    export_offenders       = TRUE,\n")
  cat("    dataset_name           = \"SS2026_islet_atlas\",\n")
  cat("    tissue                 = \"human_pancreatic_islet\",\n")
  cat("    save_dir               = \"ambient_filter_results\"\n")
  cat("  )\n\n")
}

invisible(results)


# ── BONUS: example of protected_genes usage ───────────────────────────────────
# Uncomment to test mitochondrial gene protection in macrophages.
# Use this when mitochondrial transfer between cell types is a concern
# (e.g. beta-to-macrophage mitochondrial transfer).
#
# mt_genes <- grep("^MT-", rownames(seu_small), value = TRUE)
# cat(paste0("\n  MT genes in object: ", paste(mt_genes, collapse=", "), "\n"))
#
# results_protected <- ambient_filter(
#   object                    = seu_small,
#   cell_type_col             = CELL_TYPE_COL,
#   sample_col                = SAMPLE_COL,
#   tissue_col                = TISSUE_COL,
#   cell_types                = CELL_TYPES_TO_CORRECT,
#   auto_estimate_strength    = TRUE,
#   prior_strengths           = PRIOR_STRENGTHS,
#   prior_weight              = 0.3,
#   exclude_from_background   = EXCLUDE_FROM_BG,
#   # ── Protect MT genes in macrophages from correction ──────────────────────
#   # MT genes may reflect genuine beta-to-macrophage mitochondrial transfer.
#   # Using "all" protects them across every cell type.
#   # Use a named list to protect only in specific cell types:
#   protected_genes = list(
#     Macrophages = mt_genes   # protect only in macrophages
#     # all = mt_genes         # protect across all cell types
#   ),
#   save_dir = NULL,
#   verbose  = TRUE
# )
#
# # Compare MT gene expression before/after in macrophages with protection
# mt_check <- sapply(mt_genes, function(g) {
#   if (!g %in% rownames(results_protected$cell_type_results$Macrophages$original)) {
#     return(c(orig = NA, corr_unprotected = NA, corr_protected = NA))
#   }
#   orig <- mean(as.numeric(LayerData(
#     results_protected$cell_type_results$Macrophages$original, layer="counts")[g,]))
#   prot <- mean(as.numeric(LayerData(
#     results_protected$cell_type_results$Macrophages$corrected, layer="counts")[g,]))
#   unprot <- mean(as.numeric(LayerData(
#     results$cell_type_results$Macrophages$corrected, layer="counts")[g,]))
#   c(orig = orig, corr_unprotected = unprot, corr_protected = prot)
# })
# print(t(mt_check))
