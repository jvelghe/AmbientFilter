# ── AmbientFilter — main entry point ──────────────────────────────────────────
# ────────────────────────────────────────────────────────────────────────────────

#' Per-cell-type ambient RNA correction
#'
#' Main entry point for the AmbientFilter pipeline. Estimates per-sample
#' background profiles from the full object, then iterates over each specified
#' cell type (or all cell types) applying [DropletUtils::removeAmbience()]
#' with cell-type-specific correction strengths. Corrected counts are written
#' back into a copy of the full Seurat object.
#'
#' Designed for islet and exocrine pancreas scRNA-seq datasets where beta cell
#' and alpha cell transcripts (INS, TTR, GCG) strongly contaminate the ambient
#' pool and are erroneously detected in non-expressing populations such as
#' macrophages, lymphocytes, mast cells, and ductal cells. Run this function
#' **upstream of subclustering**.
#'
#' @section Correction strength guidance:
#' | Cell type | Recommended strength | Rationale |
#' |-----------|---------------------|-----------|
#' | Lymphocytes, ductal, stellate, mast | 1.0 | No biological reason to contain ambient transcripts |
#' | Macrophages | 0.5 | Legitimately internalise ambient RNA via efferocytosis/trogocytosis |
#' | Beta cells, alpha cells | 0.0 or skip | They *are* the ambient source; correcting them removes real signal |
#'
#' @param object A Seurat object. Must have raw counts in the specified assay
#'   and cell type + sample columns in metadata.
#' @param cell_type_col Character. Metadata column containing cell type labels.
#'   Default `"cell_type"`.
#' @param sample_col Character. Metadata column containing sample/donor IDs.
#'   Default `"orig.ident"`.
#' @param tissue_col Character or NULL. Metadata column for tissue identity.
#'   When provided, background estimation is stratified per sample within
#'   tissue. Strongly recommended for multi-tissue datasets (pancreas +
#'   spleen etc.). Default `NULL`.
#' @param cell_types Character vector or NULL. Cell type labels to correct.
#'   If `NULL`, all unique values in `cell_type_col` are processed.
#'   Default `NULL`.
#' @param correction_strength Named numeric vector or single numeric.
#'   Per-cell-type correction strength in \[0, 1\]. A named vector maps cell
#'   type labels to strength values; unlisted cell types default to `1.0`.
#'   A single value applies to all cell types. See section above for guidance.
#'   Default `1.0`.
#' @param exclude_from_background Character vector or NULL. Cell type labels
#'   to exclude when computing the background estimate. Useful when one cell
#'   type is so numerically dominant that it inflates the ambient profile
#'   (e.g. exclude `"beta"` from background estimation if your sample is
#'   strongly beta-cell enriched). Note this affects the background estimate
#'   for *all* cell types, not just the excluded type itself. Default `NULL`.
#' @param assay Character. Seurat assay containing raw counts. Default `"RNA"`.
#' @param min_cells Integer. Minimum cells per cell-type/sample combination.
#'   Combinations with fewer cells are skipped. Default `10`.
#' @param features Named list or NULL. Per-cell-type gene lists to
#'   force-correct (passed to `features` in `removeAmbience`). Names should
#'   match `cell_types`. Default `NULL`.
#' @param protected_genes Named list or NULL. Per-cell-type gene lists to
#'   **exclude from correction**. Genes in this list retain their original
#'   counts regardless of ambient weight. Names should match `cell_types`,
#'   or use `"all"` to protect genes across all cell types. Examples:
#'   \cr\cr
#'   Protect mitochondrial genes in macrophages only (recommended when
#'   mitochondrial transfer between cell types is a possibility):
#'   \cr
#'   `protected_genes = list(Macrophages = grep("^MT-", rownames(seu), value = TRUE))`
#'   \cr\cr
#'   Protect MT genes across all cell types:
#'   \cr
#'   `protected_genes = list(all = grep("^MT-", rownames(seu), value = TRUE))`
#'   \cr\cr
#'   Default `NULL` (no protection — correction applied to all genes).
#' @param make_plots Logical. Generate diagnostic plots for each cell type.
#'   Default `TRUE`.
#' @param label_difference Numeric. |log2FC| threshold for labelling genes in
#'   correlation plots. Default `0.5`.
#' @param highlight_genes Character vector. Genes to always highlight in
#'   diagnostic plots (e.g. `c("INS", "TTR", "GCG")` for islet data).
#'   Default `c("INS", "TTR", "GCG", "PRSS1", "CELA3A")`.
#' @param return_full_corrected_obj Logical. If `TRUE`, writes all corrected
#'   counts back into a copy of the full input object. Memory-intensive for
#'   large objects. Default `FALSE` (only per-cell-type corrected subsets are
#'   returned). The full corrected object is always accessible at
#'   `results$updated_object` if this is `TRUE`.
#' @param save_dir Character or NULL. Directory to save plots and summary
#'   tables. Created if it does not exist. Default `NULL` (no saving).
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A list with:
#' \describe{
#'   \item{`cell_type_results`}{Named list (one entry per cell type) of
#'     correction results from [correct_cell_type()]. Each entry has
#'     `original`, `corrected`, `changed_genes`, `samples_used`,
#'     `samples_skipped`.}
#'   \item{`background`}{The per-sample background profiles from
#'     [estimate_background()].}
#'   \item{`plots`}{Named list of diagnostic plots from
#'     [plot_correction_diagnostics()] (if `make_plots = TRUE`).}
#'   \item{`summary`}{Data frame from [summarise_correction()] summarising
#'     correction statistics per cell type.}
#'   \item{`updated_object`}{The full Seurat object with corrected counts
#'     written back for all processed cell types (only if
#'     `return_full_corrected_obj = TRUE`).}
#'   \item{`params`}{List of parameters used for this run, for
#'     reproducibility.}
#' }
#'
#' @examples
#' \dontrun{
#' library(AmbientFilter)
#'
#' # Basic islet dataset — full correction for immune cells,
#' # conservative for macrophages
#' results <- ambient_filter(
#'   object         = seu,
#'   cell_type_col  = "cell_type",
#'   sample_col     = "hpap_id",
#'   tissue_col     = "tissue",
#'   cell_types     = c("macrophage", "lymphocyte", "mast_cell",
#'                      "ductal", "stellate", "endothelial"),
#'   correction_strength = c(
#'     macrophage   = 0.5,   # nibbly — conservative
#'     lymphocyte   = 1.0,
#'     mast_cell    = 1.0,
#'     ductal       = 1.0,
#'     stellate     = 1.0,
#'     endothelial  = 1.0
#'   ),
#'   exclude_from_background = NULL,
#'   highlight_genes  = c("INS", "TTR", "GCG", "PRSS1"),
#'   make_plots       = TRUE,
#'   save_dir         = "ambient_filter_results"
#' )
#'
#' # Inspect summary
#' results$summary
#'
#' # Access corrected macrophage object for subclustering
#' mac_corrected <- results$cell_type_results$macrophage$corrected
#'
#' # Compare INS expression before and after in macrophages
#' results$plots$macrophage$violin_original
#' results$plots$macrophage$violin_corrected
#' }
#'
#' @export
ambient_filter <- function(object,
                           cell_type_col              = "cell_type",
                           sample_col                 = "orig.ident",
                           tissue_col                 = NULL,
                           cell_types                 = NULL,
                           correction_strength        = NULL,
                           auto_estimate_strength     = TRUE,
                           prior_strengths            = NULL,
                           prior_weight               = 0.3,
                           top_genes_pct              = 0.3,
                           min_samples_for_estimation = 3,
                           exclude_from_background    = NULL,
                           assay                      = "RNA",
                           min_cells                  = 10,
                           features                   = NULL,
                           protected_genes            = NULL,
                           make_plots                 = TRUE,
                           label_difference           = 0.5,
                           highlight_genes            = c("INS", "TTR", "GCG",
                                                          "PRSS1", "CELA3A"),
                           run_consistency_analysis   = TRUE,
                           min_sample_fraction        = 0.5,
                           min_cell_types             = 2,
                           export_offenders           = TRUE,
                           dataset_name               = "my_dataset",
                           tissue                     = "unknown",
                           n_donors                   = NULL,
                           platform                   = "10x_Chromium",
                           return_full_corrected_obj  = FALSE,
                           save_dir                   = NULL,
                           verbose                    = TRUE) {

  # ── Validate inputs ────────────────────────────────────────────────────────
  .validate_seurat(object, assay)
  .validate_columns(object, sample_col, cell_type_col)

  if (!is.null(save_dir) && !dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
    if (verbose) cli::cli_alert_info(paste0("Created output directory: ", save_dir))
  }

  # ── Resolve cell types to process ─────────────────────────────────────────
  all_types <- unique(as.character(object@meta.data[[cell_type_col]]))
  if (is.null(cell_types)) {
    cell_types <- all_types
    if (verbose) {
      cli::cli_alert_info(
        paste0("cell_types not specified — processing all ",
               length(cell_types), " cell types found in '",
               cell_type_col, "'.")
      )
    }
  } else {
    missing_types <- setdiff(cell_types, all_types)
    if (length(missing_types) > 0) {
      rlang::warn(
        paste0("The following cell_types were not found in '", cell_type_col,
               "': ", paste(missing_types, collapse = ", "),
               ". They will be skipped.")
      )
      cell_types <- intersect(cell_types, all_types)
    }
  }

  if (length(cell_types) == 0) {
    rlang::abort("No valid cell types to process.")
  }

  # ── Step 1: Estimate background ────────────────────────────────────────────
  if (verbose) cli::cli_h1("AmbientFilter — per-cell-type ambient RNA correction")

  background <- estimate_background(
    object             = object,
    sample_col         = sample_col,
    tissue_col         = tissue_col,
    assay              = assay,
    exclude_cell_types = exclude_from_background,
    cell_type_col      = cell_type_col,
    verbose            = verbose
  )

  # ── Step 1.5: Estimate correction strength (if auto) ────────────────────
  # If correction_strength is NULL or auto_estimate_strength = TRUE,
  # learn correction strengths from the data using cross-sample
  # ambient correlation analysis. Manual values can still be supplied
  # via prior_strengths and blended via prior_weight.
  strength_estimates <- NULL

  if (auto_estimate_strength || is.null(correction_strength)) {
    if (verbose) {
      cli::cli_h2("Step 1.5: Estimating data-driven correction strengths")
      cli::cli_alert_info(
        paste0("Learning correction strengths from cross-sample ambient ",
               "correlations. Supply prior_strengths to blend with biological ",
               "knowledge (e.g. macrophage = 0.5).")
      )
    }
    strength_estimates <- estimate_correction_strength(
      object              = object,
      background          = background,
      cell_type_col       = cell_type_col,
      sample_col          = sample_col,
      tissue_col          = tissue_col,
      assay               = assay,
      cell_types          = cell_types,
      min_cells           = min_cells,
      min_samples         = min_samples_for_estimation,
      top_genes_pct       = top_genes_pct,
      compute_gene_weights = TRUE,
      prior_strengths     = prior_strengths,
      prior_weight        = prior_weight,
      make_plots          = make_plots,
      verbose             = verbose
    )
    # Use estimated strengths, overriding any manually supplied value
    correction_strength <- strength_estimates$correction_strength

    if (verbose) {
      cli::cli_alert_success("Correction strengths estimated from data.")
    }
  } else {
    if (verbose) {
      cli::cli_alert_info(
        paste0("Using manually supplied correction_strength. ",
               "Set auto_estimate_strength = TRUE to learn from data instead.")
      )
    }
  }

  # ── Step 2: Correct each cell type ────────────────────────────────────────
  cell_type_results <- vector("list", length(cell_types))
  names(cell_type_results) <- cell_types
  plots <- vector("list", length(cell_types))
  names(plots) <- cell_types

  for (ct in cell_types) {

    # Resolve per-cell-type correction strength
    ct_strength <- .resolve_correction_strength(correction_strength, ct)

    # Resolve per-cell-type feature overrides
    ct_features <- if (!is.null(features)) features[[ct]] else NULL

    # Run correction
    # Resolve protected genes for this cell type
    ct_protected <- NULL
    if (!is.null(protected_genes)) {
      if ("all" %in% names(protected_genes)) {
        ct_protected <- protected_genes[["all"]]
      }
      if (ct %in% names(protected_genes)) {
        ct_protected <- unique(c(ct_protected, protected_genes[[ct]]))
      }
      if (!is.null(ct_protected) && length(ct_protected) > 0 && verbose) {
        cli::cli_alert_info(
          paste0("  Protecting ", length(ct_protected),
                 " genes from correction in '", ct, "'.")
        )
      }
    }

    ct_result <- correct_cell_type(
      object              = object,
      background          = background,
      cell_label          = ct,
      cell_type_col       = cell_type_col,
      sample_col          = sample_col,
      tissue_col          = tissue_col,
      assay               = assay,
      correction_strength = ct_strength,
      min_cells           = min_cells,
      features            = ct_features,
      protected_genes     = ct_protected,
      verbose             = verbose
    )
    cell_type_results[[ct]] <- ct_result

    # Diagnostics
    if (make_plots) {
      ct_plots <- plot_correction_diagnostics(
        result           = ct_result,
        cell_label       = ct,
        assay            = assay,
        label_difference = label_difference,
        highlight_genes  = highlight_genes
      )
      plots[[ct]] <- ct_plots

      # Save if requested
      if (!is.null(save_dir)) {
        ct_dir <- file.path(save_dir, ct)
        dir.create(ct_dir, showWarnings = FALSE)
        ggplot2::ggsave(
          file.path(ct_dir, "correlation.pdf"),
          ct_plots$correlation, width = 8, height = 6
        )
        if (!is.null(ct_plots$top_genes_bar)) {
          ggplot2::ggsave(
            file.path(ct_dir, "top_changed_genes.pdf"),
            ct_plots$top_genes_bar, width = 7, height = 5
          )
        }
        if (!is.null(ct_plots$violin_original)) {
          ggplot2::ggsave(
            file.path(ct_dir, "violin_original.pdf"),
            ct_plots$violin_original, width = 10, height = 5
          )
        }
        if (!is.null(ct_plots$violin_corrected)) {
          ggplot2::ggsave(
            file.path(ct_dir, "violin_corrected.pdf"),
            ct_plots$violin_corrected, width = 10, height = 5
          )
        }
        # Save changed genes table
        utils::write.csv(
          ct_result$changed_genes,
          file.path(ct_dir, "changed_genes.csv"),
          row.names = FALSE
        )
      }
    }
  }

  # ── Step 3: Write corrected counts back to full object ────────────────────
  updated_object <- NULL
  if (return_full_corrected_obj) {
    if (verbose) {
      cli::cli_h2("Writing corrected counts back to full object")
      cli::cli_alert_info(
        "Note: this duplicates the Seurat object in memory."
      )
    }
    updated_object <- object
    for (ct in cell_types) {
      ct_cells     <- rownames(object@meta.data)[
        object@meta.data[[cell_type_col]] == ct
      ]
      corr_counts  <- .get_counts(cell_type_results[[ct]]$corrected, assay)
      # Match barcodes (subset may reorder)
      shared_cells <- intersect(ct_cells, colnames(corr_counts))
      if (length(shared_cells) > 0) {
        SeuratObject::LayerData(
          updated_object, assay = assay, layer = "counts"
        )[, shared_cells] <- corr_counts[, shared_cells]
      }
    }
    updated_object <- Seurat::NormalizeData(
      updated_object, assay = assay, verbose = FALSE
    )
    if (verbose) cli::cli_alert_success("Full corrected object ready.")
  }

  # ── Step 4: Build summary ─────────────────────────────────────────────────
  result_obj <- list(
    cell_type_results  = cell_type_results,
    background         = background,
    strength_estimates = strength_estimates,
    plots              = if (make_plots) plots else NULL,
    updated_object     = updated_object,
    ambient_map        = NULL,   # populated below if run_consistency_analysis
    params = list(
      cell_type_col              = cell_type_col,
      sample_col                 = sample_col,
      tissue_col                 = tissue_col,
      cell_types                 = cell_types,
      correction_strength        = correction_strength,
      auto_estimate_strength     = auto_estimate_strength,
      prior_strengths            = prior_strengths,
      prior_weight               = prior_weight,
      protected_genes            = protected_genes,
      exclude_from_background    = exclude_from_background,
      assay                      = assay,
      min_cells                  = min_cells,
      run_consistency_analysis   = run_consistency_analysis,
      min_sample_fraction        = min_sample_fraction,
      min_cell_types             = min_cell_types,
      date_run                   = Sys.time()
    )
  )
  result_obj$summary <- summarise_correction(result_obj)

  # Save summary
  if (!is.null(save_dir)) {
    utils::write.csv(
      result_obj$summary,
      file.path(save_dir, "correction_summary.csv"),
      row.names = FALSE
    )
  }

  # ── Step 5: Cross-sample consistency analysis (hybrid method) ─────────────
  if (run_consistency_analysis) {
    if (verbose) cli::cli_h1("Step 5: Cross-sample ambient consistency analysis")

    consistency_save_dir <- if (!is.null(save_dir)) {
      file.path(save_dir, "ambient_offenders")
    } else NULL

    ambient_map <- identify_consistent_ambient(
      results             = result_obj,
      object              = object,
      sample_col          = sample_col,
      tissue_col          = tissue_col,
      cell_type_col       = cell_type_col,
      assay               = assay,
      lfc_threshold       = 0.3,
      min_sample_fraction = min_sample_fraction,
      min_cell_types      = min_cell_types,
      min_cells           = min_cells,
      make_plots          = make_plots,
      save_dir            = consistency_save_dir,
      verbose             = verbose
    )
    result_obj$ambient_map <- ambient_map

    # Export shareable offenders table
    if (export_offenders && !is.null(save_dir)) {
      export_ambient_offenders(
        ambient_map  = ambient_map,
        dataset_name = dataset_name,
        tissue       = tissue,
        n_donors     = n_donors,
        platform     = platform,
        notes        = paste0("Generated by AmbientFilter::ambient_filter(). ",
                              "Cell types corrected: ",
                              paste(cell_types, collapse = ", ")),
        save_dir     = consistency_save_dir,
        save_excel   = TRUE
      )
    }
  }

  if (verbose) {
    cli::cli_h1("AmbientFilter complete")
    cli::cli_h2("Per-cell-type correction summary")
    print(result_obj$summary)
    if (run_consistency_analysis && !is.null(result_obj$ambient_map)) {
      n_cross <- nrow(result_obj$ambient_map$cross_cell_type)
      cli::cli_alert_success(
        paste0(n_cross, " cross-cell-type ambient offenders identified. ",
               "See results$ambient_map and results$ambient_map$plots.")
      )
      if (n_cross > 0) {
        cli::cli_h2("Top cross-cell-type ambient offenders")
        print(head(result_obj$ambient_map$cross_cell_type, 10))
      }
    }
  }

  result_obj
}
