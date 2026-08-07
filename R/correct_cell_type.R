# ── AmbientFilter — per-cell-type correction ──────────────────────────────────
# ────────────────────────────────────────────────────────────────────────────────

#' Correct ambient RNA contamination for a single cell type
#'
#' Subsets cells of the specified type, groups them across all samples
#' (as in dRopt — grouping avoids over-correcting cells that happen to be
#' near a particular cell type), and applies [DropletUtils::removeAmbience()]
#' using the pre-computed per-sample background profiles.
#'
#' @param object A Seurat object.
#' @param background A named list of background profiles as returned by
#'   [estimate_background()]. Names must match the values in `sample_col`
#'   (or `"<tissue>.<sample>"` if `tissue_col` was used).
#' @param cell_label Character. The cell type label to correct. Must be a
#'   value present in `cell_type_col`.
#' @param cell_type_col Character. Metadata column containing cell type labels.
#'   Default `"cell_type"`.
#' @param sample_col Character. Metadata column containing sample IDs.
#'   Default `"orig.ident"`.
#' @param tissue_col Character or NULL. Tissue metadata column. Must match
#'   what was used in [estimate_background()]. Default `NULL`.
#' @param assay Character. Assay to correct. Default `"RNA"`.
#' @param correction_strength Numeric in \[0, 1\]. Scaling factor applied to
#'   the ambient subtraction. `1.0` = full correction (appropriate for
#'   lymphocytes, ductal cells, etc.). `0.5` = half correction (recommended
#'   for macrophages, which legitimately internalise ambient RNA via
#'   efferocytosis/trogocytosis). `0` = no correction. Can be a single value
#'   applied to all samples or a named numeric vector of per-sample values.
#'   Default `1.0`.
#' @param min_cells Integer. Minimum number of cells required for a
#'   cell-type/sample combination to be included. Combinations below this
#'   threshold are skipped and a warning is raised. Default `10`.
#' @param features Character vector or NULL. Specific genes to force-correct
#'   (passed to `features` in `removeAmbience`). Use only when known
#'   non-expressed genes are not being removed automatically. Default `NULL`.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A list with elements:
#' \describe{
#'   \item{`original`}{Seurat object subset containing only `cell_label` cells,
#'     with uncorrected counts.}
#'   \item{`corrected`}{Seurat object subset with ambient-corrected counts and
#'     re-normalised data layer.}
#'   \item{`changed_genes`}{Data frame of genes whose mean expression changed
#'     by more than `log2(1.5)` between original and corrected, sorted by
#'     magnitude of change.}
#'   \item{`samples_used`}{Character vector of sample IDs that had sufficient
#'     cells and were included.}
#'   \item{`samples_skipped`}{Character vector of sample IDs that were skipped
#'     due to insufficient cells.}
#' }
#'
#' @examples
#' \dontrun{
#' bg <- estimate_background(seu, sample_col = "hpap_id", tissue_col = "tissue")
#'
#' lymphocyte_result <- correct_cell_type(
#'   object             = seu,
#'   background         = bg,
#'   cell_label         = "lymphocyte",
#'   cell_type_col      = "cell_type",
#'   sample_col         = "hpap_id",
#'   tissue_col         = "tissue",
#'   correction_strength = 1.0
#' )
#'
#' macrophage_result <- correct_cell_type(
#'   object             = seu,
#'   background         = bg,
#'   cell_label         = "macrophage",
#'   correction_strength = 0.5  # conservative — macrophages are nibbly
#' )
#' }
#'
#' @export
correct_cell_type <- function(object,
                              background,
                              cell_label,
                              cell_type_col       = "cell_type",
                              sample_col          = "orig.ident",
                              tissue_col          = NULL,
                              assay               = "RNA",
                              correction_strength = 1.0,
                              min_cells           = 10,
                              features            = NULL,
                              verbose             = TRUE) {

  .validate_seurat(object, assay)

  meta   <- object@meta.data
  counts <- .get_counts(object, assay)

  # ── Subset to target cell type ─────────────────────────────────────────────
  target_cells <- rownames(meta)[meta[[cell_type_col]] == cell_label]
  if (length(target_cells) == 0) {
    rlang::abort(
      paste0("No cells found with cell_type_col '", cell_type_col,
             "' == '", cell_label, "'.")
    )
  }

  original_subset <- subset(object, cells = target_cells)
  meta_sub        <- original_subset@meta.data
  counts_sub      <- .get_counts(original_subset, assay)

  # ── Build sample key matching background list names ────────────────────────
  if (!is.null(tissue_col)) {
    sample_key <- paste(meta_sub[[tissue_col]], meta_sub[[sample_col]], sep = ".")
  } else {
    sample_key <- meta_sub[[sample_col]]
  }

  samples_all     <- unique(sample_key)
  samples_used    <- character(0)
  samples_skipped <- character(0)

  # Corrected counts matrix — start as copy, update per sample
  corrected_counts <- counts_sub

  if (verbose) {
    cli::cli_h2(paste0("Correcting: ", cell_label))
    cli::cli_alert_info(
      paste0(length(target_cells), " cells across ",
             length(samples_all), " samples | correction_strength = ",
             correction_strength)
    )
  }

  pb <- if (verbose) {
    cli::cli_progress_bar(
      paste0("  Samples [", cell_label, "]"),
      total = length(samples_all)
    )
  } else NULL

  for (s in samples_all) {

    # Check background exists for this sample
    if (!s %in% names(background)) {
      rlang::warn(
        paste0("No background profile found for sample '", s,
               "'. Skipping.")
      )
      samples_skipped <- c(samples_skipped, s)
      if (verbose) cli::cli_progress_update(id = pb)
      next
    }

    cell_idx <- which(sample_key == s)
    n_cells  <- length(cell_idx)

    if (!.check_min_cells(n_cells, min_cells, cell_label, s)) {
      samples_skipped <- c(samples_skipped, s)
      if (verbose) cli::cli_progress_update(id = pb)
      next
    }

    # Extract counts for this cell type in this sample
    sample_counts <- counts_sub[, cell_idx, drop = FALSE]
    ambient_prof  <- background[[s]]

    # Align genes (background may have different order / superset)
    shared_genes  <- intersect(rownames(sample_counts), names(ambient_prof))
    if (length(shared_genes) == 0) {
      rlang::warn(
        paste0("No shared genes between count matrix and background for '",
               s, "'. Skipping.")
      )
      samples_skipped <- c(samples_skipped, s)
      if (verbose) cli::cli_progress_update(id = pb)
      next
    }
    ambient_aligned <- ambient_prof[shared_genes]

    # Scale by correction_strength (0 = no correction, 1 = full)
    ambient_scaled <- ambient_aligned * correction_strength

    # ── Call DropletUtils::removeAmbience ────────────────────────────────────
    # removeAmbience expects:
    #   y         = integer count matrix (genes x cells)
    #   ambient   = named numeric vector (gene proportions)
    #   features  = optional character vector of genes to force-correct
    corrected_sample <- tryCatch({
      DropletUtils::removeAmbience(
        y        = sample_counts[shared_genes, , drop = FALSE],
        ambient  = ambient_scaled,
        features = features,
        verbose  = FALSE
      )
    }, error = function(e) {
      rlang::warn(
        paste0("removeAmbience failed for '", cell_label, "' sample '", s,
               "': ", conditionMessage(e), ". Returning uncorrected counts.")
      )
      NULL
    })

    if (!is.null(corrected_sample)) {
      # Write corrected counts back into the full corrected matrix
      corrected_counts[shared_genes, cell_idx] <- corrected_sample
      samples_used <- c(samples_used, s)
    } else {
      samples_skipped <- c(samples_skipped, s)
    }

    if (verbose) cli::cli_progress_update(id = pb)
  }

  if (verbose) {
    cli::cli_progress_done(id = pb)
    cli::cli_alert_success(
      paste0("  Corrected ", length(samples_used), "/", length(samples_all),
             " samples.")
    )
  }

  # ── Build corrected Seurat subset ──────────────────────────────────────────
  corrected_subset <- .make_corrected_seurat(original_subset,
                                             corrected_counts,
                                             assay, cell_label)

  # ── Summarise changed genes ────────────────────────────────────────────────
  orig_means <- Matrix::rowMeans(
    .get_counts(original_subset, assay)
  )
  corr_means <- Matrix::rowMeans(corrected_counts)

  # log2 fold change of corrected vs original mean counts
  # add pseudocount to avoid log(0)
  lfc <- log2((corr_means + 1) / (orig_means + 1))
  changed <- data.frame(
    gene            = names(lfc),
    mean_original   = orig_means,
    mean_corrected  = corr_means,
    log2fc          = lfc,
    abs_log2fc      = abs(lfc),
    direction       = ifelse(lfc < 0, "decreased", "unchanged"),
    stringsAsFactors = FALSE
  )
  changed <- changed[order(changed$abs_log2fc, decreasing = TRUE), ]
  changed <- changed[changed$abs_log2fc > log2(1.5), ]
  rownames(changed) <- NULL

  list(
    original        = original_subset,
    corrected       = corrected_subset,
    changed_genes   = changed,
    samples_used    = samples_used,
    samples_skipped = samples_skipped
  )
}
