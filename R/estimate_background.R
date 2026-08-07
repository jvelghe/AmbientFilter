# ── AmbientFilter — background estimation ─────────────────────────────────────
# ────────────────────────────────────────────────────────────────────────────────

#' Estimate per-sample ambient RNA background profiles
#'
#' For each sample (optionally stratified by tissue), computes the proportion
#' of total counts attributable to each gene. This profile represents the
#' expected ambient contamination distribution and is passed to
#' [DropletUtils::removeAmbience()] during correction.
#'
#' The core assumption (inherited from dRopt) is that the ambient pool reflects
#' the transcriptional composition of all cells present in the sample. In islet
#' datasets this means INS and TTR dominate the ambient profile, which is
#' exactly the signal to subtract from non-expressing cell types.
#'
#' @param object A Seurat object containing all cell types. Must have a raw
#'   counts layer.
#' @param sample_col Character. Metadata column containing sample/donor IDs.
#'   Default `"orig.ident"`.
#' @param tissue_col Character or NULL. Metadata column containing tissue
#'   identity (e.g. `"tissue"` with values `"pancreas"`, `"spleen"`). When
#'   provided, background is estimated per sample *within* tissue to avoid
#'   cross-tissue contamination of the ambient profile. Strongly recommended
#'   when your object contains cells from multiple tissues. Default `NULL`.
#' @param assay Character. Assay to use for count extraction. Default `"RNA"`.
#' @param exclude_cell_types Character vector or NULL. Cell type labels
#'   (from `cell_type_col`) to exclude when computing the background. Useful
#'   if a cell type is so dominant in a sample (e.g. beta cells in islet
#'   enrichments) that it artificially inflates the background. Default `NULL`
#'   (use all cells).
#' @param cell_type_col Character. Metadata column containing cell type labels.
#'   Required only if `exclude_cell_types` is specified. Default `"cell_type"`.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#'
#' @return A named list of named numeric vectors. The outer name is the
#'   sample ID (or `"<tissue>.<sample>"` when `tissue_col` is supplied).
#'   Each inner vector gives the proportion of reads per gene for that sample.
#'   Proportions sum to 1. Genes with zero counts across the sample are
#'   included with proportion 0.
#'
#' @examples
#' \dontrun{
#' bg <- estimate_background(
#'   object      = seu,
#'   sample_col  = "hpap_id",
#'   tissue_col  = "tissue",
#'   assay       = "RNA"
#' )
#'
#' # Inspect the dominant ambient genes for one sample
#' head(sort(bg[["pancreas.HPAP022"]], decreasing = TRUE), 20)
#' }
#'
#' @export
estimate_background <- function(object,
                                sample_col   = "orig.ident",
                                tissue_col   = NULL,
                                assay        = "RNA",
                                exclude_cell_types = NULL,
                                cell_type_col      = "cell_type",
                                verbose      = TRUE) {

  .validate_seurat(object, assay)
  .validate_columns(object, sample_col,
                    if (!is.null(exclude_cell_types)) cell_type_col else sample_col)

  meta   <- object@meta.data
  counts <- .get_counts(object, assay)

  # ── Build grouping key ──────────────────────────────────────────────────────
  if (!is.null(tissue_col)) {
    if (!tissue_col %in% colnames(meta)) {
      rlang::abort(
        paste0("tissue_col '", tissue_col, "' not found in metadata.")
      )
    }
    group_key <- paste(meta[[tissue_col]], meta[[sample_col]], sep = ".")
  } else {
    group_key <- meta[[sample_col]]
  }

  samples <- unique(group_key)

  if (verbose) {
    cli::cli_h2("Estimating background profiles")
    cli::cli_alert_info(
      paste0("Found ", length(samples), " sample",
             if (length(samples) != 1) "s", " to process",
             if (!is.null(tissue_col))
               paste0(" (stratified by '", tissue_col, "')") else "")
    )
  }

  # ── Compute per-sample ambient profile ────────────────────────────────────
  bg_list <- vector("list", length(samples))
  names(bg_list) <- samples

  pb <- if (verbose) {
    cli::cli_progress_bar("Estimating background", total = length(samples))
  } else NULL

  for (s in samples) {
    cells_in_sample <- which(group_key == s)

    # Optionally exclude dominant cell types from background
    if (!is.null(exclude_cell_types) && cell_type_col %in% colnames(meta)) {
      keep <- !meta[[cell_type_col]][cells_in_sample] %in% exclude_cell_types
      cells_in_sample <- cells_in_sample[keep]
      if (length(cells_in_sample) == 0) {
        rlang::warn(
          paste0("After excluding specified cell types, sample '", s,
                 "' has no cells left. Using all cells for background.")
        )
        cells_in_sample <- which(group_key == s)
      }
    }

    # Sum counts across all cells in this sample group
    sample_counts <- counts[, cells_in_sample, drop = FALSE]
    gene_totals   <- Matrix::rowSums(sample_counts)
    total_counts  <- sum(gene_totals)

    if (total_counts == 0) {
      rlang::warn(
        paste0("Sample '", s, "' has zero total counts. ",
               "Setting all proportions to 0.")
      )
      bg_list[[s]] <- setNames(rep(0, nrow(counts)), rownames(counts))
    } else {
      bg_list[[s]] <- as.numeric(gene_totals) / total_counts
      names(bg_list[[s]]) <- rownames(counts)
    }

    if (verbose) cli::cli_progress_update(id = pb)
  }

  if (verbose) {
    cli::cli_progress_done(id = pb)
    cli::cli_alert_success("Background profiles estimated for all samples.")
  }

  bg_list
}
