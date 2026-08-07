# ── AmbientFilter — internal utilities ────────────────────────────────────────
# ────────────────────────────────────────────────────────────────────────────────

#' @keywords internal
.validate_seurat <- function(object, assay) {
  if (!methods::is(object, "Seurat")) {
    rlang::abort("`object` must be a Seurat object.")
  }
  if (!assay %in% SeuratObject::Assays(object)) {
    rlang::abort(
      paste0("Assay '", assay, "' not found in object. ",
             "Available assays: ", paste(SeuratObject::Assays(object), collapse = ", "))
    )
  }
  # Seurat v5: check counts layer exists
  counts_check <- tryCatch(
    SeuratObject::LayerData(object, assay = assay, layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(counts_check)) {
    rlang::abort(
      paste0("No 'counts' layer found in assay '", assay, "'. ",
             "Run JoinLayers() or ensure raw counts are present.")
    )
  }
  invisible(TRUE)
}

#' @keywords internal
.validate_columns <- function(object, sample_col, cell_type_col) {
  meta <- object@meta.data
  if (!sample_col %in% colnames(meta)) {
    rlang::abort(
      paste0("sample_col '", sample_col, "' not found in metadata. ",
             "Available columns: ", paste(colnames(meta), collapse = ", "))
    )
  }
  if (!cell_type_col %in% colnames(meta)) {
    rlang::abort(
      paste0("cell_type_col '", cell_type_col, "' not found in metadata. ",
             "Available columns: ", paste(colnames(meta), collapse = ", "))
    )
  }
  invisible(TRUE)
}

#' @keywords internal
.get_counts <- function(object, assay) {
  SeuratObject::LayerData(object, assay = assay, layer = "counts")
}

#' @keywords internal
.set_counts <- function(object, assay, counts_mat) {
  SeuratObject::LayerData(object, assay = assay, layer = "counts") <- counts_mat
  # Re-normalise data layer
  object <- Seurat::NormalizeData(object, assay = assay, verbose = FALSE)
  object
}

#' @keywords internal
.get_cell_barcodes <- function(object, cell_type_col, cell_label,
                               sample_col, sample_id) {
  meta <- object@meta.data
  idx <- meta[[cell_type_col]] == cell_label &
         meta[[sample_col]]   == sample_id
  rownames(meta)[idx]
}

#' @keywords internal
.resolve_correction_strength <- function(correction_strength, cell_label) {
  if (is.null(correction_strength)) return(1.0)
  if (is.numeric(correction_strength) && length(correction_strength) == 1) {
    return(correction_strength)
  }
  if (is.list(correction_strength) || is.numeric(correction_strength)) {
    val <- correction_strength[[cell_label]]
    if (is.null(val)) {
      rlang::warn(
        paste0("No correction_strength entry for '", cell_label,
               "'. Defaulting to 1.0 (full correction).")
      )
      return(1.0)
    }
    if (val < 0 || val > 1) {
      rlang::abort(
        paste0("correction_strength values must be between 0 and 1. ",
               "Got ", val, " for '", cell_label, "'.")
      )
    }
    return(val)
  }
  rlang::abort(
    "`correction_strength` must be NULL, a single numeric, or a named list/vector."
  )
}

#' @keywords internal
.check_min_cells <- function(n_cells, min_cells, cell_label, sample_id) {
  if (n_cells < min_cells) {
    rlang::warn(
      paste0("Cell type '", cell_label, "' in sample '", sample_id,
             "' has only ", n_cells, " cells (min_cells = ", min_cells,
             "). Skipping this sample/cell-type combination.")
    )
    return(FALSE)
  }
  TRUE
}

#' @keywords internal
.make_corrected_seurat <- function(original_subset, corrected_counts,
                                   assay, cell_label) {
  corrected <- original_subset
  SeuratObject::LayerData(corrected, assay = assay, layer = "counts") <- corrected_counts
  corrected <- Seurat::NormalizeData(corrected, assay = assay, verbose = FALSE)
  corrected
}
