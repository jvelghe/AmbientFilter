# ── AmbientFilter — cross-sample ambient consistency analysis ──────────────────
# ────────────────────────────────────────────────────────────────────────────────

#' Identify consistently ambient genes across samples and cell types
#'
#' After per-sample ambient correction with [ambient_filter()], this function
#' performs a pooled cross-sample analysis to identify genes that are
#' *consistently* corrected across donors — the signature of true ambient
#' contamination rather than sample-specific processing artefacts.
#'
#' The logic follows from the core tension between per-sample and pooled
#' correction: per-sample correction captures sample-specific ambient soup
#' composition (which varies by dissociation quality, cell viability, and
#' processing batch), while pooled analysis provides the statistical power
#' to distinguish genes that are *consistently* ambient across most donors
#' from those that are only noisy in one or two samples.
#'
#' For each cell type and each sample, this function records which genes
#' were significantly reduced by correction (log2FC < -`lfc_threshold`,
#' the "corrected direction"). A gene is classified as a **consistent ambient
#' offender** in a given cell type if it appears in the corrected gene list
#' in at least `min_sample_fraction` of samples. Genes that are consistently
#' ambient across multiple cell types are flagged as **cross-cell-type offenders**
#' — the strongest evidence for true soup contamination.
#'
#' @param results The list returned by [ambient_filter()].
#' @param object The original Seurat object (needed for per-sample per-cell-type
#'   count extraction).
#' @param sample_col Character. Metadata column containing sample IDs.
#'   Default `"orig.ident"`.
#' @param tissue_col Character or NULL. Tissue metadata column. Default `NULL`.
#' @param cell_type_col Character. Cell type metadata column.
#'   Default `"cell_type"`.
#' @param assay Character. Assay to use. Default `"RNA"`.
#' @param lfc_threshold Numeric. Minimum |log2FC| (corrected vs original) to
#'   consider a gene as corrected in a given sample. Default `0.3`.
#' @param min_sample_fraction Numeric in (0, 1\]. Minimum fraction of samples
#'   in which a gene must be corrected to be classified as a consistent ambient
#'   offender for that cell type. Default `0.5` (present in ≥ 50% of samples).
#' @param min_cell_types Integer. Minimum number of cell types in which a gene
#'   must be a consistent offender to be flagged as a cross-cell-type ambient
#'   gene. Default `2`.
#' @param min_cells Integer. Minimum cells per cell-type/sample combination to
#'   include in the per-sample analysis. Default `5`.
#' @param make_plots Logical. Generate heatmap and summary plots. Default `TRUE`.
#' @param save_dir Character or NULL. Directory to save outputs. Default `NULL`.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A list with:
#' \describe{
#'   \item{`per_cell_type`}{Named list (one per cell type). Each element is a
#'     data frame with columns `gene`, `n_samples_corrected`,
#'     `n_samples_total`, `fraction_corrected`, `mean_lfc_corrected`,
#'     `is_consistent_ambient`. Sorted by `fraction_corrected` descending.}
#'   \item{`cross_cell_type`}{Data frame of genes that are consistent ambient
#'     offenders in ≥ `min_cell_types` cell types. Columns: `gene`,
#'     `n_cell_types_affected`, `cell_types_affected`, `mean_fraction_corrected`.
#'     This is your publishable ambient contamination gene list.}
#'   \item{`consistency_matrix`}{Matrix (genes × cell types) of
#'     `fraction_corrected` values. Suitable for heatmap visualisation.
#'     Only includes genes that are consistent offenders in ≥ 1 cell type.}
#'   \item{`recommended_features`}{Named list of character vectors, one per
#'     cell type, containing the consistent ambient genes for that cell type.
#'     Can be passed directly to the `features` argument of a second-pass
#'     [ambient_filter()] run to force-correct these genes.}
#'   \item{`plots`}{List of ggplot objects (if `make_plots = TRUE`):
#'     `heatmap`, `cross_cell_type_bar`, `per_sample_consistency`.}
#' }
#'
#' @section Two-pass workflow:
#' The recommended workflow is:
#' 1. Run [ambient_filter()] with `correction_strength` tuned per cell type.
#' 2. Run [identify_consistent_ambient()] on the results.
#' 3. Inspect `cross_cell_type` — these are your true ambient offenders.
#' 4. Optionally re-run [ambient_filter()] with `features = recommended_features`
#'    to force-correct these genes even in samples where they weren't
#'    automatically removed.
#'
#' @examples
#' \dontrun{
#' # Step 1: run per-sample correction
#' results <- ambient_filter(seu, ...)
#'
#' # Step 2: identify consistent ambient offenders across samples
#' ambient_map <- identify_consistent_ambient(
#'   results             = results,
#'   object              = seu,
#'   sample_col          = "hpap_id",
#'   tissue_col          = "tissue",
#'   min_sample_fraction = 0.5,
#'   min_cell_types      = 2
#' )
#'
#' # Inspect cross-cell-type offenders
#' ambient_map$cross_cell_type
#'
#' # View the ambient contamination heatmap
#' ambient_map$plots$heatmap
#'
#' # Step 3 (optional): re-run with force-correction of consistent offenders
#' results_v2 <- ambient_filter(
#'   seu,
#'   features = ambient_map$recommended_features,
#'   ...
#' )
#' }
#'
#' @export
identify_consistent_ambient <- function(results,
                                        object,
                                        sample_col          = "orig.ident",
                                        tissue_col          = NULL,
                                        cell_type_col       = "cell_type",
                                        assay               = "RNA",
                                        lfc_threshold       = 0.3,
                                        min_sample_fraction = 0.5,
                                        min_cell_types      = 2,
                                        min_cells           = 5,
                                        make_plots          = TRUE,
                                        save_dir            = NULL,
                                        verbose             = TRUE) {

  if (!is.null(save_dir) && !dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  cell_types <- names(results$cell_type_results)
  meta       <- object@meta.data

  # Build sample key (matching estimate_background)
  if (!is.null(tissue_col) && tissue_col %in% colnames(meta)) {
    sample_key <- paste(meta[[tissue_col]], meta[[sample_col]], sep = ".")
  } else {
    sample_key <- meta[[sample_col]]
  }

  if (verbose) {
    cli::cli_h1("Identifying consistent ambient genes across samples")
    cli::cli_alert_info(
      paste0("Analysing ", length(cell_types), " cell types | ",
             "min_sample_fraction = ", min_sample_fraction, " | ",
             "min_cell_types = ", min_cell_types)
    )
  }

  # ── Per-cell-type per-sample gene-level analysis ──────────────────────────
  per_ct_results <- vector("list", length(cell_types))
  names(per_ct_results) <- cell_types

  pb <- if (verbose) {
    cli::cli_progress_bar("Analysing cell types", total = length(cell_types))
  } else NULL

  for (ct in cell_types) {

    ct_res  <- results$cell_type_results[[ct]]
    orig    <- ct_res$original
    corr    <- ct_res$corrected

    # Get the sample key for cells in this subset
    ct_cells     <- colnames(orig)
    ct_meta      <- meta[ct_cells, , drop = FALSE]
    if (!is.null(tissue_col) && tissue_col %in% colnames(ct_meta)) {
      ct_sample_key <- paste(ct_meta[[tissue_col]],
                             ct_meta[[sample_col]], sep = ".")
    } else {
      ct_sample_key <- ct_meta[[sample_col]]
    }

    samples_in_ct <- unique(ct_sample_key)

    # Per-sample LFC for each gene
    # Collect: for each sample, which genes were reduced (lfc < -threshold)?
    sample_gene_lfc <- vector("list", length(samples_in_ct))
    names(sample_gene_lfc) <- samples_in_ct

    for (s in samples_in_ct) {
      s_cells <- which(ct_sample_key == s)
      if (length(s_cells) < min_cells) next

      orig_counts <- .get_counts(orig, assay)[, s_cells, drop = FALSE]
      corr_counts <- .get_counts(corr, assay)[, s_cells, drop = FALSE]

      orig_means <- Matrix::rowMeans(orig_counts)
      corr_means <- Matrix::rowMeans(corr_counts)

      lfc <- log2((corr_means + 1) / (orig_means + 1))
      sample_gene_lfc[[s]] <- lfc
    }

    # Remove NULLs (samples with too few cells)
    sample_gene_lfc <- Filter(Negate(is.null), sample_gene_lfc)
    n_samples        <- length(sample_gene_lfc)

    if (n_samples == 0) {
      per_ct_results[[ct]] <- data.frame(
        gene                 = character(0),
        n_samples_corrected  = integer(0),
        n_samples_total      = integer(0),
        fraction_corrected   = numeric(0),
        mean_lfc_corrected   = numeric(0),
        is_consistent_ambient = logical(0),
        stringsAsFactors     = FALSE
      )
      if (verbose) cli::cli_progress_update(id = pb)
      next
    }

    # All genes present across samples
    all_genes <- Reduce(union, lapply(sample_gene_lfc, names))

    # For each gene: in how many samples was it reduced below -lfc_threshold?
    gene_summary <- lapply(all_genes, function(g) {
      lfcs_for_gene <- sapply(sample_gene_lfc, function(s_lfc) {
        if (g %in% names(s_lfc)) s_lfc[[g]] else NA_real_
      })
      valid      <- lfcs_for_gene[!is.na(lfcs_for_gene)]
      n_corrected <- sum(valid < -lfc_threshold)
      data.frame(
        gene                  = g,
        n_samples_corrected   = n_corrected,
        n_samples_total       = length(valid),
        fraction_corrected    = if (length(valid) > 0) n_corrected / length(valid) else 0,
        mean_lfc_corrected    = if (length(valid) > 0) mean(valid) else NA_real_,
        stringsAsFactors      = FALSE
      )
    })
    gene_df <- do.call(rbind, gene_summary)
    gene_df$is_consistent_ambient <- gene_df$fraction_corrected >= min_sample_fraction
    gene_df <- gene_df[order(gene_df$fraction_corrected, decreasing = TRUE), ]
    rownames(gene_df) <- NULL

    per_ct_results[[ct]] <- gene_df
    if (verbose) cli::cli_progress_update(id = pb)
  }

  if (verbose) cli::cli_progress_done(id = pb)

  # ── Cross-cell-type offender analysis ────────────────────────────────────
  if (verbose) cli::cli_h2("Computing cross-cell-type ambient offenders")

  # Collect consistent ambient genes per cell type
  consistent_per_ct <- lapply(per_ct_results, function(df) {
    if (nrow(df) == 0) return(character(0))
    df$gene[df$is_consistent_ambient]
  })

  # Which genes appear as consistent offenders in >= min_cell_types cell types?
  all_consistent <- unlist(consistent_per_ct)
  gene_ct_counts <- table(all_consistent)

  cross_ct_genes <- names(gene_ct_counts)[gene_ct_counts >= min_cell_types]

  cross_ct_df <- do.call(rbind, lapply(cross_ct_genes, function(g) {
    cts_affected <- names(consistent_per_ct)[
      sapply(consistent_per_ct, function(genes) g %in% genes)
    ]
    # Mean fraction corrected across affected cell types
    mean_frac <- mean(sapply(cts_affected, function(ct) {
      df <- per_ct_results[[ct]]
      df$fraction_corrected[df$gene == g]
    }), na.rm = TRUE)

    data.frame(
      gene                  = g,
      n_cell_types_affected = length(cts_affected),
      cell_types_affected   = paste(cts_affected, collapse = ", "),
      mean_fraction_corrected = round(mean_frac, 3),
      stringsAsFactors      = FALSE
    )
  }))

  if (!is.null(cross_ct_df) && nrow(cross_ct_df) > 0) {
    cross_ct_df <- cross_ct_df[order(cross_ct_df$n_cell_types_affected,
                                     cross_ct_df$mean_fraction_corrected,
                                     decreasing = TRUE), ]
    rownames(cross_ct_df) <- NULL
  } else {
    cross_ct_df <- data.frame(
      gene = character(0), n_cell_types_affected = integer(0),
      cell_types_affected = character(0),
      mean_fraction_corrected = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  if (verbose) {
    cli::cli_alert_success(
      paste0("Found ", nrow(cross_ct_df),
             " cross-cell-type ambient offenders (present in >= ",
             min_cell_types, " cell types).")
    )
  }

  # ── Build consistency matrix (genes x cell types) ─────────────────────────
  # Only genes that are consistent in >= 1 cell type
  all_consistent_unique <- unique(unlist(consistent_per_ct))

  if (length(all_consistent_unique) > 0) {
    cons_mat <- do.call(cbind, lapply(cell_types, function(ct) {
      df <- per_ct_results[[ct]]
      frac <- setNames(df$fraction_corrected, df$gene)
      frac[all_consistent_unique]
    }))
    colnames(cons_mat)  <- cell_types
    rownames(cons_mat)  <- all_consistent_unique
    cons_mat[is.na(cons_mat)] <- 0
    # Sort rows by mean across cell types
    cons_mat <- cons_mat[order(rowMeans(cons_mat), decreasing = TRUE), ,
                         drop = FALSE]
  } else {
    cons_mat <- matrix(nrow = 0, ncol = length(cell_types),
                       dimnames = list(NULL, cell_types))
  }

  # ── Recommended features for second-pass correction ───────────────────────
  recommended_features <- lapply(per_ct_results, function(df) {
    if (nrow(df) == 0) return(character(0))
    df$gene[df$is_consistent_ambient]
  })

  # ── Plots ─────────────────────────────────────────────────────────────────
  plots <- NULL
  if (make_plots) {
    plots <- .make_ambient_plots(
      per_ct_results    = per_ct_results,
      cross_ct_df       = cross_ct_df,
      cons_mat          = cons_mat,
      cell_types        = cell_types,
      min_sample_fraction = min_sample_fraction,
      save_dir          = save_dir
    )
  }

  # ── Save outputs ──────────────────────────────────────────────────────────
  if (!is.null(save_dir)) {
    utils::write.csv(
      cross_ct_df,
      file.path(save_dir, "cross_celltype_ambient_genes.csv"),
      row.names = FALSE
    )
    for (ct in cell_types) {
      utils::write.csv(
        per_ct_results[[ct]],
        file.path(save_dir, paste0(ct, "_ambient_consistency.csv")),
        row.names = FALSE
      )
    }
    if (nrow(cons_mat) > 0) {
      utils::write.csv(
        as.data.frame(cons_mat),
        file.path(save_dir, "consistency_matrix.csv")
      )
    }
    if (verbose) {
      cli::cli_alert_success(paste0("Results saved to: ", save_dir))
    }
  }

  if (verbose) {
    cli::cli_h2("Cross-cell-type ambient offenders")
    if (nrow(cross_ct_df) > 0) {
      print(head(cross_ct_df, 20))
    } else {
      cli::cli_alert_info("No cross-cell-type offenders found at current thresholds.")
    }
  }

  list(
    per_cell_type         = per_ct_results,
    cross_cell_type       = cross_ct_df,
    consistency_matrix    = cons_mat,
    recommended_features  = recommended_features,
    plots                 = plots,
    params = list(
      lfc_threshold       = lfc_threshold,
      min_sample_fraction = min_sample_fraction,
      min_cell_types      = min_cell_types,
      min_cells           = min_cells,
      date_run            = Sys.time()
    )
  )
}


# ── Internal plot helpers ─────────────────────────────────────────────────────

#' @keywords internal
.make_ambient_plots <- function(per_ct_results, cross_ct_df, cons_mat,
                                cell_types, min_sample_fraction, save_dir) {

  plots <- list()

  # ── 1. Consistency heatmap (genes x cell types) ───────────────────────────
  if (nrow(cons_mat) > 0) {
    # Limit to top 40 genes for readability
    plot_mat <- head(cons_mat, 40)

    heatmap_df <- reshape2::melt(
      data.frame(gene = rownames(plot_mat), plot_mat, check.names = FALSE),
      id.vars      = "gene",
      variable.name = "cell_type",
      value.name   = "fraction_corrected"
    )
    heatmap_df$fraction_corrected <- as.numeric(heatmap_df$fraction_corrected)
    heatmap_df$gene <- factor(heatmap_df$gene,
                              levels = rev(rownames(plot_mat)))

    heatmap_plot <- ggplot2::ggplot(
      heatmap_df,
      ggplot2::aes(x = cell_type, y = gene, fill = fraction_corrected)
    ) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
      ggplot2::scale_fill_gradient2(
        low      = "#F1EFE8",
        mid      = "#97D0A0",
        high     = "#1D4E2A",
        midpoint = 0.5,
        limits   = c(0, 1),
        name     = "Fraction of\nsamples corrected"
      ) +
      ggplot2::geom_hline(
        yintercept = nrow(plot_mat) - cumsum(
          table(rowSums(plot_mat >= min_sample_fraction))
        ) + 0.5,
        colour = "grey40", linewidth = 0.3, linetype = "dashed"
      ) +
      ggplot2::labs(
        title    = "Ambient gene consistency across cell types and samples",
        subtitle = paste0("Fraction of samples in which each gene was corrected ",
                          "(dark = consistently ambient)"),
        x        = "Cell type",
        y        = "Gene"
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
        axis.text.y      = ggplot2::element_text(face  = "italic", size = 8),
        plot.title       = ggplot2::element_text(face  = "bold"),
        panel.grid       = ggplot2::element_blank(),
        legend.position  = "right"
      )
    plots$heatmap <- heatmap_plot
  }

  # ── 2. Cross-cell-type offender bar chart ─────────────────────────────────
  if (nrow(cross_ct_df) > 0) {
    top_cross <- head(cross_ct_df, 30)
    cross_bar <- ggplot2::ggplot(
      top_cross,
      ggplot2::aes(
        x    = reorder(gene, n_cell_types_affected + mean_fraction_corrected),
        y    = mean_fraction_corrected,
        fill = n_cell_types_affected
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_gradient(
        low  = "#97D0A0",
        high = "#1D4E2A",
        name = "# cell types\naffected"
      ) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title    = "Cross-cell-type ambient offenders",
        subtitle = "Mean fraction of samples corrected across affected cell types",
        x        = NULL,
        y        = "Mean fraction of samples corrected"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    plots$cross_cell_type_bar <- cross_bar
  }

  # ── 3. Per-sample consistency dot plot (one panel per cell type) ──────────
  per_sample_df <- do.call(rbind, lapply(names(per_ct_results), function(ct) {
    df <- per_ct_results[[ct]]
    if (nrow(df) == 0) return(NULL)
    top_genes <- head(df$gene[df$is_consistent_ambient], 15)
    if (length(top_genes) == 0) return(NULL)
    df_top <- df[df$gene %in% top_genes, ]
    df_top$cell_type <- ct
    df_top
  }))

  if (!is.null(per_sample_df) && nrow(per_sample_df) > 0) {
    dot_plot <- ggplot2::ggplot(
      per_sample_df,
      ggplot2::aes(
        x    = cell_type,
        y    = reorder(gene, fraction_corrected),
        size = fraction_corrected,
        colour = mean_lfc_corrected
      )
    ) +
      ggplot2::geom_point() +
      ggplot2::scale_size_continuous(
        range  = c(1, 6),
        name   = "Fraction\ncorrected"
      ) +
      ggplot2::scale_colour_gradient2(
        low      = "#1D4E2A",
        mid      = "#97D0A0",
        high     = "#F1EFE8",
        midpoint = -0.3,
        name     = "Mean\nlog2FC"
      ) +
      ggplot2::labs(
        title    = "Per-cell-type ambient consistency",
        subtitle = "Top consistently corrected genes | dot size = fraction of samples | colour = mean log2FC",
        x        = "Cell type",
        y        = "Gene"
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.text.y = ggplot2::element_text(face = "italic", size = 8),
        plot.title  = ggplot2::element_text(face = "bold")
      )
    plots$per_sample_consistency <- dot_plot
  }

  # Save plots
  if (!is.null(save_dir)) {
    if (!is.null(plots$heatmap)) {
      ggplot2::ggsave(
        file.path(save_dir, "ambient_consistency_heatmap.pdf"),
        plots$heatmap,
        width = max(8, length(names(per_ct_results)) * 1.2),
        height = min(20, max(6, nrow(cons_mat) * 0.35))
      )
    }
    if (!is.null(plots$cross_cell_type_bar)) {
      ggplot2::ggsave(
        file.path(save_dir, "cross_celltype_ambient_bar.pdf"),
        plots$cross_cell_type_bar, width = 8, height = 6
      )
    }
    if (!is.null(plots$per_sample_consistency)) {
      ggplot2::ggsave(
        file.path(save_dir, "per_sample_consistency_dot.pdf"),
        plots$per_sample_consistency, width = 10, height = 8
      )
    }
  }

  plots
}


# ── export_ambient_offenders ──────────────────────────────────────────────────

#' Export a shareable ambient offenders reference table
#'
#' Produces a clean, publication-ready table of ambient RNA offenders identified
#' by [identify_consistent_ambient()], suitable for sharing with the field as a
#' community reference resource. Saves as CSV, RDS, and optionally Excel.
#'
#' The ambient offenders table captures which genes are consistently detected
#' as ambient contaminants — across samples *and* across cell types — in a
#' given tissue/dataset context. For islet scRNA-seq, this is a practically
#' important resource because the INS/TTR ambient problem is universal and
#' reproducible community references for it do not yet exist.
#'
#' @param ambient_map The list returned by [identify_consistent_ambient()].
#' @param dataset_name Character. Short name for your dataset (e.g.
#'   `"HPAP_islet_atlas"`). Used in the output metadata header.
#' @param tissue Character. Tissue of origin (e.g. `"human_pancreatic_islet"`).
#' @param n_donors Integer or NULL. Number of donors in the dataset. If NULL,
#'   inferred from the results where possible. Default `NULL`.
#' @param platform Character. Sequencing platform (e.g. `"10x_Chromium_v3"`).
#'   Default `"10x_Chromium"`.
#' @param notes Character or NULL. Free-text notes about the dataset or
#'   processing. Default `NULL`.
#' @param save_dir Character. Directory to save outputs. Required.
#' @param save_excel Logical. Also save as `.xlsx` (requires `openxlsx`).
#'   Default `TRUE`.
#'
#' @return A data frame — the ambient offenders table — invisibly. Also writes:
#' \describe{
#'   \item{`ambient_offenders.csv`}{Full offenders table with metadata header.}
#'   \item{`ambient_offenders.rds`}{R object for direct loading in other
#'     analyses.}
#'   \item{`ambient_offenders.xlsx`}{Excel file with two sheets:
#'     `cross_cell_type` (genes ambient in multiple cell types) and
#'     `per_cell_type` (full per-cell-type gene lists). Only if
#'     `save_excel = TRUE`.}
#' }
#'
#' @examples
#' \dontrun{
#' ambient_map <- identify_consistent_ambient(results, seu, ...)
#'
#' offenders <- export_ambient_offenders(
#'   ambient_map  = ambient_map,
#'   dataset_name = "HPAP_islet_atlas_2026",
#'   tissue       = "human_pancreatic_islet",
#'   n_donors     = 87,
#'   platform     = "10x_Chromium_v3",
#'   notes        = "CD45-enriched + whole pancreas samples; SS2026 atlas",
#'   save_dir     = "ambient_filter_results"
#' )
#' }
#'
#' @export
export_ambient_offenders <- function(ambient_map,
                                     dataset_name = "my_dataset",
                                     tissue       = "unknown",
                                     n_donors     = NULL,
                                     platform     = "10x_Chromium",
                                     notes        = NULL,
                                     save_dir,
                                     save_excel   = TRUE) {

  if (missing(save_dir) || is.null(save_dir)) {
    rlang::abort("`save_dir` is required for export_ambient_offenders().")
  }
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  cross_ct  <- ambient_map$cross_cell_type
  per_ct    <- ambient_map$per_cell_type
  cons_mat  <- ambient_map$consistency_matrix
  params    <- ambient_map$params

  # ── Build the full offenders table ─────────────────────────────────────────
  # Start with cross-cell-type offenders as the primary list
  if (nrow(cross_ct) > 0) {
    offenders <- cross_ct
  } else {
    # Fall back to any gene consistent in at least one cell type
    all_consistent <- do.call(rbind, lapply(names(per_ct), function(ct) {
      df <- per_ct[[ct]]
      df_yes <- df[df$is_consistent_ambient, ]
      if (nrow(df_yes) == 0) return(NULL)
      data.frame(
        gene                    = df_yes$gene,
        n_cell_types_affected   = 1L,
        cell_types_affected     = ct,
        mean_fraction_corrected = df_yes$fraction_corrected,
        stringsAsFactors        = FALSE
      )
    }))
    offenders <- if (!is.null(all_consistent)) all_consistent else
      data.frame(gene = character(0), n_cell_types_affected = integer(0),
                 cell_types_affected = character(0),
                 mean_fraction_corrected = numeric(0),
                 stringsAsFactors = FALSE)
  }

  # Add per-cell-type fraction details from consistency matrix
  if (nrow(cons_mat) > 0 && nrow(offenders) > 0) {
    mat_genes   <- rownames(cons_mat)
    ct_cols     <- colnames(cons_mat)
    frac_detail <- as.data.frame(cons_mat[
      intersect(mat_genes, offenders$gene), , drop = FALSE
    ])
    colnames(frac_detail) <- paste0("frac_corrected__", ct_cols)
    frac_detail$gene      <- rownames(frac_detail)
    offenders <- merge(offenders, frac_detail, by = "gene", all.x = TRUE)
  }

  offenders <- offenders[order(offenders$n_cell_types_affected,
                               offenders$mean_fraction_corrected,
                               decreasing = TRUE), ]
  rownames(offenders) <- NULL

  # ── Metadata header ────────────────────────────────────────────────────────
  metadata <- data.frame(
    field = c(
      "dataset_name", "tissue", "n_donors", "platform",
      "date_generated", "ambientfilter_version",
      "lfc_threshold", "min_sample_fraction", "min_cell_types",
      "n_offenders_cross_celltype", "n_cell_types_analysed", "notes"
    ),
    value = c(
      dataset_name,
      tissue,
      as.character(n_donors %||% "unknown"),
      platform,
      as.character(Sys.Date()),
      as.character(utils::packageVersion("AmbientFilter")),
      as.character(params$lfc_threshold),
      as.character(params$min_sample_fraction),
      as.character(params$min_cell_types),
      as.character(nrow(cross_ct)),
      as.character(length(per_ct)),
      as.character(notes %||% "")
    ),
    stringsAsFactors = FALSE
  )

  # ── Save CSV (metadata block + data) ─────────────────────────────────────
  csv_path <- file.path(save_dir, "ambient_offenders.csv")

  # Write metadata as commented header lines, then the data
  meta_lines <- apply(metadata, 1, function(r) paste0("# ", r[1], ": ", r[2]))
  writeLines(c(
    "# AmbientFilter — Ambient RNA Offenders Reference Table",
    paste0("# Generated: ", Sys.time()),
    "#",
    "# This table reports genes identified as consistent ambient RNA",
    "# contaminants across donors and cell types. Genes listed here",
    "# appear as significantly corrected (reduced) by removeAmbience()",
    "# in >= min_sample_fraction of donors for >= min_cell_types cell types.",
    "# Share this file to help other labs clean similar datasets.",
    "#",
    meta_lines,
    "#"
  ), csv_path)
  suppressWarnings(
    utils::write.table(offenders, csv_path, append = TRUE, sep = ",",
                       row.names = FALSE, col.names = TRUE, quote = TRUE)
  )

  # ── Save RDS ──────────────────────────────────────────────────────────────
  rds_obj <- list(
    offenders    = offenders,
    per_cell_type = per_ct,
    metadata     = metadata,
    params       = params
  )
  saveRDS(rds_obj, file.path(save_dir, "ambient_offenders.rds"))

  # ── Save Excel ────────────────────────────────────────────────────────────
  if (save_excel) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      rlang::warn(
        "openxlsx not installed — skipping Excel export. ",
        "Install with: install.packages('openxlsx')"
      )
    } else {
      wb <- openxlsx::createWorkbook()

      # Sheet 1: metadata
      openxlsx::addWorksheet(wb, "metadata")
      openxlsx::writeData(wb, "metadata", metadata)
      openxlsx::setColWidths(wb, "metadata", cols = 1:2, widths = c(30, 50))

      # Sheet 2: cross-cell-type offenders (the headline list)
      openxlsx::addWorksheet(wb, "cross_celltype_offenders")
      openxlsx::writeData(wb, "cross_celltype_offenders", offenders)
      # Highlight header
      header_style <- openxlsx::createStyle(
        fgFill = "#4472C4", fontColour = "#FFFFFF",
        textDecoration = "bold", halign = "left"
      )
      openxlsx::addStyle(wb, "cross_celltype_offenders",
                         style = header_style, rows = 1,
                         cols = 1:ncol(offenders), gridExpand = TRUE)
      openxlsx::setColWidths(wb, "cross_celltype_offenders",
                             cols = 1:ncol(offenders), widths = "auto")

      # Sheet 3: consistency matrix
      if (nrow(cons_mat) > 0) {
        openxlsx::addWorksheet(wb, "consistency_matrix")
        mat_df <- as.data.frame(cons_mat)
        mat_df <- cbind(gene = rownames(mat_df), mat_df)
        openxlsx::writeData(wb, "consistency_matrix", mat_df)
        openxlsx::addStyle(wb, "consistency_matrix",
                           style = header_style, rows = 1,
                           cols = 1:ncol(mat_df), gridExpand = TRUE)
        # Conditional formatting — green scale on fraction values
        openxlsx::conditionalFormatting(
          wb, "consistency_matrix",
          cols = 2:ncol(mat_df),
          rows = 2:(nrow(mat_df) + 1),
          style = c("#FFFFFF", "#1D4E2A"),
          rule  = c(0, 1),
          type  = "colourScale"
        )
        openxlsx::setColWidths(wb, "consistency_matrix",
                               cols = 1:ncol(mat_df), widths = "auto")
      }

      # Sheet 4-N: per-cell-type tables
      for (ct in names(per_ct)) {
        df <- per_ct[[ct]]
        if (nrow(df) == 0) next
        # Sanitise sheet name (Excel max 31 chars, no special chars)
        sheet_name <- substr(gsub("[^A-Za-z0-9_]", "_", ct), 1, 31)
        openxlsx::addWorksheet(wb, sheet_name)
        openxlsx::writeData(wb, sheet_name, df)
        openxlsx::addStyle(wb, sheet_name,
                           style = header_style, rows = 1,
                           cols = 1:ncol(df), gridExpand = TRUE)
        # Highlight consistent ambient rows
        if (any(df$is_consistent_ambient)) {
          ambient_style <- openxlsx::createStyle(fgFill = "#D4EDD8")
          ambient_rows  <- which(df$is_consistent_ambient) + 1
          openxlsx::addStyle(wb, sheet_name, style = ambient_style,
                             rows = ambient_rows, cols = 1:ncol(df),
                             gridExpand = TRUE)
        }
        openxlsx::setColWidths(wb, sheet_name,
                               cols = 1:ncol(df), widths = "auto")
      }

      openxlsx::saveWorkbook(
        wb, file.path(save_dir, "ambient_offenders.xlsx"),
        overwrite = TRUE
      )
    }
  }

  cli::cli_alert_success(
    paste0("Ambient offenders exported to: ", save_dir, "/ambient_offenders.*")
  )
  cli::cli_alert_info(
    paste0(nrow(offenders), " ambient offender genes identified | ",
           "share ambient_offenders.csv + .rds with the field!")
  )

  invisible(offenders)
}
