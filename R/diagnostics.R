# ── AmbientFilter — diagnostics ───────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────────────────────

#' Plot correction diagnostics for a single cell type
#'
#' Produces three diagnostic visualisations for a corrected cell type:
#' a correlation plot of original vs corrected mean expression, violin plots
#' for the most-changed genes, and a bar chart of the top changed genes
#' ranked by log2FC magnitude.
#'
#' @param result A list as returned by [correct_cell_type()], containing
#'   `original`, `corrected`, and `changed_genes` elements.
#' @param cell_label Character. Label for plot titles.
#' @param assay Character. Assay to use. Default `"RNA"`.
#' @param label_difference Numeric. Minimum |log2FC| to label a gene in the
#'   correlation plot. Default `0.5`.
#' @param n_violin Integer. Number of top-changed genes to show as violin
#'   plots. Default `6`.
#' @param highlight_genes Character vector or NULL. Specific genes to always
#'   highlight in plots regardless of their log2FC (e.g. `c("INS", "TTR",
#'   "GCG")` for islet datasets). Default `NULL`.
#'
#' @return A named list of ggplot objects:
#' \describe{
#'   \item{`correlation`}{Scatter plot of original vs corrected mean expression.}
#'   \item{`violin_original`}{Violin plots of the top-changed genes in the
#'     original (uncorrected) object.}
#'   \item{`violin_corrected`}{Violin plots of the top-changed genes in the
#'     corrected object.}
#'   \item{`top_genes_bar`}{Bar chart of top changed genes by |log2FC|.}
#' }
#'
#' @export
plot_correction_diagnostics <- function(result,
                                        cell_label       = "cell type",
                                        assay            = "RNA",
                                        label_difference = 0.5,
                                        n_violin         = 6,
                                        highlight_genes  = NULL) {

  changed   <- result$changed_genes
  orig_obj  <- result$original
  corr_obj  <- result$corrected

  orig_means <- Matrix::rowMeans(.get_counts(orig_obj, assay))
  corr_means <- Matrix::rowMeans(.get_counts(corr_obj, assay))

  shared <- intersect(names(orig_means), names(corr_means))
  plot_df <- data.frame(
    gene           = shared,
    mean_original  = as.numeric(orig_means[shared]),
    mean_corrected = as.numeric(corr_means[shared]),
    stringsAsFactors = FALSE
  )
  plot_df$log2fc    <- log2((plot_df$mean_corrected + 1) /
                            (plot_df$mean_original  + 1))
  plot_df$abs_log2fc <- abs(plot_df$log2fc)

  # Genes to label
  label_genes <- unique(c(
    plot_df$gene[plot_df$abs_log2fc >= label_difference],
    highlight_genes
  ))
  plot_df$label <- ifelse(plot_df$gene %in% label_genes,
                          plot_df$gene, NA_character_)

  # ── Correlation plot ─────────────────────────────────────────────────────
  corr_plot <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = mean_original, y = mean_corrected)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = abs_log2fc),
      size = 0.6, alpha = 0.5
    ) +
    ggplot2::scale_colour_gradient(
      low  = "grey70",
      high = "#C94040",
      name = "|log2FC|"
    ) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    ggplot2::scale_x_continuous(trans = "log1p") +
    ggplot2::scale_y_continuous(trans = "log1p") +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      size = 2.5, hjust = -0.1, vjust = 0.5,
      na.rm = TRUE, colour = "#333333"
    ) +
    ggplot2::labs(
      title    = paste0(cell_label, " — ambient correction correlation"),
      subtitle = paste0("Points below the line = reduced by correction | ",
                        "labelled: |log2FC| >= ", label_difference),
      x        = "Mean counts (original, log1p)",
      y        = "Mean counts (corrected, log1p)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )

  # ── Top changed genes bar chart ──────────────────────────────────────────
  top_changed <- head(changed, 20)
  bar_plot <- ggplot2::ggplot(
    top_changed,
    ggplot2::aes(
      x    = reorder(gene, abs_log2fc),
      y    = log2fc,
      fill = ifelse(log2fc < 0, "Reduced", "Increased")
    )
  ) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::scale_fill_manual(
      values = c("Reduced" = "#5FAD6E", "Increased" = "#C94040")
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title    = paste0(cell_label, " — top changed genes"),
      subtitle = "log2FC of corrected vs original mean counts",
      x        = NULL,
      y        = "log2FC (corrected / original)"
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  # ── Violin plots ─────────────────────────────────────────────────────────
  top_violin_genes <- head(changed$gene, n_violin)

  .make_violin <- function(obj, genes, title_suffix) {
    if (length(genes) == 0 || ncol(obj) == 0) return(NULL)
    # Only plot genes present in object
    genes <- genes[genes %in% rownames(obj)]
    if (length(genes) == 0) return(NULL)

    Seurat::VlnPlot(
      obj,
      features = genes,
      assay    = assay,
      ncol     = min(3, length(genes)),
      pt.size  = 0
    ) +
      patchwork::plot_annotation(
        title = paste0(cell_label, " — ", title_suffix)
      )
  }

  vln_orig <- .make_violin(orig_obj, top_violin_genes, "original (uncorrected)")
  vln_corr <- .make_violin(corr_obj, top_violin_genes, "corrected")

  list(
    correlation      = corr_plot,
    violin_original  = vln_orig,
    violin_corrected = vln_corr,
    top_genes_bar    = bar_plot
  )
}

#' Summarise ambient correction across all cell types
#'
#' Generates a summary data frame of correction statistics across all
#' cell types processed by [ambient_filter()].
#'
#' @param results The list returned by [ambient_filter()].
#'
#' @return A data frame with one row per cell type, columns:
#'   `cell_type`, `n_cells`, `samples_used`, `samples_skipped`,
#'   `n_genes_changed`, `top_reduced_gene`, `top_reduced_log2fc`.
#'
#' @export
summarise_correction <- function(results) {
  cell_types <- names(results$cell_type_results)

  rows <- lapply(cell_types, function(ct) {
    r <- results$cell_type_results[[ct]]
    changed <- r$changed_genes
    reduced <- changed[changed$log2fc < 0, ]

    data.frame(
      cell_type         = ct,
      n_cells           = ncol(r$corrected),
      samples_used      = length(r$samples_used),
      samples_skipped   = length(r$samples_skipped),
      n_genes_changed   = nrow(changed),
      top_reduced_gene  = if (nrow(reduced) > 0) reduced$gene[1] else NA_character_,
      top_reduced_log2fc= if (nrow(reduced) > 0) round(reduced$log2fc[1], 3) else NA_real_,
      stringsAsFactors  = FALSE
    )
  })

  do.call(rbind, rows)
}
