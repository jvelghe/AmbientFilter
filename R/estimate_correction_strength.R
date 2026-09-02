# ── AmbientFilter — data-driven correction strength estimation ─────────────────
# ────────────────────────────────────────────────────────────────────────────────

#' Estimate per-cell-type correction strength from data
#'
#' Instead of manually specifying correction strengths (e.g. `macrophage = 0.5`,
#' `lymphocyte = 1.0`), this function learns appropriate correction strengths
#' from the data itself by exploiting a key property of ambient contamination:
#' **ambient genes track the background profile across samples, while genuinely
#' expressed genes do not.**
#'
#' For each cell type, the function computes — across all samples — the
#' Spearman correlation between each gene's mean expression in that cell type
#' and its proportion in the per-sample ambient background. A gene whose
#' expression in macrophages rises and falls in lockstep with the ambient
#' INS proportion across samples is almost certainly ambient. A gene whose
#' expression is independent of the ambient profile is almost certainly real.
#'
#' The **ambient correlation score** for a cell type is the mean of the top-gene
#' correlations — i.e. the degree to which the highest-expressed genes track
#' the background. This is used as the correction strength: a cell type where
#' most expression tracks the background gets a high strength (aggressive
#' correction); a cell type where expression is largely independent of the
#' background gets a low strength (conservative correction).
#'
#' Macrophages, which genuinely internalise ambient RNA via efferocytosis, will
#' automatically receive a lower estimated strength because their
#' background-independent expression (their own marker genes) will dilute the
#' ambient correlation signal. Lymphocytes, whose "islet" signal is purely
#' ambient, will receive a higher estimated strength.
#'
#' @section Method details:
#' For each cell type `c` and sample `s`:
#' 1. Compute mean expression vector `mu_cs` (genes x 1).
#' 2. Retrieve ambient profile `pi_s` (genes x 1, gene proportions).
#' 3. Compute Spearman correlation `rho_cs = cor(mu_cs, pi_s, method="spearman")`.
#'
#' Across all samples for cell type `c`:
#' 4. `ambient_correlation_c = median(rho_cs)` across samples (median is more
#'    robust than mean to outlier samples with unusual processing).
#' 5. `correction_strength_c = clamp(ambient_correlation_c, lower, upper)` —
#'    optionally clamped to `[min_strength, max_strength]` to prevent
#'    degenerate values.
#'
#' Optional refinement — **gene-level weighting**:
#' Rather than one correlation per sample, compute per-gene across-sample
#' correlations and use the distribution to build a gene-level correction
#' weight vector. This is returned in `gene_weights` and can be passed to
#' a future gene-level version of `correct_cell_type()`.
#'
#' @param object A Seurat object.
#' @param background A named list of background profiles as returned by
#'   [estimate_background()].
#' @param cell_type_col Character. Metadata column containing cell type labels.
#'   Default `"cell_type"`.
#' @param sample_col Character. Metadata column containing sample IDs.
#'   Default `"orig.ident"`.
#' @param tissue_col Character or NULL. Tissue metadata column. Must match
#'   what was used in [estimate_background()]. Default `NULL`.
#' @param assay Character. Assay to use. Default `"RNA"`.
#' @param cell_types Character vector or NULL. Cell types to estimate strength
#'   for. If NULL, all cell types in `cell_type_col`. Default `NULL`.
#' @param min_cells Integer. Minimum cells per cell-type/sample combination
#'   to include. Default `10`.
#' @param min_samples Integer. Minimum number of samples required for a
#'   reliable correlation estimate. Cell types with fewer samples receive a
#'   fallback strength. Default `3`.
#' @param min_strength Numeric in \[0,1\]. Floor for estimated correction
#'   strength. Prevents over-correction even when correlation is high. A value
#'   of `0.3` means even "pure ambient" cell types are corrected at most 100%
#'   but at least 30% — preserving some signal. Default `0.2`.
#' @param max_strength Numeric in \[0,1\]. Ceiling for estimated correction
#'   strength. Default `1.0`.
#' @param fallback_strength Numeric. Strength to assign when a cell type has
#'   fewer than `min_samples` samples. Default `0.5` (conservative).
#' @param top_genes_pct Numeric in (0,1\]. Fraction of highest-expressed genes
#'   to use for the correlation calculation. Using all genes can dilute the
#'   signal with many near-zero genes; using the top 20-50% focuses the
#'   estimate on expressed genes. Default `0.3`.
#' @param compute_gene_weights Logical. Whether to compute per-gene correction
#'   weights (more expensive but enables gene-level correction). Default `TRUE`.
#' @param prior_strengths Named numeric vector or NULL. Optional manual priors.
#'   When provided, the estimated strength is blended with the prior:
#'   `final = (1 - prior_weight) * estimated + prior_weight * prior`.
#'   Useful for macrophages where you want the estimate to be pulled toward
#'   your biological knowledge. Default `NULL`.
#' @param prior_weight Numeric in \[0,1\]. Weight given to `prior_strengths`
#'   in the blend. `0` = fully data-driven; `1` = fully manual (same as not
#'   using this function). Default `0.3`.
#' @param make_plots Logical. Generate diagnostic plots. Default `TRUE`.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A list with:
#' \describe{
#'   \item{`correction_strength`}{Named numeric vector of estimated correction
#'     strengths, one per cell type. Drop this directly into `ambient_filter()`
#'     as the `correction_strength` argument.}
#'   \item{`per_sample_correlations`}{Named list (one per cell type) of numeric
#'     vectors giving the per-sample Spearman correlations used to derive the
#'     estimate.}
#'   \item{`gene_weights`}{Named list (one per cell type) of named numeric
#'     vectors giving per-gene correction weights in \[0,1\], derived from
#'     cross-sample gene-background correlations. Higher weight = more
#'     confidently ambient. `NULL` if `compute_gene_weights = FALSE`.}
#'   \item{`summary`}{Data frame summarising the estimate per cell type:
#'     `cell_type`, `n_samples`, `median_correlation`, `estimated_strength`,
#'     `prior_strength` (if used), `final_strength`, `method`.}
#'   \item{`plots`}{List of ggplot objects (if `make_plots = TRUE`).}
#' }
#'
#' @examples
#' \dontrun{
#' bg <- estimate_background(seu, sample_col = "hpap_id", tissue_col = "tissue")
#'
#' # Fully data-driven
#' strength_est <- estimate_correction_strength(
#'   object     = seu,
#'   background = bg,
#'   sample_col = "hpap_id",
#'   tissue_col = "tissue"
#' )
#'
#' # Inspect what the model learned
#' strength_est$summary
#'
#' # With a biological prior for macrophages
#' strength_est <- estimate_correction_strength(
#'   object          = seu,
#'   background      = bg,
#'   prior_strengths = c(macrophage = 0.5),
#'   prior_weight    = 0.4   # 40% prior, 60% data
#' )
#'
#' # Use estimated strengths in the main pipeline
#' results <- ambient_filter(
#'   object              = seu,
#'   background          = bg,
#'   correction_strength = strength_est$correction_strength,
#'   ...
#' )
#' }
#'
#' @export
estimate_correction_strength <- function(object,
                                         background,
                                         cell_type_col       = "cell_type",
                                         sample_col          = "orig.ident",
                                         tissue_col          = NULL,
                                         assay               = "RNA",
                                         cell_types          = NULL,
                                         min_cells           = 10,
                                         min_samples         = 3,
                                         min_strength        = 0.2,
                                         max_strength        = 1.0,
                                         fallback_strength   = 0.5,
                                         top_genes_pct       = 0.3,
                                         compute_gene_weights = TRUE,
                                         prior_strengths     = NULL,
                                         prior_weight        = 0.3,
                                         make_plots          = TRUE,
                                         verbose             = TRUE) {

  .validate_seurat(object, assay)

  meta   <- object@meta.data
  counts <- .get_counts(object, assay)

  # ── Build sample key ───────────────────────────────────────────────────────
  if (!is.null(tissue_col) && tissue_col %in% colnames(meta)) {
    sample_key <- paste(meta[[tissue_col]], meta[[sample_col]], sep = ".")
  } else {
    sample_key <- meta[[sample_col]]
  }

  # ── Resolve cell types ─────────────────────────────────────────────────────
  all_types <- unique(as.character(meta[[cell_type_col]]))
  if (is.null(cell_types)) cell_types <- all_types
  cell_types <- intersect(cell_types, all_types)

  if (verbose) {
    cli::cli_h1("Estimating data-driven correction strengths")
    cli::cli_alert_info(
      paste0("Method: cross-sample ambient correlation | ",
             length(cell_types), " cell types | ",
             "top_genes_pct = ", top_genes_pct)
    )
  }

  n_genes <- nrow(counts)
  n_top   <- max(10L, floor(n_genes * top_genes_pct))

  # ── Per-cell-type correlation analysis ────────────────────────────────────
  per_ct_corr    <- vector("list", length(cell_types))
  gene_weights   <- vector("list", length(cell_types))
  names(per_ct_corr)  <- cell_types
  names(gene_weights) <- cell_types

  pb <- if (verbose) {
    cli::cli_progress_bar("Computing ambient correlations", total = length(cell_types))
  } else NULL

  for (ct in cell_types) {

    ct_cells <- rownames(meta)[meta[[cell_type_col]] == ct]
    if (length(ct_cells) == 0) {
      per_ct_corr[[ct]] <- numeric(0)
      if (verbose) cli::cli_progress_update(id = pb)
      next
    }

    ct_meta      <- meta[ct_cells, , drop = FALSE]
    ct_counts    <- counts[, ct_cells, drop = FALSE]
    ct_sample_key <- if (!is.null(tissue_col) && tissue_col %in% colnames(ct_meta)) {
      paste(ct_meta[[tissue_col]], ct_meta[[sample_col]], sep = ".")
    } else {
      ct_meta[[sample_col]]
    }

    samples <- intersect(unique(ct_sample_key), names(background))

    # Per-sample correlation: mean expression in this cell type vs ambient profile
    sample_corrs       <- numeric(length(samples))
    names(sample_corrs) <- samples

    # For gene-level weights: collect per-gene x per-sample expression
    gene_sample_mat <- matrix(
      NA_real_,
      nrow = n_genes,
      ncol = length(samples),
      dimnames = list(rownames(counts), samples)
    )

    for (i_s in seq_along(samples)) {
      s       <- samples[i_s]
      s_cells <- which(ct_sample_key == s)

      if (length(s_cells) < min_cells) {
        sample_corrs[i_s] <- NA_real_
        next
      }

      s_means <- Matrix::rowMeans(ct_counts[, s_cells, drop = FALSE])

      # Align with background
      bg_vec  <- background[[s]]
      shared  <- intersect(names(s_means), names(bg_vec))
      if (length(shared) < 20) {
        sample_corrs[i_s] <- NA_real_
        next
      }

      s_means_shared <- as.numeric(s_means[shared])
      bg_shared      <- as.numeric(bg_vec[shared])

      # Focus on top expressed genes to avoid zero-inflation diluting correlation
      top_idx <- order(s_means_shared, decreasing = TRUE)[seq_len(
        min(n_top, length(shared))
      )]
      rho <- tryCatch(
        stats::cor(s_means_shared[top_idx], bg_shared[top_idx],
                   method = "spearman", use = "complete.obs"),
        error = function(e) NA_real_
      )
      sample_corrs[i_s] <- rho

      # Store full gene means for gene-weight computation
      gene_sample_mat[shared, s] <- s_means_shared
    }

    per_ct_corr[[ct]] <- sample_corrs[!is.na(sample_corrs)]

    # ── Gene-level weights ─────────────────────────────────────────────────
    if (compute_gene_weights) {
      valid_samples <- samples[!is.na(sample_corrs)]
      if (length(valid_samples) >= 2) {

        # For each gene: across-sample correlation of mean expression vs
        # ambient proportion — high correlation = confidently ambient
        gene_corrs <- vapply(rownames(counts), function(g) {
          gene_expr <- gene_sample_mat[g, valid_samples]
          bg_props  <- vapply(valid_samples, function(s) {
            if (g %in% names(background[[s]])) background[[s]][[g]] else 0
          }, numeric(1))
          if (sum(!is.na(gene_expr)) < 2 || stats::sd(gene_expr, na.rm=TRUE) == 0) {
            return(0)
          }
          tryCatch(
            stats::cor(gene_expr, bg_props,
                       method = "spearman", use = "complete.obs"),
            error = function(e) 0
          )
        }, numeric(1))

        # Rescale to [0, 1]: negative correlations → 0 (confidently NOT ambient)
        gene_corrs[is.na(gene_corrs)] <- 0
        gene_corrs_pos <- pmax(gene_corrs, 0)
        # Normalise to [0, 1]
        if (max(gene_corrs_pos) > 0) {
          gene_weights[[ct]] <- gene_corrs_pos / max(gene_corrs_pos)
        } else {
          gene_weights[[ct]] <- setNames(rep(0, n_genes), rownames(counts))
        }
      } else {
        gene_weights[[ct]] <- NULL
      }
    }

    if (verbose) cli::cli_progress_update(id = pb)
  }

  if (verbose) cli::cli_progress_done(id = pb)

  # ── Derive correction strengths from correlations ─────────────────────────
  summary_rows <- lapply(cell_types, function(ct) {
    corrs    <- per_ct_corr[[ct]]
    n_samp   <- length(corrs)
    method   <- "correlation"

    if (n_samp < min_samples) {
      median_cor  <- NA_real_
      est_strength <- fallback_strength
      method       <- "fallback"
    } else {
      median_cor   <- stats::median(corrs, na.rm = TRUE)
      # Map correlation to [min_strength, max_strength]
      # cor = 1 → max_strength (fully ambient)
      # cor = 0 → min_strength (fully independent)
      # cor < 0 → min_strength (expression anti-correlated with ambient)
      est_strength <- min_strength + pmax(median_cor, 0) *
                      (max_strength - min_strength)
      est_strength <- min(max(est_strength, min_strength), max_strength)
    }

    # Apply prior blend if provided
    prior_str <- if (!is.null(prior_strengths) && ct %in% names(prior_strengths)) {
      prior_strengths[[ct]]
    } else NA_real_

    if (!is.na(prior_str) && !is.null(prior_strengths)) {
      final_str <- (1 - prior_weight) * est_strength + prior_weight * prior_str
      method    <- paste0(method, "+prior(w=", prior_weight, ")")
    } else {
      final_str <- est_strength
    }

    data.frame(
      cell_type          = ct,
      n_samples          = n_samp,
      median_correlation = round(median_cor, 3),
      estimated_strength = round(est_strength, 3),
      prior_strength     = round(prior_str, 3),
      final_strength     = round(final_str, 3),
      method             = method,
      stringsAsFactors   = FALSE
    )
  })

  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  # Extract final correction strength vector
  correction_strength <- setNames(summary_df$final_strength,
                                  summary_df$cell_type)

  if (verbose) {
    cli::cli_h2("Estimated correction strengths")
    print(summary_df[, c("cell_type", "n_samples", "median_correlation",
                          "final_strength", "method")])
    cli::cli_alert_info(
      paste0("Use results$correction_strength directly in ambient_filter(). ",
             "Cell types with low correlation are conservatively corrected (biological signal preserved).")
    )
  }

  # ── Plots ─────────────────────────────────────────────────────────────────
  plots <- NULL
  if (make_plots) {
    plots <- .make_strength_plots(summary_df, per_ct_corr)
  }

  list(
    correction_strength   = correction_strength,
    per_sample_correlations = per_ct_corr,
    gene_weights          = if (compute_gene_weights) gene_weights else NULL,
    summary               = summary_df,
    plots                 = plots,
    params = list(
      top_genes_pct    = top_genes_pct,
      min_strength     = min_strength,
      max_strength     = max_strength,
      fallback_strength = fallback_strength,
      prior_strengths  = prior_strengths,
      prior_weight     = prior_weight,
      method           = "spearman_correlation",
      date_run         = Sys.time()
    )
  )
}


# ── Internal plot helpers ─────────────────────────────────────────────────────

#' @keywords internal
.make_strength_plots <- function(summary_df, per_ct_corr) {

  plots <- list()

  # ── 1. Estimated strength bar chart ───────────────────────────────────────
  bar_df <- summary_df[order(summary_df$final_strength, decreasing = TRUE), ]

  strength_bar <- ggplot2::ggplot(
    bar_df,
    ggplot2::aes(
      x    = reorder(cell_type, final_strength),
      y    = final_strength,
      fill = final_strength
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
                        colour = "grey40", linewidth = 0.5) +
    ggplot2::scale_fill_gradient2(
      low      = "#2D7D42",
      mid      = "#97D0A0",
      high     = "#C94040",
      midpoint = 0.5,
      limits   = c(0, 1),
      name     = "Correction\nstrength"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", final_strength)),
      hjust = -0.1, size = 3.2, colour = "grey20"
    ) +
    ggplot2::coord_flip(ylim = c(0, 1.1)) +
    ggplot2::labs(
      title    = "Data-driven correction strength estimates",
      subtitle = paste0("Derived from cross-sample Spearman correlation of\n",
                        "cell-type expression vs ambient background profile"),
      x        = NULL,
      y        = "Estimated correction strength (0 = none, 1 = full)"
    ) +
    ggplot2::annotate(
      "text", x = 0.6, y = 0.52,
      label = "conservative\n(biologically expressed)",
      hjust = 0, size = 2.8, colour = "grey50", fontface = "italic"
    ) +
    ggplot2::annotate(
      "text", x = 0.6, y = 0.98,
      label = "aggressive\n(likely ambient)",
      hjust = 1, size = 2.8, colour = "grey50", fontface = "italic"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )
  plots$strength_bar <- strength_bar

  # ── 2. Per-sample correlation violin plot ─────────────────────────────────
  corr_df <- do.call(rbind, lapply(names(per_ct_corr), function(ct) {
    corrs <- per_ct_corr[[ct]]
    if (length(corrs) == 0) return(NULL)
    data.frame(
      cell_type   = ct,
      correlation = corrs,
      stringsAsFactors = FALSE
    )
  }))

  if (!is.null(corr_df) && nrow(corr_df) > 0) {
    # Order cell types by median correlation
    ct_order <- names(sort(
      tapply(corr_df$correlation, corr_df$cell_type, stats::median),
      decreasing = TRUE
    ))
    corr_df$cell_type <- factor(corr_df$cell_type, levels = ct_order)

    corr_violin <- ggplot2::ggplot(
      corr_df,
      ggplot2::aes(x = cell_type, y = correlation, fill = cell_type)
    ) +
      ggplot2::geom_violin(scale = "width", trim = TRUE,
                           linewidth = 0.3, show.legend = FALSE) +
      ggplot2::geom_boxplot(width = 0.1, outlier.size = 0.5,
                            fill = "white", linewidth = 0.3) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          colour = "grey40", linewidth = 0.4) +
      ggplot2::scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.9) +
      ggplot2::labs(
        title    = "Per-sample ambient correlation by cell type",
        subtitle = paste0("Spearman correlation of cell-type mean expression vs ",
                          "ambient background profile per sample\n",
                          "Higher = expression tracks background = more likely ambient"),
        x        = "Cell type",
        y        = "Spearman correlation (expression vs ambient)"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title    = ggplot2::element_text(face = "bold")
      )
    plots$correlation_violin <- corr_violin
  }

  # ── 3. Scatter: median correlation vs final strength ─────────────────────
  valid_df <- summary_df[!is.na(summary_df$median_correlation), ]
  if (nrow(valid_df) > 0) {
    scatter <- ggplot2::ggplot(
      valid_df,
      ggplot2::aes(
        x     = median_correlation,
        y     = final_strength,
        label = cell_type,
        colour = grepl("prior", method)
      )
    ) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           linetype = "dotted", colour = "grey70") +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_text(hjust = -0.15, size = 3, show.legend = FALSE) +
      ggplot2::scale_colour_manual(
        values = c("FALSE" = "#2D7D42", "TRUE" = "#4472C4"),
        labels = c("FALSE" = "data-driven", "TRUE" = "prior-blended"),
        name   = "Method"
      ) +
      ggplot2::xlim(c(-0.1, 1.1)) +
      ggplot2::ylim(c(0, 1.1)) +
      ggplot2::labs(
        title    = "Median ambient correlation vs final correction strength",
        subtitle = "Diagonal = perfect mapping | deviation = prior blending or clamping",
        x        = "Median Spearman correlation (expression vs ambient)",
        y        = "Final correction strength applied"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    plots$correlation_vs_strength <- scatter
  }

  plots
}
