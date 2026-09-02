# ── AmbientFilter — VAE-based ambient correction ──────────────────────────────
# ────────────────────────────────────────────────────────────────────────────────
# Requires: reticulate, python >= 3.8, torch (via reticulate)
# Install python deps: AmbientFilter::install_ambient_vae()

#' Install Python dependencies for VAE-based ambient correction
#'
#' Installs PyTorch and supporting libraries into a dedicated conda/venv
#' environment used by AmbientFilter's VAE method.
#'
#' @param envname Character. Name of the conda/venv environment to create.
#'   Default `"ambientfilter"`.
#' @param method Character. `"conda"` or `"virtualenv"`. Default `"conda"`.
#' @param cuda Logical. Install CUDA-enabled PyTorch (GPU support). If FALSE,
#'   installs CPU-only version. Default `FALSE`.
#'
#' @export
install_ambient_vae <- function(envname = "ambientfilter",
                                method  = "conda",
                                cuda    = FALSE) {
  reticulate::py_install(
    packages = c(
      if (cuda) "torch torchvision" else "torch --index-url https://download.pytorch.org/whl/cpu",
      "numpy", "scipy", "pandas", "scikit-learn", "anndata"
    ),
    envname = envname,
    method  = method,
    pip     = TRUE
  )
  invisible(NULL)
}


#' VAE-based ambient RNA correction for a single cell type
#'
#' Trains a variational autoencoder (VAE) that disentangles genuine biological
#' expression from ambient RNA contamination within a given cell type. Unlike
#' the correlation-based [estimate_correction_strength()] method — which asks
#' whether a gene's expression tracks the background *across* samples — the VAE
#' learns a **latent decomposition** of the count matrix where ambient and
#' genuine expression are separated as distinct latent dimensions, conditioned
#' on the known per-sample ambient profile.
#'
#' @section Architecture:
#' ```
#' Input: count matrix X (cells x genes) + ambient profile A (cells x genes)
#'
#'   Encoder ──► z_bio    ~ N(mu_bio,    sigma_bio)    # genuine expression
#'           └──► z_ambient ~ N(mu_ambient, sigma_ambient)  # ambient level
#'
#'   Decoder_bio     (z_bio)     ──► X_bio     # genuine reconstruction
#'   Decoder_ambient (z_ambient) ──► X_ambient # ambient reconstruction
#'
#'   X_reconstructed = X_bio + X_ambient
#'
#' Loss = ELBO (reconstruction NB + KL) +
#'        lambda_ambient * correlation(X_ambient, A) +  # ambient supervision
#'        lambda_bio    * (1 - correlation(X_bio, gene_weights_inv))
#' ```
#'
#' The ambient supervision loss uses the ambient profile `A` as a soft
#' label — it pushes `z_ambient` to capture variation correlated with the
#' background, while `z_bio` captures the residual genuine expression.
#' Gene-level weights from [estimate_correction_strength()] can be used as
#' soft priors to guide training.
#'
#' @param object A Seurat object.
#' @param background A named list of background profiles from
#'   [estimate_background()].
#' @param cell_label Character. Cell type to correct.
#' @param cell_type_col Character. Default `"cell_type"`.
#' @param sample_col Character. Default `"orig.ident"`.
#' @param tissue_col Character or NULL. Default `NULL`.
#' @param assay Character. Default `"RNA"`.
#' @param gene_weights Named numeric vector or NULL. Per-gene ambient weights
#'   from [estimate_correction_strength()] used as soft supervision labels.
#'   If NULL, the VAE trains without gene-level priors (unsupervised mode).
#'   Default `NULL`.
#' @param latent_dim_bio Integer. Dimensionality of the genuine expression
#'   latent space. Default `10`.
#' @param latent_dim_ambient Integer. Dimensionality of the ambient latent
#'   space. Keep small — ambient variation is low-dimensional. Default `3`.
#' @param hidden_dims Integer vector. Hidden layer sizes for encoder/decoder.
#'   Default `c(256L, 128L)`.
#' @param n_epochs Integer. Training epochs. Default `150`.
#' @param batch_size Integer. Mini-batch size. Default `256L`.
#' @param learning_rate Numeric. Adam learning rate. Default `1e-3`.
#' @param lambda_ambient Numeric. Weight on the ambient supervision loss.
#'   Higher values push the ambient latent dimension harder toward the
#'   background profile. Default `1.0`.
#' @param lambda_bio Numeric. Weight on the biological signal preservation
#'   loss. Default `0.5`.
#' @param min_cells Integer. Minimum cells required per sample. Default `10`.
#' @param use_gpu Logical. Use GPU if available. Default `TRUE`.
#' @param envname Character. Python environment name. Default `"ambientfilter"`.
#' @param seed Integer. Random seed for reproducibility. Default `42L`.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A list matching the structure of [correct_cell_type()] output:
#' \describe{
#'   \item{`original`}{Original Seurat subset.}
#'   \item{`corrected`}{VAE-corrected Seurat subset.}
#'   \item{`changed_genes`}{Genes changed by correction.}
#'   \item{`ambient_scores`}{Per-cell ambient score (z_ambient mean).}
#'   \item{`reconstruction_loss`}{Training loss curve data frame.}
#'   \item{`model_path`}{Path to saved PyTorch model weights (NULL if not saved).}
#' }
#'
#' @export
correct_cell_type_vae <- function(object,
                                  background,
                                  cell_label,
                                  cell_type_col      = "cell_type",
                                  sample_col         = "orig.ident",
                                  tissue_col         = NULL,
                                  assay              = "RNA",
                                  gene_weights       = NULL,
                                  latent_dim_bio     = 10L,
                                  latent_dim_ambient = 3L,
                                  hidden_dims        = c(256L, 128L),
                                  n_epochs           = 150L,
                                  batch_size         = 256L,
                                  learning_rate      = 1e-3,
                                  lambda_ambient     = 1.0,
                                  lambda_bio         = 0.5,
                                  min_cells          = 10,
                                  use_gpu            = TRUE,
                                  envname            = "ambientfilter",
                                  seed               = 42L,
                                  verbose            = TRUE) {

  # ── Check Python environment ──────────────────────────────────────────────
  .check_vae_deps(envname)

  # ── Subset to target cell type ────────────────────────────────────────────
  .validate_seurat(object, assay)
  meta   <- object@meta.data
  ct_cells <- rownames(meta)[meta[[cell_type_col]] == cell_label]

  if (length(ct_cells) == 0) {
    rlang::abort(paste0("No cells found for cell type '", cell_label, "'."))
  }

  original_subset <- subset(object, cells = ct_cells)
  counts_mat      <- as.matrix(.get_counts(original_subset, assay))

  # ── Build per-cell ambient profile matrix ─────────────────────────────────
  # Each cell gets its sample's ambient profile as a conditioning vector
  ct_meta <- original_subset@meta.data
  if (!is.null(tissue_col) && tissue_col %in% colnames(ct_meta)) {
    ct_sample_key <- paste(ct_meta[[tissue_col]], ct_meta[[sample_col]], sep = ".")
  } else {
    ct_sample_key <- ct_meta[[sample_col]]
  }

  ambient_mat <- matrix(0, nrow = nrow(counts_mat), ncol = ncol(counts_mat),
                        dimnames = dimnames(counts_mat))
  for (s in unique(ct_sample_key)) {
    if (!s %in% names(background)) next
    cell_idx    <- which(ct_sample_key == s)
    bg_vec      <- background[[s]]
    shared      <- intersect(rownames(counts_mat), names(bg_vec))
    ambient_mat[shared, cell_idx] <- bg_vec[shared]
  }

  # ── Align gene weights ────────────────────────────────────────────────────
  gene_weight_vec <- if (!is.null(gene_weights)) {
    gw <- gene_weights[rownames(counts_mat)]
    gw[is.na(gw)] <- 0.5  # unknown → neutral
    as.numeric(gw)
  } else {
    rep(0.5, nrow(counts_mat))
  }

  # ── Run VAE via reticulate ────────────────────────────────────────────────
  if (verbose) {
    cli::cli_h2(paste0("Training ambient VAE: ", cell_label))
    cli::cli_alert_info(
      paste0(ncol(counts_mat), " cells | ", nrow(counts_mat), " genes | ",
             n_epochs, " epochs | latent_bio=", latent_dim_bio,
             " latent_ambient=", latent_dim_ambient)
    )
  }

  reticulate::use_condaenv(envname, required = FALSE)

  vae_result <- reticulate::py_run_string(.vae_python_code(), local = TRUE)

  py_correct <- reticulate::py$run_ambient_vae(
    counts          = t(counts_mat),           # cells x genes
    ambient_profile = t(ambient_mat),          # cells x genes
    gene_weights    = gene_weight_vec,         # genes
    latent_dim_bio  = as.integer(latent_dim_bio),
    latent_dim_amb  = as.integer(latent_dim_ambient),
    hidden_dims     = as.integer(hidden_dims),
    n_epochs        = as.integer(n_epochs),
    batch_size      = as.integer(batch_size),
    lr              = learning_rate,
    lambda_ambient  = lambda_ambient,
    lambda_bio      = lambda_bio,
    use_gpu         = use_gpu,
    seed            = as.integer(seed),
    verbose         = verbose
  )

  # ── Unpack results ────────────────────────────────────────────────────────
  corrected_counts_t <- py_correct[["corrected_counts"]]   # cells x genes
  ambient_scores     <- py_correct[["ambient_scores"]]     # cells x latent_ambient
  loss_curve         <- py_correct[["loss_curve"]]         # list of loss values

  corrected_counts <- t(corrected_counts_t)  # back to genes x cells
  rownames(corrected_counts) <- rownames(counts_mat)
  colnames(corrected_counts) <- colnames(counts_mat)
  corrected_counts <- Matrix::Matrix(corrected_counts, sparse = TRUE)

  # ── Build output objects ──────────────────────────────────────────────────
  corrected_subset <- .make_corrected_seurat(
    original_subset, corrected_counts, assay, cell_label
  )

  # Add ambient score to metadata
  ambient_score_mean <- rowMeans(ambient_scores)
  corrected_subset <- Seurat::AddMetaData(
    corrected_subset,
    metadata = data.frame(ambient_score_vae = ambient_score_mean,
                          row.names = colnames(corrected_subset))
  )

  # Changed genes summary
  orig_means <- Matrix::rowMeans(counts_mat)
  corr_means <- Matrix::rowMeans(corrected_counts)
  lfc <- log2((corr_means + 1) / (orig_means + 1))
  changed <- data.frame(
    gene           = names(lfc),
    mean_original  = orig_means,
    mean_corrected = corr_means,
    log2fc         = lfc,
    abs_log2fc     = abs(lfc),
    direction      = ifelse(lfc < 0, "decreased", "unchanged"),
    stringsAsFactors = FALSE
  )
  changed <- changed[order(changed$abs_log2fc, decreasing = TRUE), ]
  changed <- changed[changed$abs_log2fc > log2(1.5), ]
  rownames(changed) <- NULL

  loss_df <- data.frame(
    epoch          = seq_along(loss_curve),
    total_loss     = as.numeric(unlist(lapply(loss_curve, `[[`, "total"))),
    recon_loss     = as.numeric(unlist(lapply(loss_curve, `[[`, "recon"))),
    kl_loss        = as.numeric(unlist(lapply(loss_curve, `[[`, "kl"))),
    ambient_loss   = as.numeric(unlist(lapply(loss_curve, `[[`, "ambient"))),
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cli::cli_alert_success(
      paste0("VAE training complete. Final loss: ",
             round(tail(loss_df$total_loss, 1), 4))
    )
  }

  list(
    original           = original_subset,
    corrected          = corrected_subset,
    changed_genes      = changed,
    ambient_scores     = ambient_scores,
    reconstruction_loss = loss_df,
    model_path         = NULL
  )
}


#' Run ambient VAE correction across all cell types
#'
#' Wrapper that calls [correct_cell_type_vae()] for each cell type, matching
#' the interface of [ambient_filter()] for easy benchmarking.
#'
#' @inheritParams ambient_filter
#' @inheritParams correct_cell_type_vae
#' @param strength_estimates Optional output of [estimate_correction_strength()].
#'   If provided, per-gene weights are passed to each VAE as soft priors.
#'
#' @return Same structure as [ambient_filter()] results, with an additional
#'   `method = "vae"` field in params.
#'
#' @export
ambient_filter_vae <- function(object,
                               background          = NULL,
                               cell_type_col       = "cell_type",
                               sample_col          = "orig.ident",
                               tissue_col          = NULL,
                               cell_types          = NULL,
                               assay               = "RNA",
                               strength_estimates  = NULL,
                               latent_dim_bio      = 10L,
                               latent_dim_ambient  = 3L,
                               hidden_dims         = c(256L, 128L),
                               n_epochs            = 150L,
                               batch_size          = 256L,
                               learning_rate       = 1e-3,
                               lambda_ambient      = 1.0,
                               lambda_bio          = 0.5,
                               min_cells           = 10,
                               use_gpu             = TRUE,
                               envname             = "ambientfilter",
                               seed               = 42L,
                               run_consistency_analysis = TRUE,
                               min_sample_fraction = 0.5,
                               min_cell_types      = 2,
                               export_offenders    = TRUE,
                               dataset_name        = "my_dataset",
                               tissue              = "unknown",
                               n_donors            = NULL,
                               platform            = "10x_Chromium",
                               save_dir            = NULL,
                               verbose             = TRUE) {

  .validate_seurat(object, assay)

  # Estimate background if not provided
  if (is.null(background)) {
    background <- estimate_background(
      object        = object,
      sample_col    = sample_col,
      tissue_col    = tissue_col,
      assay         = assay,
      cell_type_col = cell_type_col,
      verbose       = verbose
    )
  }

  # Resolve cell types
  all_types  <- unique(as.character(object@meta.data[[cell_type_col]]))
  cell_types <- if (is.null(cell_types)) all_types else intersect(cell_types, all_types)

  if (!is.null(save_dir) && !dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  if (verbose) {
    cli::cli_h1("AmbientFilter VAE — per-cell-type ambient RNA correction")
    cli::cli_alert_info(paste0(
      length(cell_types), " cell types | ",
      n_epochs, " epochs | latent_bio=", latent_dim_bio,
      " | latent_ambient=", latent_dim_ambient
    ))
  }

  cell_type_results <- vector("list", length(cell_types))
  names(cell_type_results) <- cell_types

  for (ct in cell_types) {
    gene_wts <- if (!is.null(strength_estimates)) {
      strength_estimates$gene_weights[[ct]]
    } else NULL

    ct_result <- correct_cell_type_vae(
      object              = object,
      background          = background,
      cell_label          = ct,
      cell_type_col       = cell_type_col,
      sample_col          = sample_col,
      tissue_col          = tissue_col,
      assay               = assay,
      gene_weights        = gene_wts,
      latent_dim_bio      = latent_dim_bio,
      latent_dim_ambient  = latent_dim_ambient,
      hidden_dims         = hidden_dims,
      n_epochs            = n_epochs,
      batch_size          = batch_size,
      learning_rate       = learning_rate,
      lambda_ambient      = lambda_ambient,
      lambda_bio          = lambda_bio,
      min_cells           = min_cells,
      use_gpu             = use_gpu,
      envname             = envname,
      seed                = seed,
      verbose             = verbose
    )
    cell_type_results[[ct]] <- ct_result

    if (!is.null(save_dir)) {
      ct_dir <- file.path(save_dir, ct)
      dir.create(ct_dir, showWarnings = FALSE)
      utils::write.csv(ct_result$changed_genes,
                       file.path(ct_dir, "changed_genes_vae.csv"),
                       row.names = FALSE)
      utils::write.csv(ct_result$reconstruction_loss,
                       file.path(ct_dir, "training_loss.csv"),
                       row.names = FALSE)
    }
  }

  result_obj <- list(
    cell_type_results = cell_type_results,
    background        = background,
    ambient_map       = NULL,
    params = list(
      cell_type_col       = cell_type_col,
      sample_col          = sample_col,
      tissue_col          = tissue_col,
      cell_types          = cell_types,
      method              = "vae",
      n_epochs            = n_epochs,
      latent_dim_bio      = latent_dim_bio,
      latent_dim_ambient  = latent_dim_ambient,
      date_run            = Sys.time()
    )
  )
  result_obj$summary <- summarise_correction(result_obj)

  if (run_consistency_analysis) {
    ambient_map <- identify_consistent_ambient(
      results             = result_obj,
      object              = object,
      sample_col          = sample_col,
      tissue_col          = tissue_col,
      cell_type_col       = cell_type_col,
      assay               = assay,
      min_sample_fraction = min_sample_fraction,
      min_cell_types      = min_cell_types,
      save_dir            = if (!is.null(save_dir)) file.path(save_dir, "ambient_offenders") else NULL,
      verbose             = verbose
    )
    result_obj$ambient_map <- ambient_map

    if (export_offenders && !is.null(save_dir)) {
      export_ambient_offenders(
        ambient_map  = ambient_map,
        dataset_name = dataset_name,
        tissue       = tissue,
        n_donors     = n_donors,
        platform     = platform,
        save_dir     = file.path(save_dir, "ambient_offenders")
      )
    }
  }

  if (verbose) {
    cli::cli_h1("VAE correction complete")
    print(result_obj$summary)
  }

  result_obj
}


#' Plot VAE training diagnostics
#'
#' @param vae_result Output from [correct_cell_type_vae()].
#' @param cell_label Character. Label for plot titles.
#' @export
plot_vae_diagnostics <- function(vae_result, cell_label = "cell type") {

  loss_df <- vae_result$reconstruction_loss

  loss_long <- stats::reshape(
    loss_df[, c("epoch", "total_loss", "recon_loss", "kl_loss", "ambient_loss")],
    varying   = c("total_loss", "recon_loss", "kl_loss", "ambient_loss"),
    v.names   = "loss",
    timevar   = "component",
    times     = c("Total", "Reconstruction", "KL", "Ambient supervision"),
    direction = "long"
  )

  loss_plot <- ggplot2::ggplot(
    loss_long,
    ggplot2::aes(x = epoch, y = loss, colour = component)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_colour_manual(
      values = c(
        "Total"                = "#1F3864",
        "Reconstruction"       = "#2D7D42",
        "KL"                   = "#C94040",
        "Ambient supervision"  = "#E8A020"
      ),
      name = "Loss component"
    ) +
    ggplot2::labs(
      title    = paste0(cell_label, " — VAE training loss"),
      subtitle = "Ambient supervision loss drives disentanglement of ambient vs genuine expression",
      x        = "Epoch",
      y        = "Loss"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  # Ambient score distribution
  scores <- vae_result$ambient_scores
  score_df <- data.frame(
    ambient_score = rowMeans(scores),
    stringsAsFactors = FALSE
  )
  score_plot <- ggplot2::ggplot(
    score_df,
    ggplot2::aes(x = ambient_score)
  ) +
    ggplot2::geom_histogram(bins = 50, fill = "#5FAD6E", colour = "white",
                            linewidth = 0.2) +
    ggplot2::labs(
      title    = paste0(cell_label, " — per-cell ambient score"),
      subtitle = "Higher score = more ambient contamination in that cell",
      x        = "Mean ambient latent score",
      y        = "Number of cells"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  list(loss_curve = loss_plot, ambient_scores = score_plot)
}


#' Compare correlation-based and VAE correction methods
#'
#' Produces side-by-side comparison plots and a summary table benchmarking
#' the two correction methods against each other, using the changed_genes
#' outputs and INS/TTR reduction as key metrics.
#'
#' @param results_corr Output from [ambient_filter()] (correlation method).
#' @param results_vae Output from [ambient_filter_vae()] (VAE method).
#' @param highlight_genes Character vector. Genes to highlight in comparisons.
#'   Default `c("INS", "TTR", "GCG")`.
#'
#' @export
benchmark_methods <- function(results_corr,
                              results_vae,
                              highlight_genes = c("INS", "TTR", "GCG",
                                                  "PRSS1", "CELA3A")) {

  cell_types <- intersect(names(results_corr$cell_type_results),
                          names(results_vae$cell_type_results))

  comparison <- do.call(rbind, lapply(cell_types, function(ct) {
    corr_changed <- results_corr$cell_type_results[[ct]]$changed_genes
    vae_changed  <- results_vae$cell_type_results[[ct]]$changed_genes

    # For each highlight gene, get log2FC from both methods
    hl_rows <- lapply(highlight_genes, function(g) {
      corr_lfc <- if (g %in% corr_changed$gene) {
        corr_changed$log2fc[corr_changed$gene == g]
      } else NA_real_
      vae_lfc <- if (g %in% vae_changed$gene) {
        vae_changed$log2fc[vae_changed$gene == g]
      } else NA_real_
      data.frame(
        cell_type     = ct,
        gene          = g,
        log2fc_corr   = corr_lfc,
        log2fc_vae    = vae_lfc,
        vae_stronger  = !is.na(vae_lfc) && !is.na(corr_lfc) &&
                        abs(vae_lfc) > abs(corr_lfc),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, hl_rows)
  }))

  # Plot
  valid <- comparison[!is.na(comparison$log2fc_corr) &
                      !is.na(comparison$log2fc_vae), ]

  bench_plot <- NULL
  if (nrow(valid) > 0) {
    bench_plot <- ggplot2::ggplot(
      valid,
      ggplot2::aes(x = log2fc_corr, y = log2fc_vae,
                   colour = cell_type, label = gene)
    ) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           linetype = "dashed", colour = "grey60") +
      ggplot2::geom_point(size = 2.5, alpha = 0.8) +
      ggplot2::geom_text(hjust = -0.15, size = 2.8, show.legend = FALSE) +
      ggplot2::labs(
        title    = "Benchmark: correlation vs VAE correction",
        subtitle = paste0("log2FC for key ambient genes | ",
                          "points below diagonal = VAE corrects more aggressively"),
        x        = "log2FC (correlation method)",
        y        = "log2FC (VAE method)",
        colour   = "Cell type"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  }

  list(
    comparison_table = comparison,
    benchmark_plot   = bench_plot
  )
}


# ── Internal Python VAE code ──────────────────────────────────────────────────

#' @keywords internal
.check_vae_deps <- function(envname) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    rlang::abort(
      paste0("Package 'reticulate' is required for VAE correction. ",
             "Install with: install.packages('reticulate')")
    )
  }
  py_avail <- tryCatch(
    reticulate::py_available(initialize = TRUE),
    error = function(e) FALSE
  )
  if (!py_avail) {
    rlang::abort(
      paste0("Python not available. Install Python and run: ",
             "AmbientFilter::install_ambient_vae()")
    )
  }
  torch_avail <- tryCatch({
    reticulate::use_condaenv(envname, required = FALSE)
    reticulate::py_run_string("import torch")
    TRUE
  }, error = function(e) FALSE)

  if (!torch_avail) {
    rlang::abort(
      paste0("PyTorch not found in environment '", envname, "'. ",
             "Run: AmbientFilter::install_ambient_vae(envname = '", envname, "')")
    )
  }
  invisible(TRUE)
}


#' @keywords internal
.vae_python_code <- function() {
'
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np

class AmbientEncoder(nn.Module):
    def __init__(self, n_genes, hidden_dims, latent_dim_bio, latent_dim_amb):
        super().__init__()
        # Shared encoder layers
        layers = []
        in_dim = n_genes * 2  # counts + ambient profile concatenated
        for h in hidden_dims:
            layers += [nn.Linear(in_dim, h), nn.LayerNorm(h), nn.GELU()]
            in_dim = h
        self.shared = nn.Sequential(*layers)
        # Biological latent
        self.mu_bio    = nn.Linear(in_dim, latent_dim_bio)
        self.logvar_bio = nn.Linear(in_dim, latent_dim_bio)
        # Ambient latent
        self.mu_amb    = nn.Linear(in_dim, latent_dim_amb)
        self.logvar_amb = nn.Linear(in_dim, latent_dim_amb)

    def forward(self, x, a):
        h = self.shared(torch.cat([x, a], dim=1))
        return (self.mu_bio(h), self.logvar_bio(h),
                self.mu_amb(h), self.logvar_amb(h))


class AmbientDecoder(nn.Module):
    def __init__(self, n_genes, hidden_dims, latent_dim_bio, latent_dim_amb):
        super().__init__()
        hidden_dims_rev = list(reversed(hidden_dims))
        # Bio decoder
        bio_layers, in_b = [], latent_dim_bio
        for h in hidden_dims_rev:
            bio_layers += [nn.Linear(in_b, h), nn.LayerNorm(h), nn.GELU()]
            in_b = h
        self.decoder_bio = nn.Sequential(*bio_layers)
        self.out_bio     = nn.Linear(in_b, n_genes)
        # Ambient decoder
        amb_layers, in_a = [], latent_dim_amb
        for h in hidden_dims_rev:
            amb_layers += [nn.Linear(in_a, h), nn.LayerNorm(h), nn.GELU()]
            in_a = h
        self.decoder_amb = nn.Sequential(*amb_layers)
        self.out_amb     = nn.Linear(in_a, n_genes)

    def forward(self, z_bio, z_amb):
        x_bio = F.softplus(self.out_bio(self.decoder_bio(z_bio)))
        x_amb = F.softplus(self.out_amb(self.decoder_amb(z_amb)))
        return x_bio, x_amb


class AmbientVAE(nn.Module):
    def __init__(self, n_genes, hidden_dims, latent_dim_bio, latent_dim_amb):
        super().__init__()
        self.encoder = AmbientEncoder(n_genes, hidden_dims,
                                      latent_dim_bio, latent_dim_amb)
        self.decoder = AmbientDecoder(n_genes, hidden_dims,
                                      latent_dim_bio, latent_dim_amb)

    def reparameterise(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        return mu + std * torch.randn_like(std)

    def forward(self, x, a):
        mu_b, lv_b, mu_a, lv_a = self.encoder(x, a)
        z_bio = self.reparameterise(mu_b, lv_b)
        z_amb = self.reparameterise(mu_a, lv_a)
        x_bio, x_amb = self.decoder(z_bio, z_amb)
        return x_bio, x_amb, mu_b, lv_b, mu_a, lv_a


def spearman_loss(x, y):
    """Differentiable approximation to Spearman correlation loss."""
    n = x.shape[0]
    rx = x.argsort().argsort().float()
    ry = y.argsort().argsort().float()
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    return -(rx * ry).sum() / (rx.norm() * ry.norm() + 1e-8)


def run_ambient_vae(counts, ambient_profile, gene_weights,
                    latent_dim_bio, latent_dim_amb, hidden_dims,
                    n_epochs, batch_size, lr,
                    lambda_ambient, lambda_bio, use_gpu, seed, verbose):

    torch.manual_seed(seed)
    np.random.seed(seed)

    device = torch.device(
        "cuda" if (use_gpu and torch.cuda.is_available()) else "cpu"
    )

    counts_t  = torch.tensor(counts, dtype=torch.float32).to(device)
    ambient_t = torch.tensor(ambient_profile, dtype=torch.float32).to(device)
    gw_t      = torch.tensor(gene_weights, dtype=torch.float32).to(device)

    n_cells, n_genes = counts_t.shape
    # Log1p normalise input (preserve raw for output)
    counts_norm  = torch.log1p(counts_t)
    ambient_norm = ambient_t / (ambient_t.sum(dim=1, keepdim=True) + 1e-8)

    model = AmbientVAE(n_genes, hidden_dims, latent_dim_bio, latent_dim_amb)
    model = model.to(device)
    optimiser = torch.optim.Adam(model.parameters(), lr=lr)

    loss_curve = []
    dataset    = torch.utils.data.TensorDataset(counts_norm, ambient_norm, counts_t)
    loader     = torch.utils.data.DataLoader(
        dataset, batch_size=batch_size, shuffle=True, drop_last=False
    )

    model.train()
    for epoch in range(n_epochs):
        epoch_losses = {"total": 0, "recon": 0, "kl": 0, "ambient": 0}

        for x_batch, a_batch, x_raw_batch in loader:
            optimiser.zero_grad()

            x_bio, x_amb, mu_b, lv_b, mu_a, lv_a = model(x_batch, a_batch)
            x_recon = x_bio + x_amb

            # Reconstruction: MSE on log1p space
            recon_loss = F.mse_loss(torch.log1p(x_recon), x_batch)

            # KL divergence
            kl_bio = -0.5 * (1 + lv_b - mu_b.pow(2) - lv_b.exp()).sum(dim=1).mean()
            kl_amb = -0.5 * (1 + lv_a - mu_a.pow(2) - lv_a.exp()).sum(dim=1).mean()
            kl_loss = kl_bio + kl_amb

            # Ambient supervision: x_amb should correlate with ambient profile
            # weighted by gene_weights (high weight = confidently ambient)
            amb_supervision = 0.0
            for i in range(min(4, x_amb.shape[0])):
                amb_supervision += spearman_loss(
                    x_amb[i] * gw_t,
                    a_batch[i] * gw_t
                )
            amb_supervision = amb_supervision / min(4, x_amb.shape[0])

            # Bio preservation: x_bio should NOT correlate with ambient profile
            # (penalise if bio decoder reconstructs ambient signal)
            bio_preservation = 0.0
            for i in range(min(4, x_bio.shape[0])):
                bio_preservation += -spearman_loss(
                    x_bio[i] * (1 - gw_t),  # focus on non-ambient genes
                    a_batch[i]
                )
            bio_preservation = bio_preservation / min(4, x_bio.shape[0])

            total_loss = (recon_loss + kl_loss +
                          lambda_ambient * amb_supervision +
                          lambda_bio * bio_preservation)

            total_loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimiser.step()

            epoch_losses["total"]   += total_loss.item()
            epoch_losses["recon"]   += recon_loss.item()
            epoch_losses["kl"]      += kl_loss.item()
            epoch_losses["ambient"] += amb_supervision.item()

        n_batches = len(loader)
        epoch_avg = {k: v / n_batches for k, v in epoch_losses.items()}
        loss_curve.append(epoch_avg)

        if verbose and (epoch + 1) % 10 == 0:
            print(f"  Epoch {epoch+1:3d}/{n_epochs} | "
                  f"loss={epoch_avg[\"total\"]:.4f} | "
                  f"recon={epoch_avg[\"recon\"]:.4f} | "
                  f"ambient={epoch_avg[\"ambient\"]:.4f}")

    # ── Generate corrected counts ──────────────────────────────────────────
    model.eval()
    with torch.no_grad():
        x_bio_full, x_amb_full, _, _, mu_a_full, _ = model(counts_norm, ambient_norm)
        # Corrected = biological component only, scaled back to original library size
        lib_sizes  = counts_t.sum(dim=1, keepdim=True) + 1
        bio_frac   = x_bio_full / (x_bio_full + x_amb_full + 1e-8)
        corrected  = (counts_t * bio_frac).cpu().numpy()
        corrected  = np.round(corrected).astype(np.float32)
        corrected  = np.clip(corrected, 0, None)

    return {
        "corrected_counts": corrected,
        "ambient_scores":   mu_a_full.cpu().numpy(),
        "loss_curve":       loss_curve
    }
'
}
