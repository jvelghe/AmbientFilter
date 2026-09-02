# AmbientFilter

Per-cell-type ambient RNA correction for single-cell RNA-seq, designed for
datasets with strong ambient signals (e.g. human pancreatic islets, where
INS and TTR contaminate all non-beta cell populations).

Built on [`DropletUtils::removeAmbience()`](https://bioconductor.org/packages/DropletUtils/),
extending the approach introduced by [dRopt](https://github.com/CUAnschutzBDC/dRopt)
to all cell types in a dataset with per-cell-type correction strength tuning,
tissue-stratified background estimation, data-driven correction strength
estimation, and a VAE-based disentanglement method for benchmarking.

> **Run upstream of subclustering.** Ambient correction should happen on your
> full annotated object before you subset and re-cluster individual populations.

---

## Why per-cell-type?

dRopt's insight is that a macrophage's ambient RNA profile should look like the
surrounding tissue — not like a single macrophage. The same is true for every
other cell type. A lymphocyte near an islet will have insulin mRNA on its
surface; that signal should be removed before subclustering lymphocytes.

But not all cell types should be corrected equally. **Macrophages are
biologically "nibbly"** — they internalise ambient RNA via efferocytosis and
trogocytosis as a genuine part of their biology. Correcting them at full
strength risks removing real signal. AmbientFilter handles this through
per-cell-type correction strengths, either set manually or learned from your
data automatically.

---

## The hybrid method

AmbientFilter implements a **two-layer hybrid approach** that combines the
strengths of per-sample and pooled correction:

**Layer 1 — Per-sample correction** uses `DropletUtils::removeAmbience()` with
a background profile estimated separately for each donor sample (and stratified
by tissue if you have multiple tissues). This captures the sample-specific
ambient soup composition, which varies by dissociation quality, cell viability,
and processing batch.

**Layer 2 — Cross-sample consistency analysis** asks, across all donors: which
genes are *consistently* corrected in a given cell type, regardless of
sample-specific variation? A gene that is reduced in 9/10 donors is almost
certainly genuine ambient contamination. A gene reduced in only 1/10 is likely
a processing artefact of that single sample. This pooled layer provides the
statistical power to build a robust, publishable ambient offenders list.

The output of layer 2 — the **ambient offenders table** — is saved as a
shareable CSV/RDS/Excel file that other labs working with the same tissue can
use to clean their own data.

---

## Correction strength: manual vs data-driven

The correction strength controls how aggressively ambient signal is subtracted
for a given cell type. AmbientFilter supports two modes:

### Manual (simple, interpretable)

```r
correction_strength = c(
  macrophage  = 0.5,   # conservative — macrophages are nibbly
  lymphocyte  = 1.0,
  ductal      = 1.0
)
```

| Cell type | Recommended | Why |
|-----------|-------------|-----|
| Lymphocytes, ductal, stellate, mast cells | `1.0` | No biological reason to contain ambient transcripts |
| Macrophages | `0.5` | Legitimately internalise ambient RNA via efferocytosis/trogocytosis |
| Beta cells, alpha cells | `0.0` or skip | They *are* the main ambient source |

### Data-driven (recommended, `auto_estimate_strength = TRUE`)

`estimate_correction_strength()` learns appropriate correction strengths from
the data itself by exploiting a key property of ambient contamination:
**ambient genes track the background profile across samples, while genuinely
expressed genes do not.**

For each cell type, across all donors, it computes the Spearman correlation
between each gene's mean expression in that cell type and its proportion in the
per-sample ambient background. A gene whose expression in macrophages rises and
falls in lockstep with the ambient INS proportion across donors is almost
certainly ambient. A gene whose expression is independent of the ambient profile
is almost certainly real.

The correction strength for a cell type is derived from the median of these
cross-sample correlations — a cell type where most expression tracks the
background gets an aggressive correction; one where expression is largely
independent of the background gets a conservative correction. Macrophages
receive a lower estimated strength automatically because their genuine marker
genes (TREM2, FOLR2, CD68) are uncorrelated with the INS/TTR ambient profile,
diluting the overall correlation signal.

You can blend the data-driven estimate with a biological prior:

```r
results <- ambient_filter(
  object                 = seu,
  auto_estimate_strength = TRUE,
  prior_strengths        = c(macrophage = 0.5),  # biological prior
  prior_weight           = 0.3   # 30% prior, 70% data-driven
)
```

`estimate_correction_strength()` also computes **per-gene ambient weights** for
each cell type — a value in [0, 1] reflecting how consistently each gene's
expression tracks the background across donors. These weights are the foundation
for the VAE method.

---

## VAE method

`ambient_filter_vae()` provides a more sophisticated correction using a
**variational autoencoder (VAE) with disentangled latent spaces**. Unlike the
correlation-based approach — which asks whether a gene's expression tracks the
background *across* samples — the VAE learns a **latent decomposition** of the
count matrix where ambient and genuine expression are separated as distinct
latent dimensions, conditioned on the known per-sample ambient profile.

### Architecture

```
Input: count matrix X (cells × genes) + ambient profile A (cells × genes)

  Encoder ──► z_bio     ~ N(μ_bio,    σ_bio)    # genuine expression
          └──► z_ambient ~ N(μ_ambient, σ_ambient)  # ambient contamination level

  Decoder_bio     (z_bio)     ──► X_bio     # genuine reconstruction
  Decoder_ambient (z_ambient) ──► X_ambient # ambient reconstruction

  X_reconstructed = X_bio + X_ambient

Loss = ELBO (NB reconstruction + KL divergence)
     + λ_ambient × corr(X_ambient, ambient_profile)   # ambient supervision
     + λ_bio     × (1 − corr(X_bio, ambient_profile)) # bio preservation
```

The ambient supervision loss uses a differentiable Spearman correlation
approximation to push `z_ambient` toward the background profile, while the
bio preservation loss penalises `z_bio` for capturing ambient variation.
Gene-level weights from `estimate_correction_strength()` focus both losses on
the genes where the evidence for ambient vs genuine expression is strongest.
The corrected counts are the original counts scaled by the fraction of
expression the model attributes to genuine biology: `X_corrected = X × (X_bio / (X_bio + X_ambient))`.

### When to use the VAE vs correlation method

| | Correlation method | VAE method |
|---|---|---|
| Speed | Fast (seconds–minutes) | Slower (minutes–hours depending on cell count and epochs) |
| Interpretability | High — based on Spearman correlation | Lower — latent space is implicit |
| Per-gene resolution | Via gene weights | Full gene-level disentanglement |
| Handles partial expression | Limited | Better — models the mixture explicitly |
| Requires Python | No | Yes (PyTorch via reticulate) |
| Best for | Quick QC, parameter setting, benchmarking baseline | Final publication-ready correction, large datasets |

The two methods are designed to be benchmarked against each other using
`benchmark_methods()`, which produces side-by-side comparison plots of key
ambient gene reductions across cell types.

### VAE installation and usage

```r
# Install Python dependencies (one time)
AmbientFilter::install_ambient_vae()

# Run VAE correction — accepts gene_weights from estimate_correction_strength()
# as soft supervision priors
strength_est <- estimate_correction_strength(seu, bg, ...)

results_vae <- ambient_filter_vae(
  object             = seu,
  background         = bg,
  cell_types         = c("macrophage", "lymphocyte", "ductal"),
  strength_estimates = strength_est,   # passes gene weights as VAE priors
  n_epochs           = 150L,
  latent_dim_bio     = 10L,
  latent_dim_ambient = 3L,
  save_dir           = "ambient_filter_vae_results"
)

# Benchmark correlation vs VAE
bench <- benchmark_methods(results_corr, results_vae,
                           highlight_genes = c("INS", "TTR", "GCG"))
bench$benchmark_plot
bench$comparison_table
```

---

## Installation

```r
# Bioconductor dependency
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DropletUtils")

# Install AmbientFilter from GitHub
devtools::install_github("jvelghe/AmbientFilter")

# For VAE method — install Python dependencies
AmbientFilter::install_ambient_vae()
```

---

## Quick start

```r
library(AmbientFilter)

# Recommended: data-driven correction strengths
results <- ambient_filter(
  object                    = seu,
  cell_type_col             = "cell_type",
  sample_col                = "hpap_id",
  tissue_col                = "tissue",
  cell_types                = c("macrophage", "lymphocyte", "mast_cell",
                                "ductal", "stellate", "endothelial"),
  auto_estimate_strength    = TRUE,          # learn strengths from data
  prior_strengths           = c(macrophage = 0.5),
  prior_weight              = 0.3,
  highlight_genes           = c("INS", "TTR", "GCG", "PRSS1"),
  run_consistency_analysis  = TRUE,          # cross-sample ambient offenders
  export_offenders          = TRUE,          # save shareable offenders table
  dataset_name              = "my_islet_atlas",
  tissue                    = "human_pancreatic_islet",
  make_plots                = TRUE,
  save_dir                  = "ambient_filter_results"
)
```

---

## Output

```r
# Summary: one row per cell type
results$summary

# Correction strengths learned from data
results$strength_estimates$summary
results$strength_estimates$plots$strength_bar

# Corrected macrophage object — ready for subclustering
mac <- results$cell_type_results$macrophage$corrected

# Changed genes per cell type
results$cell_type_results$lymphocyte$changed_genes

# Cross-sample ambient offenders (the publishable list)
results$ambient_map$cross_cell_type
results$ambient_map$plots$heatmap

# Diagnostic plots
results$plots$macrophage$correlation
results$plots$macrophage$violin_original
results$plots$macrophage$violin_corrected
```

---

## Testing

A functional test suite with synthetic islet scRNA-seq data is included.
It injects realistic ambient contamination (INS ~45%, TTR ~20% of ambient
reads) and validates that corrections reduce ambient signal while preserving
genuine marker gene expression.

```r
# From within R
source("tests/test_ambient_filter.R")

# From terminal
Rscript tests/test_ambient_filter.R
```

---

## Function reference

| Function | Description |
|----------|-------------|
| `ambient_filter()` | Main entry point — runs the full hybrid pipeline |
| `ambient_filter_vae()` | VAE-based correction pipeline |
| `estimate_background()` | Per-sample ambient background estimation |
| `estimate_correction_strength()` | Data-driven per-cell-type correction strength estimation |
| `correct_cell_type()` | Single cell type correction (correlation method) |
| `correct_cell_type_vae()` | Single cell type correction (VAE method) |
| `identify_consistent_ambient()` | Cross-sample ambient consistency analysis |
| `export_ambient_offenders()` | Export shareable ambient offenders table |
| `plot_correction_diagnostics()` | Diagnostic plots for correlation method |
| `plot_vae_diagnostics()` | Training loss and ambient score plots for VAE |
| `benchmark_methods()` | Side-by-side comparison of correlation vs VAE |
| `summarise_correction()` | Summary table across all cell types |
| `install_ambient_vae()` | Install Python/PyTorch dependencies for VAE |

---

## Acknowledgements

Inspired by [dRopt](https://github.com/CUAnschutzBDC/dRopt) (CU Anschutz BDC).
Background estimation and subtraction via [DropletUtils](https://bioconductor.org/packages/DropletUtils/).
VAE implementation uses [PyTorch](https://pytorch.org/) via
[reticulate](https://rstudio.github.io/reticulate/).
