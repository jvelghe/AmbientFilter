# AmbientFilter

Per-cell-type ambient RNA correction for single-cell RNA-seq, designed for
datasets with strong ambient signals (e.g. human pancreatic islets, where
INS and TTR contaminate all non-beta cell populations).

Built on [`DropletUtils::removeAmbience()`](https://bioconductor.org/packages/DropletUtils/),
extending the approach introduced by [dRopt](https://github.com/CUAnschutzBDC/dRopt)
to all cell types in a dataset with per-cell-type correction strength tuning
and tissue-stratified background estimation.

> Run **upstream of subclustering**. Ambient correction should happen on your
> full annotated object before you subset and re-cluster individual populations.

---

## Why per-cell-type?

dRopt's insight is that a macrophage's ambient RNA profile should look like the
surrounding tissue — not like a single macrophage. The same is true for every
other cell type. A lymphocyte near an islet will have insulin mRNA on its
surface; that signal should be removed before subclustering lymphocytes.

The key extension here is **correction strength per cell type**:

| Cell type | Recommended strength | Why |
|-----------|---------------------|-----|
| Lymphocytes, ductal, stellate, mast cells | `1.0` | No biological reason to contain ambient transcripts |
| Macrophages | `0.5` | Legitimately internalise ambient RNA via efferocytosis/trogocytosis |
| Beta cells, alpha cells | `0.0` or skip | They *are* the main ambient source — correcting them removes real signal |

---

## Installation

```r
# Install dependencies first
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DropletUtils")

# Install from GitHub
devtools::install_github("jvelghe/AmbientFilter")
```

---

## Quick start

```r
library(AmbientFilter)

results <- ambient_filter(
  object         = seu,
  cell_type_col  = "cell_type",
  sample_col     = "hpap_id",
  tissue_col     = "tissue",          # stratify background per tissue
  cell_types     = c("macrophage", "lymphocyte", "mast_cell",
                     "ductal", "stellate", "endothelial"),
  correction_strength = c(
    macrophage   = 0.5,               # conservative — macrophages are nibbly
    lymphocyte   = 1.0,
    mast_cell    = 1.0,
    ductal       = 1.0,
    stellate     = 1.0,
    endothelial  = 1.0
  ),
  highlight_genes  = c("INS", "TTR", "GCG", "PRSS1"),
  make_plots       = TRUE,
  save_dir         = "ambient_filter_results"
)
```

---

## Output

```r
# Summary table — one row per cell type
results$summary

# Corrected macrophage object — ready for subclustering
mac <- results$cell_type_results$macrophage$corrected

# Changed genes per cell type
results$cell_type_results$lymphocyte$changed_genes

# Diagnostic plots
results$plots$macrophage$correlation
results$plots$macrophage$violin_original
results$plots$macrophage$violin_corrected

# Full corrected Seurat object (only if return_full_corrected_obj = TRUE)
results$updated_object
```

---

## Function reference

| Function | Description |
|----------|-------------|
| `ambient_filter()` | Main entry point — runs the full pipeline |
| `estimate_background()` | Per-sample background estimation |
| `correct_cell_type()` | Single cell type correction |
| `plot_correction_diagnostics()` | Diagnostic plots for one cell type |
| `summarise_correction()` | Summary table across all cell types |

---

## Acknowledgements

Inspired by [dRopt](https://github.com/CUAnschutzBDC/dRopt) (CU Anschutz BDC).
Background estimation and subtraction via [DropletUtils](https://bioconductor.org/packages/DropletUtils/).
