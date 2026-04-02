# GeoMx — Self-Study Guide

> Work through this guide at your own pace after the workshop. Each notebook is independent — pre-computed results are in `GeoMx/results/`, so you can open any notebook directly without having run the ones before it. The numbers suggest a logical progression, not a strict dependency.

---

## Dataset

> 📌 **Placeholder** — Replace with your dataset description.
>
> *Example fields to include:*
> - Source / publication
> - Tissue type and condition (e.g., disease vs. healthy)
> - Number of slides, ROIs/segments, AOI types
> - Panel used (PKC file name, approximate gene count)
> - Any preprocessing or simulation notes

---

## Quick Reference

| Topic | Notebook |
|---|---|
| Data loading — DCC, PKC, annotation Excel | `codes/1_geomx_setup_qc.ipynb` |
| Segment QC — reads, saturation, nuclei, area, NTC | `codes/1_geomx_setup_qc.ipynb` |
| Probe QC — low ratio and Grubbs outlier filtering | `codes/1_geomx_setup_qc.ipynb` |
| Limit of Quantification (LOQ) and gene-level filtering | `codes/1_geomx_setup_qc.ipynb` |
| Q3 normalization and PCA exploration | `codes/2_geomx_limma_norm_DE.ipynb` |
| Differential expression — limma-voom + duplicateCorrelation | `codes/2_geomx_limma_norm_DE.ipynb` |
| Volcano plot and DEG heatmaps | `codes/2_geomx_limma_norm_DE.ipynb` |
| Background estimation from negative probes | `codes/3_geomx_SpatialDecon.ipynb` |
| Cell-type deconvolution — SpatialDecon + reference matrix | `codes/3_geomx_SpatialDecon.ipynb` |
| Cell abundance visualization — barplot and heatmap | `codes/3_geomx_SpatialDecon.ipynb` |

---

## Environment & Launch

**Environment:** `geomx_env` (R-based)

```bash
sbatch --mem=50GB ./launch_jupyter_conda.sh \
  /gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2/conda_envs/geomx_env
```

Launch Jupyter from the **workshop root folder** — all paths inside the notebooks are relative to that location.

---

## Notebook 1 — Dataset Setup & QC

**`GeoMx/codes/1_geomx_setup_qc.ipynb`**

### Context

GeoMx data is loaded from three file types: per-segment DCC expression files, a PKC probe annotation file, and an Excel segment annotation file. These are assembled into a `NanoStringGeoMxSet` object using `readNanoStringGeoMxSet()`. The notebook then runs segment-level QC, probe-level QC, computes the Limit of Quantification (LOQ) per segment, and filters both low-quality segments and lowly-detected genes before saving a tidy dataset for downstream analysis.

Key metadata columns kept in `meta_cols`: `Sample_ID`, `slide_name`, `region`, `segment`, `class`, `aoi`, `roi`, `area`, `nuclei`, `pathology`, `ROI_Coordinate_X`, `ROI_Coordinate_Y`.

### Tasks

- [ ] Load the object and explore its structure. Understand the difference between `pData()` (segment metadata), `sData()` (extended protocol metadata including QC flags), and `fData()` (probe/gene metadata). What does each contain?
- [ ] Review the segment QC parameters in `QC_params` and compare them to the commented defaults:
  - `minSegmentReads = 1000`, `percentTrimmed/Stitched/Aligned = 80%`, `percentSaturation = 50%`
  - `minNegativeCount = 1.2` (default: 10 — deliberately relaxed here)
  - `minNuclei = 50` (default: 100 — relaxed), `minArea = 5000`, `maxNTCCount = 10000`

  Why might relaxing `minNuclei` and `minNegativeCount` be appropriate for this dataset?
- [ ] Walk through each QC histogram, faceted by `segment` type. Do trimmed, stitched, and aligned read percentages differ across segment types? Does the NTC distribution follow the expected pattern?
- [ ] Understand how NTC and negative probe percentage are computed from the expression matrix and stored in `protocolData`. Why is the NegGeoMean removed from `pData` before calling `aggregateCounts()`?
- [ ] After segment filtering, note how the dataset goes from 235 → 211 segments. Check the QC flag table written to `results/segment_metadata_QC.tsv` — the `Flag` column explains each failure reason.
- [ ] Walk through probe QC with `setBioProbeQCFlags()` using `minProbeRatio = 0.1` and `percentFailGrubbs = 20`. Two flags are applied: `LowProbeRatio` and `GlobalGrubbsOutlier`. How many probes does each flag catch?
- [ ] Follow the LOQ calculation:

  $$LOQ = \text{GeoMean(NegProbe)} \times \text{GeoSD(NegProbe)}^2 \quad (\text{minimum} = 2)$$

  The resulting `LOQ_Mat` (genes × segments boolean matrix) records whether each gene exceeds the per-segment noise floor. How is it used to define gene detection rates?
- [ ] Segments with < 10% genes detected above LOQ are removed. Then genes detected in < 10% of remaining segments are filtered out, while `NegProbe-WTX` is explicitly retained. Confirm the final dimensions: the notebook targets ~10,028 genes and 201 segments.
- [ ] Note all outputs in `GeoMx/results/`: `segment_metadata_QC.tsv`, `raw_counts.tsv`, `raw_features.tsv`, `raw_metadata.tsv`, `tidy_counts.tsv`, `tidy_features.tsv`, `tidy_metadata.tsv`, `tidy_geomx_obj.rds`, `tidy_spe_obj.rds`.

### Reflect

> - Why is geometric mean used for the negative probe summary rather than arithmetic mean?
> - The LOQ is computed per segment — why does this matter compared to applying a single global expression cutoff?
> - How would you adjust the 10% gene detection threshold if your tissue had highly heterogeneous cell types where key markers are only expressed in a small fraction of segments?

### 💡 Inspect Instead

```r
library(GeomxTools); library(tidyverse)

target_data <- readRDS("GeoMx/results/tidy_geomx_obj.rds")
dim(target_data)
pData(target_data) %>% select(region, segment, class, GeneDetectionRate) %>% head()

# QC summary
read_tsv("GeoMx/results/segment_metadata_QC.tsv") %>% count(QCStatus, Flag)
```

---

## Notebook 2 — Normalization & Differential Expression

**`GeoMx/codes/2_geomx_limma_norm_DE.ipynb`**

### Context

Starting from the tidy TSVs produced in notebook 1, this notebook performs Q3 (upper-quartile) normalization via `edgeR::calcNormFactors()`, explores the dataset with PCA, and runs differential expression using `limma-voom` with `duplicateCorrelation` to account for the slide batch effect.

**DE question:** Within each kidney cell type, are there expression differences between DKD and normal samples?

The design matrix `~0 + CellType + class:CellType` produces interaction coefficients (e.g., `glomeruli_DKD`) representing the disease effect within each cell type. The notebook works through glomeruli as the primary example.

### Tasks

- [ ] Load the counts and metadata from `GeoMx/results/`. Note how `CellType` is recoded from `segment`: `Geometric Segment` → `glomeruli`, `PanCK−` → `DistalTubules`, `PanCK+` → `ProximalTubules`. Why recode rather than use the raw column?
- [ ] Examine the PCA biplots colored by `slide_name` and by `CellType + class`. Do slides cluster separately (visible batch effect)? Do the three cell types separate on PC1/PC2? Check PC3/PC4 as well.
- [ ] Run `colnames(fit$coefficients)` and list all model terms. Identify which coefficient represents DKD vs. normal within glomeruli (`glomeruli_DKD`).
- [ ] Understand why `duplicateCorrelation` is used instead of adding `slide_name` directly to the design. What does `corfit$consensus.correlation` measure? Values > 0.5 would call for a second round on residuals — check the actual value in this dataset.
- [ ] Run `topTable()` for `glomeruli_DKD` with `number=Inf`. Apply cutoffs `adj.P.Val < 0.05` and `|logFC| ≥ log2(1.5)`. How many genes pass? Do any recognized disease-related genes appear in the top hits?
- [ ] Interpret the volcano plot: dashed lines mark `logFC = ±0.59` (1.5-fold) and `-log10(adj.P.Val) = 1.3` (FDR 5%). Significant genes are labeled using `ggrepel`.
- [ ] Compare the two heatmaps of the top 20 DEGs: raw `logCPM` vs. z-scored values. Which better reveals the expression pattern across samples?
- [ ] Results are saved to `GeoMx/results/DE_glomeruli_DKDvsNormal.tsv`. Identify the three lines you would change to run DE for `ProximalTubules_DKD` or `DistalTubules_DKD`.

### Reflect

> - Why is Q3/upper-quartile normalization more appropriate for GeoMx than library-size (CPM) normalization?
> - `duplicateCorrelation` treats slide as a random effect. When would you prefer a fixed-effect model? When would DESeq2 be the better choice?
> - The design uses `~0 + CellType + class:CellType`. What is the advantage of this parameterization when you want to test the disease effect within one cell type at a time?
> - If a gene is significantly DE in glomeruli but not in tubules, does that prove it is unaffected in tubules?

### 💡 Inspect Instead

```r
library(tidyverse)

de <- read_tsv("GeoMx/results/DE_glomeruli_DKDvsNormal.tsv")

de %>%
  filter(adj.P.Val < 0.05, abs(logFC) >= log2(1.5)) %>%
  arrange(desc(abs(logFC))) %>%
  select(Gene, logFC, adj.P.Val) %>%
  head(20)
```

---

## Notebook 3 — Cell-Type Deconvolution with SpatialDecon

**`GeoMx/codes/3_geomx_SpatialDecon.ipynb`**

### Context

Because GeoMx segments contain mixtures of cells, cell-type composition must be estimated by deconvolving bulk-like expression profiles. This notebook uses the `SpatialDecon` package with the pre-built **Human Kidney Atlas** (`Kidney_HCA`) reference matrix from the Nanostring CellProfileLibrary. Background is estimated from `NegProbe-WTX` counts. Nuclei counts (simulated for this dataset) are optionally passed to scale the deconvolution output.

### Tasks

- [ ] Walk through Q3 normalization applied to the `NanoStringGeoMxSet` object (stored as the `q_norm` assay). This normalized matrix — not raw counts — is the input to deconvolution.
- [ ] Follow the background calculation with `derive_GeoMx_background()`, using `NegProbe-WTX` to estimate per-segment, per-module background. Why is this step described as optional?
- [ ] Inspect the `Kidney_HCA` reference matrix (genes × cell types) via its `Heatmap()`. Which cell types have the most distinct expression signatures? How many cell types does it include?
- [ ] Compare the two deconvolution calls in the notebook:
  - `runspatialdecon()` — wrapper taking the full `NanoStringGeoMxSet` object
  - `spatialdecon()` — direct call with normalized matrix, background, and `cell_counts = sData(geomx)$nuclei`

  How does including nuclei counts change the `beta` output? What do the rows and columns of `beta` represent?
- [ ] Examine the `TIL_barplot()` output. Which cell types dominate each AOI type? Are proportions consistent within the same compartment?
- [ ] In the `ComplexHeatmap`, columns are split by `celltype` (Glomeruli, Proximal Tubules, Distal Tubules) and annotated by `class` (DKD vs. normal). Do any cell types shift in abundance between disease and healthy samples?
- [ ] Connect back to notebook 2: do the top DE genes in glomeruli correspond to known markers of the cell types showing the largest DKD-related abundance shifts?

### Reflect

> - SpatialDecon uses non-negative least squares (NNLS). What are the key assumptions, and in what situations would it produce unreliable estimates?
> - How would results change if you used a custom reference matrix derived from a published scRNA-seq dataset of the same tissue instead of the pre-built atlas?
> - GeoMx deconvolution is entirely expression-based — it ignores the physical position of each segment. What spatial information does it miss that CosMx provides directly?
> - When would you trust a high estimated abundance for a rare cell type, and when would you be skeptical?

### 💡 Inspect Instead

```r
library(GeomxTools); library(SpatialDecon); library(tidyverse)

geomx <- readRDS("GeoMx/results/tidy_geomx_obj.rds")
geomx <- normalize(geomx, fromElt = "exprs", norm_method = "quant",
                   desiredQuantile = 0.75, toElt = "q_norm")
# From here, run derive_GeoMx_background() and spatialdecon() as in the notebook,
# or load a saved beta matrix from GeoMx/results/ if available.
```

---

## Connecting the Analyses

Once you have worked through all three notebooks, consider these integrative questions:

> - The segment QC in notebook 1 uses nuclei count and area — how would you expect these to correlate with the cell-type proportions estimated in notebook 3?
> - The DE analysis in notebook 2 identifies genes changing between DKD and normal. The deconvolution in notebook 3 identifies cell types changing in abundance. These are complementary views — how would you distinguish a gene changing because of cell-type composition shift vs. a genuine transcriptional response within a stable cell population?
> - If you could add one more analysis step after these three notebooks, what would it be and why?
