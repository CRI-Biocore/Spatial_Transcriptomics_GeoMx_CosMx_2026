# GeoMx — Self-Study Guide

> Work through this guide at your own pace after the workshop. Each notebook is independent — pre-computed results are in `GeoMx/results/`, so you can open any notebook directly without having run the ones before it. The numbers suggest a logical progression, not a strict dependency.

---

## Dataset

The dataset is the same used for the GeomxTools tutorial and is used throughout all the scripts for this workshop. FFPE and FF tissue sections from 4 diabetic kidney disease (DKD) and 3 healthy kidney sections, processed using the GeoMx Digital Spatial Profiler (DSP) platform. The ROIs were selected to focus on tubules or glomeruli regions.

 - Glomeruli: Each ROI defined as glomeruli contains a single glomerulus.

 - Tubular: Each ROI contains multiple tubules and are further classified into distal (PanCK+) or proximal (PanCK-) AOIs.

For this workshop, we simulated 3 columns: The number of nuclei and X and Y coordinates per ROI, since the initial data did not have these fields.


---

## Analysis Workflow

Use the link to explore each step. The notebooks will show code, plots and notes to describe the process.

| Topic | Notebook |
|---|---|
| GeoMx obj from DCC, PKC data, QC and filtering with GeomxTools | [1_geomx_setup_qc.ipynb](GeoMx/codes/1_geomx_setup_qc.ipynb) |
| Normalization, PCA and DE analysis with limma-voom | [2_geomx_limma_norm_DE.ipynb](GeoMx/codes/2_geomx_limma_norm_DE.ipynb) |
| Cell-type deconvolution — SpatialDecon + Kidney HCA reference | [3_geomx_SpatialDecon.ipynb](GeoMx/codes/3_geomx_SpatialDecon.ipynb) |


---

## Notes

All the content for this section was created using the `geomx_env` (R-based) conda environment. All the analysis can be executed with **10G** of memory or less. If you want to try the code locally, this is possible. 


## Step 1 — Dataset Setup & QC

[1_geomx_setup_qc.ipynb](GeoMx/codes/1_geomx_setup_qc.ipynb)

### Context

GeoMx data is loaded from three file types: per-segment DCC expression files, a PKC probe annotation file, and an Excel segment annotation file. These are assembled into a `NanoStringGeoMxSet` object using `readNanoStringGeoMxSet()`. The notebook then runs segment-level QC, probe-level QC, computes the Limit of Quantification (LOQ) per segment, and filters both low-quality segments and lowly-detected genes before saving a tidy dataset for downstream analysis.

Key metadata columns kept in `meta_cols`: `Sample_ID`, `slide_name`, `region`, `segment`, `class`, `aoi`, `roi`, `area`, `nuclei`, `pathology`, `ROI_Coordinate_X`, `ROI_Coordinate_Y`.

### FAQS

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

### Reflect

> - Why is geometric mean used for the negative probe summary rather than arithmetic mean?
> - The LOQ is computed per segment — why does this matter compared to applying a single global expression cutoff?
> - How would you adjust the 10% gene detection threshold if your tissue had highly heterogeneous cell types where key markers are only expressed in a small fraction of segments?


---

## Step 2 — Normalization & Differential Expression

[2_geomx_limma_norm_DE.ipynb](codes/2_geomx_limma_norm_DE.ipynb)


### Context

Starting from the tidy TSVs produced in notebook 1, this notebook performs Q3 (upper-quartile) normalization via `edgeR::calcNormFactors()`, explores the dataset with PCA, and runs differential expression using `limma-voom` with `duplicateCorrelation` to account for the slide batch effect.

**DE question:** Within each kidney cell type, are there expression differences between DKD and normal samples?

The design matrix `~0 + CellType + class:CellType` produces interaction coefficients (e.g., `glomeruli_DKD`) representing the disease effect within each cell type. The notebook works through glomeruli as the primary example.

### Tasks

- [ ] Examine the PCA biplots colored by `slide_name` and by `CellType + class`. Do slides cluster separately (visible batch effect)? Do the three cell types separate on PC1/PC2? Check PC3/PC4 as well.
- [ ] Run `colnames(fit$coefficients)` and list all model terms. Identify which coefficient represents DKD vs. normal within glomeruli (`glomeruli_DKD`).
- [ ] Understand why `duplicateCorrelation` is used instead of adding `slide_name` directly to the design. What does `corfit$consensus.correlation` measure?  Check the value in this dataset.
- [ ] Compare the two heatmaps of the top 20 DEGs: raw `logCPM` vs. z-scored values. Which better reveals the expression pattern across samples?
- [ ] Results are saved to `GeoMx/results/DE_glomeruli_DKDvsNormal.tsv`. Update the code to run DE for `ProximalTubules_DKD` or `DistalTubules_DKD`.

### Reflect

> - Why is Q3/upper-quartile normalization more appropriate for GeoMx than library-size (CPM) normalization?
> - `duplicateCorrelation` treats slide as a random effect. When would you prefer a fixed-effect model? When would DESeq2 be the better choice?


---

## Step 3: Cell-Type Deconvolution with SpatialDecon

[3_geomx_SpatialDecon.ipynb](codes/3_geomx_SpatialDecon.ipynb)


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


### Reflect

> - How would results change if you used a custom reference matrix derived from a published scRNA-seq dataset of the same tissue instead of the pre-built atlas?
> - GeoMx deconvolution is entirely expression-based — it ignores the physical position of each segment. What spatial information does it miss that CosMx provides directly?
> - When would you trust a high estimated abundance for a rare cell type, and when would you be skeptical?

## Relevant Links

- [GeoMx Bioconductor Workflow](https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html)
- [SpatialDecon Bioconductor Vignette](https://bioconductor.org/packages/release/bioc/vignettes/SpatialDecon/inst/doc/SpatialDecon_vignette_NSCLC.html)