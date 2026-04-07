# GeoMx — Self-Study Guide

> Work through this guide at your own pace after the workshop. Each notebook is independent — pre-computed results are in `GeoMx/results/`, so you can open any notebook directly without having run the ones before it. The numbers suggest a logical progression, not a strict dependency.

---

## Dataset

The dataset is the same used for the GeomxTools tutorial and is used throughout all the scripts for this workshop. FFPE and FF tissue sections from 4 diabetic kidney disease (DKD) and 3 healthy kidney sections, processed using the GeoMx Digital Spatial Profiler (DSP) platform. The ROIs were selected to focus on tubules or glomeruli regions.

 - Glomeruli: Each ROI defined as glomeruli contains a single glomerulus.

 - Tubular: Each ROI contains multiple tubules and are further classified into distal (PanCK+) or proximal (PanCK-) AOIs.



---

## Analysis Workflow

Use the link to explore each step. The notebooks will show code, plots and notes to describe the process.

| Topic | Notebook |
|---|---|
| GeoMx obj from DCC, PKC data, QC and filtering with GeomxTools | [1_geomx_setup_qc.ipynb](codes/1_geomx_setup_qc.ipynb) |
| Normalization, PCA and DE analysis with limma-voom | [2_geomx_limma_norm_DE.ipynb](codes/2_geomx_limma_norm_DE.ipynb) |
| Cell-type deconvolution — SpatialDecon + Kidney HCA reference | [3_geomx_SpatialDecon.ipynb](codes/3_geomx_SpatialDecon.ipynb) |


---

## Notes

All the content for this section was created using the `geomx_env` (R-based) conda environment. All the analysis can be executed with **10G** of memory or less. If you want to try the code locally, this is possible. 


## Step 1 — Dataset Setup & QC

[1_geomx_setup_qc.ipynb](codes/1_geomx_setup_qc.ipynb)

### Context

GeoMx data is loaded from three file types: per-segment DCC expression files, a PKC probe annotation file, and an Excel segment annotation file. These are assembled into a `NanoStringGeoMxSet` object using `readNanoStringGeoMxSet()`. The notebook then runs segment-level QC, probe-level QC, computes the Limit of Quantification (LOQ) per segment, and filters both low-quality segments and lowly-detected genes before saving a tidy dataset for downstream analysis.

Key metadata columns kept in `meta_cols`: `Sample_ID`, `slide_name`, `region`, `segment`, `class`, `aoi`, `roi`, `area`, `nuclei`, `pathology`, `ROI_Coordinate_X`, `ROI_Coordinate_Y`.


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

The design matrix `~ 0 + CellType + class:CellType` produces interaction coefficients (e.g., `glomeruli_DKD`) representing the disease effect within each cell type. The notebook works through glomeruli as the primary example.

### Tasks

- [ ] Examine the PCA biplots colored by `slide_name` and by `CellType + class`. Do slides cluster separately (visible batch effect)? Do the three cell types separate on PC1/PC2? Check PC3/PC4 as well.
- [ ] Compare the two heatmaps of the top 20 DEGs: raw `logCPM` vs. z-scored values. Which better reveals the expression pattern across samples?


### Reflect

> - Why is Q3/upper-quartile normalization more appropriate for GeoMx than library-size (CPM) normalization?
> - `duplicateCorrelation` treats slide as a random effect. When would you prefer a fixed-effect model? When would DESeq2 be the better choice?


---

## Step 3: Cell-Type Deconvolution with SpatialDecon

[3_geomx_SpatialDecon.ipynb](codes/3_geomx_SpatialDecon.ipynb)


### Context

Because GeoMx segments contain mixtures of cells, cell-type composition must be estimated by deconvolving bulk-like expression profiles. This notebook uses the `SpatialDecon` package with the pre-built **Human Kidney Atlas** (`Kidney_HCA`) reference matrix from the Nanostring CellProfileLibrary. Background is estimated from `NegProbe-WTX` counts. 

### Tasks

- [ ] Walk through Q3 normalization applied to the `NanoStringGeoMxSet` object (stored as the `q_norm` assay). This normalized matrix along with raw counts is the input to deconvolution.
- [ ] Inspect the `Kidney_HCA` reference matrix (genes × cell types) via its `Heatmap()`. Which cell types have the most distinct expression signatures? How many cell types does it include?
- [ ] Compare the two deconvolution calls in the notebook:
  - `runspatialdecon()` — wrapper taking the full `NanoStringGeoMxSet` object
  - `spatialdecon()` — direct call with normalized matrix, background, and `cell_counts = sData(geomx)$nuclei`
- [ ] In the `ComplexHeatmap`, columns are split by `celltype` (Glomeruli, Proximal Tubules, Distal Tubules) and annotated by `class` (DKD vs. normal). Do any cell types shift in abundance between disease and healthy samples?


## Relevant Links

- [GeoMx Bioconductor Workflow](https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html)
- [SpatialDecon Bioconductor Vignette](https://bioconductor.org/packages/release/bioc/vignettes/SpatialDecon/inst/doc/SpatialDecon_vignette_NSCLC.html)