# CosMx — Post-Class Self-Study Guide

> Work through this guide at your own pace after the workshop. Each notebook is independent — pre-computed results are in `CosMx/results/`, so you can open any notebook directly without having run the ones before it. The numbers suggest a logical progression, not a strict dependency.

---

## Dataset

For this section we use a dataset from FFPE human breast tissue, tested with the Whole-transcriptome CosMx panel as well as a 64 proteomics panel. The dataset was downloaded from [Bruker, multiomic data](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/cosmx-human-multiomic-breast-dataset/)
For this workshop, we did focus only in the transcriptomics panel results. 


---

## Quick Reference

This folder includes key steps from a spatial sc-RNA-seq dataset. 

| Topic | Notebook / Resource |
|---|---|
| AnnData obj from flat files, QC and filtering | `CosMx/codes/1_cosmx_preprocessing_qc.ipynb` |
| Obj subsetting an MTX export | `CosMx/codes/Subset_dataset.ipynb` |
| Normalization, dimension reduction, clustering and annotation  | `CosMx/codes/2_cosmx_cell_annotation.ipynb` |
| Cell-cell communication with CellPhoneDB v5 | `CosMx/codes/3_cosmx_cell_comunication.ipynb` |
| Spatially-aware Cell-cell communication with SpatialCellChat — compute (heavy, HPC script) | `CosMx/codes/4_cosmx_spatial_cell_com.R` |
| Spatially-aware Cell-cell communication with SpatialCellChat — notebook with visualizations | `CosMx/codes/4_cosmx_spatial_cell_com.ipynb` |


---

## Environments & Launch

| Notebook | Environment | Memory |
|---|---|---|
| Notebook 1, Subset | `scverse_spatial` (Python) | 100–180 GB |
| Notebooks 2–3 | `scverse_spatial` (Python) | 50 GB |
| Notebook 4 (visualization) | `seurat_spatial` (R) | 50 GB |
| Notebook 4 (compute, `.R` script) | `seurat_spatial` (R) | 180 GB, interactive node |



Launch from the **workshop root folder** — all paths inside the notebooks are relative to that location.

---

## Notebook 1 — Ingestion, Preprocessing & QC

**`CosMx/codes/1_cosmx_preprocessing_qc.ipynb`** · `scverse_spatial` · ~180 GB

### Context

CosMx flat files are loaded manually: the expression matrix (`*exprMat_file.csv.gz`) is read with `pyarrow` and immediately converted to a sparse `csr_matrix` to halve peak RAM usage, then combined with cell metadata (`*metadata_file.csv.gz`) into an `AnnData` object. Global pixel coordinates (`CenterX_global_px`, `CenterY_global_px`) are stored in `obsm["spatial"]` for `squidpy` compatibility. QC uses negative probe percentage as the primary noise metric alongside transcript counts, gene counts, and cell area.

### Tasks

- [ ] Trace the memory-efficient loading strategy: the dense expression DataFrame is converted to sparse format and deleted (`del expr; gc.collect()`) *before* loading the metadata file. Why does this order matter for a dataset of this size (384 FOVs)?
- [ ] Verify that `obsm["spatial"]` stores global pixel coordinates, not FOV-local ones. Print the coordinate range — what are the approximate x and y extents of the slide?
- [ ] Identify the negative probe columns (names contain `"NegPrb"`). `sc.pp.calculate_qc_metrics()` is called with `qc_vars=["negative_probe"]`, adding `pct_counts_negative_probe` to `adata.obs`. Why is this metric preferred over mitochondrial gene percentage for CosMx data?
- [ ] Review the FOV-level QC plots (average negative probe signal and average counts per FOV). What is the purpose of this step before applying cell-level filters?
- [ ] Review the filtering thresholds and the per-filter removal counts printed before applying the combined mask:
  - `MIN_COUNTS = 200`, `MAX_COUNTS = 7500`
  - `MIN_GENES = 100`
  - `MAX_NEG_PCT = 1.0%`
  - `MIN_AREA = 50`, `MAX_AREA = 15000`

  Which individual filter removes the most cells? What does `MAX_COUNTS = 7500` guard against?
- [ ] Compare the spatial scatter plots of `total_counts` and `pct_counts_negative_probe` before and after filtering. Is the noise spatially random, or concentrated in particular regions or FOVs?
- [ ] After filtering, negative probe columns are removed from the feature set. What are the final cell and gene counts saved to `CosMx/results/BreastCancer_filtered.h5ad`?

### Reflect

> - Why is `pct_counts_negative_probe` a better noise metric here than mitochondrial gene percentage?
> - What biological or tissue-specific considerations should guide `MIN_AREA` and `MAX_AREA`? How would these change for a tissue with very large cells (e.g., neurons)?
> - Why store global pixel coordinates in `obsm["spatial"]` rather than FOV-local coordinates?
> - How would you adjust `MIN_COUNTS` for a smaller targeted panel (e.g., 1,000 genes) vs. this larger panel?

### 💡 Inspect Instead

```python
import scanpy as sc
adata = sc.read_h5ad("CosMx/results/BreastCancer_filtered.h5ad")
print(adata)
adata.obs[["total_counts", "n_genes_by_counts",
           "pct_counts_negative_probe", "Area"]].describe()
```

---

## Notebook 1b — Subsetting the Dataset

**`CosMx/codes/Subset_dataset.ipynb`** · `scverse_spatial` · ~100 GB

> ⚠️ **Runs after notebook 1.** Reads `BreastCancer_filtered.h5ad` and produces the subset used by all downstream notebooks. Subset files are already in `CosMx/results/` — you do not need to rerun this unless you want to define a different region.

### Context

The full filtered dataset (384 FOVs) is too large for interactive annotation and communication analysis. This notebook selects a spatial rectangle targeting ~20,000–30,000 cells by identifying whole FOVs within defined coordinate bounds. The subset is saved as an `h5ad` file and also exported in 10x-compatible MTX format for loading into Seurat in the SpatialCellChat workflow.

### Tasks

- [ ] Understand the subsetting strategy: whole FOVs within a spatial rectangle are selected rather than randomly sampling individual cells. Why is FOV-based selection preferable to random cell sampling for spatial analysis?
- [ ] The coordinate bounds select the middle third in x and a central band in y:
  ```python
  x: [x_min + 1/3*(x_max-x_min),  x_min + 2/3*(x_max-x_min)]
  y: [y_min + 2/9*(y_max-y_min),  y_min + 6/9*(y_max-y_min)]
  ```
  Look at the spatial scatter plot of retained FOVs. Does the selected region cover a representative area of the tissue?
- [ ] Check `subset_adata` after subsetting: how many cells and FOVs are retained? Is it within the 20,000–30,000 cell target?
- [ ] Four files are written to `CosMx/results/Subset_MM/`: `matrix.mtx`, `barcodes.tsv`, `features.tsv`, `metadata.csv`. The count matrix is transposed (`.X.T`) for Seurat compatibility — why does Seurat expect genes × cells rather than cells × genes?

### 💡 Inspect Instead

```python
import scanpy as sc
subset = sc.read_h5ad("CosMx/results/BreastCancer_subset.h5ad")
print(f"Cells : {subset.n_obs:,}")
print(f"Genes : {subset.n_vars:,}")
print(f"FOVs  : {subset.obs['fov'].nunique()}")
```

---

## Notebook 2 — Clustering & Cell Type Annotation

**`CosMx/codes/2_cosmx_cell_annotation.ipynb`** · `scverse_spatial` · ~50 GB

### Context

Starting from `BreastCancer_subset.h5ad` (31,999 cells), this notebook follows a standard scanpy workflow: normalization to the median library size, HVG selection (top 3,000 genes, `seurat_v3` flavor on raw counts), PCA (50 components, `arpack` solver), neighbor graph (15 neighbors, 17 PCs), UMAP, and Leiden clustering at resolutions 0.2, 0.4, and 0.6.

Annotation uses **CellTypist** with three pre-downloaded models from `data/celltypist_models/`:
- `Cells_Adult_Breast.pkl` — breast-specific cell types
- `Immune_All_Low.pkl` — broad immune categories
- `Immune_All_High.pkl` — fine-grained immune types

All three are run with `majority_voting=True` over `leiden_0.4` clusters. A consensus annotation selects the best model per cluster based on mean prediction confidence.

### Tasks

- [ ] Review the normalization: `sc.pp.normalize_total()` uses `target_sum = median_counts` (the subset's median library size) rather than the conventional 10,000. Why is the median library size a better target for this probe-based dataset?
- [ ] Note that HVGs are computed with `flavor="seurat_v3"` on the raw `"counts"` layer. Why does this flavor specifically require raw counts?
- [ ] Examine the elbow plot from `sc.pl.pca_variance_ratio()`. The notebook uses `n_pcs = 17` — where is the elbow, and does the choice seem well-supported?
- [ ] Look at the UMAP colored by `Mean.PanCK`, `Mean.CD68`, `Mean.CD45`, and `Mean.DAPI` (morphological markers from segmentation). Which populations do these highlight before any clustering label is applied?
- [ ] Compare spatial scatter plots for Leiden resolutions 0.2, 0.4, and 0.6. Which resolution produces the most spatially coherent clusters? The notebook selects `leiden_0.4` — do you agree?
- [ ] Review the `rank_genes_groups` dotplot for `leiden_0.4` (Wilcoxon test on raw counts). Do the top 5 markers per cluster give an intuition for cell identity before running CellTypist?
- [ ] Walk through the three CellTypist runs. Prediction probability maxima per cell are stored as confidence scores. Why track confidence alongside predicted labels?
- [ ] Inspect `consensus_df` and `conf_summary`: most-frequent label and mean confidence per cluster per model. The final annotation uses `Immune_All_High` labels for most clusters but switches to `Cells_Adult_Breast` for clusters 0, 2, 3, and 7 (luminal epithelial). Why are two models needed?
- [ ] Review the luminal marker dotplot (`ESR1`, `PGR`, `FOXA1`, `KRT8`, `MKI67`, `TOP2A`, `KRT5`, `KRT14`, etc.) used to characterize epithelial clusters before finalizing annotation.
- [ ] Verify the `final_annotation` spatial scatter plot. Do luminal, immune, and stromal populations occupy distinct tissue regions?
- [ ] Check outputs saved to `CosMx/results/`: `BreastCancer_subset_annotated.h5ad` and the `Subset_MM/` tables (`spatial_coordinates.csv`, `pca_embeddings.csv`, `umap_embeddings.csv`, `metadata_annotation.csv`).

### Reflect

> - The notebook uses three CellTypist models and builds a consensus. What are the risks of using only one model for a complex tissue like breast cancer?
> - Majority voting at the cluster level differs from per-cell annotation. When would per-cell annotation be preferable, and when would it be less reliable?
> - The luminal marker dotplot suggests clusters 0, 2, 3, and 7 could be further subdivided. What biological distinction would you look for to justify splitting them?
> - Does the spatial organization of `final_annotation` support the expression-derived clustering, or do you see unexpected mismatches?

### 💡 Inspect Instead

```python
import scanpy as sc, squidpy as sq

adata = sc.read_h5ad("CosMx/results/BreastCancer_subset_annotated.h5ad")
print(adata.obs["final_annotation"].value_counts())
sq.pl.spatial_scatter(adata, color="final_annotation", shape=None, size=0.5)
```

---

## Notebook 3 — Cell–Cell Communication (Expression-Based)

**`CosMx/codes/3_cosmx_cell_comunication.ipynb`** · `scverse_spatial` · ~100 GB

### Context

This notebook uses **CellPhoneDB v5** (database at `data/cellphone_db_v5/cellphonedb.zip`, pre-downloaded) to infer ligand–receptor interactions from expression alone — without spatial constraints. CellPhoneDB requires two inputs: a lognorm expression `h5ad` (with the `lognorm` layer promoted to `X`) and a metadata TSV with columns `barcode_sample` and `cell_type`. Both the basic analysis method and the statistical permutation-based method are run. Results are visualized with `ktplotspy`.

This analysis intentionally ignores spatial coordinates — it establishes a baseline of putative interactions to compare against the proximity-constrained SpatialCellChat results in notebook 4.

### Tasks

- [ ] Follow the data preparation: `adata.layers["lognorm"]` is promoted to `X` in a new `AnnData` before saving as the CellPhoneDB counts input. Why does CellPhoneDB require lognorm rather than raw counts?
- [ ] Check the metadata TSV: exactly two columns — `barcode_sample` (cell ID) and `cell_type` (`final_annotation`). CellPhoneDB groups cells by `cell_type` when computing mean expression per pair.
- [ ] The **basic method** (`cpdb_analysis_method.call()`) runs without permutation testing. Key parameters: `threshold = 0.1`, `threads = 5`, `score_interactions = True`. What does the interaction score add beyond the mean expression value?
- [ ] The **statistical method** (`cpdb_statistical_analysis_method.call()`) uses `iterations = 1000` label-shuffling permutations. The p-value is the fraction of shuffled means ≥ the observed mean. What biological assumption does this null model make?
- [ ] Explore the three result tables: `means`, `significant_means` (p < 0.05), and `interaction_scores`. How many interactions survive the significance filter? Which cell-type pairs have the most significant interactions?
- [ ] Use `kpy.plot_cpdb_heatmap()` for a global view of interaction counts per cell-type pair. Which pairs are most communicative?
- [ ] The dot plot (`kpy.plot_cpdb()`) is filtered to VEGF family genes for endothelial and luminal cell-type pairs. Do you see evidence of VEGF signaling between expected partners?
- [ ] Note the key limitation: CellPhoneDB does not know which cells are physically adjacent. High-scoring pairs here may never co-localize in the tissue. Keep this in mind when comparing to notebook 4.

### Reflect

> - What is the practical difference between the basic and statistical CellPhoneDB methods? When would you use each?
> - The `threshold = 0.1` parameter excludes genes expressed in < 10% of a cell type's cells. How would changing this affect the number of testable interactions?
> - If an interaction scores highly here but is absent from notebook 4, what are the possible biological and technical explanations?
> - CellPhoneDB treats all cells of a given type as one pool regardless of tissue location. What biology might this miss in a heterogeneous tumor?

### 💡 Inspect Instead

```python
import pandas as pd

results_dir = "CosMx/results/CellPhoneDB/"
sig_means = pd.read_csv(f"{results_dir}significant_means.txt", sep="\t")
scores    = pd.read_csv(f"{results_dir}interaction_scores.txt", sep="\t")

print(f"Significant interactions: {sig_means.shape[0]}")
scores.sort_values("interaction_score", ascending=False).head(10)
```

---

## Notebook 4 — Spatially-Aware Cell–Cell Communication

**Visualization:** `CosMx/codes/4_cosmx_spatial_cell_com.ipynb` · `seurat_spatial` · ~50 GB  
**Compute:** `CosMx/codes/4_cosmx_spatial_cell_com.R` · `seurat_spatial` · 180 GB, interactive node

### Context and Workflow

**SpatialCellChat** infers cell–cell communication constrained by physical proximity. Because the compute steps are memory-intensive and slow, the workflow is intentionally split:

| Task | Where to run |
|---|---|
| Build Seurat object, create SpatialCellChat object, identify spatially over-expressed genes and interactions, compute communication probabilities, filter, save `.rds` | **`.R` script** on an HPC interactive node |
| Load the saved `.rds` and generate all visualizations | **Jupyter notebook** with the `seurat_spatial` R kernel |

The pre-computed SpatialCellChat object is already at `CosMx/results/SpatialCellChat/SpatialCellChat.rds`.

### Running the Compute Script (optional)

```bash
# Interactive node — no job time limit to worry about
srun --time=04:00:00 --mem=180G -p tier2q --pty bash -l

conda activate /gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2/conda_envs/seurat_spatial
Rscript CosMx/codes/4_cosmx_spatial_cell_com.R
```

The script prints timestamps at each major step so you can track progress.

### Tasks — Understanding the Compute Script

Work through the logic of `4_cosmx_spatial_cell_com.R` even if you do not rerun it:

- [ ] The script loads MTX files from `CosMx/results/Subset_MM/` and builds a Seurat object. A spatial filter (`CenterX_global_px > 35000 & CenterY_global_px < 98000`) keeps a focused sub-region. Why is an additional spatial filter applied here on top of the subset from notebook 1b?
- [ ] Understand the coordinate conversion: CosMx pixel coordinates are converted to micrometers using `conversion.factor = 0.12028` (120.28 nm/pixel). The tolerance factor `tol` is set to half the minimum centroid-to-centroid distance. According to the SpatialCellChat FAQ, why is an exact tolerance factor less critical here than it might seem?
- [ ] The `CellChatDB.human` database is filtered to **Cell-Cell Contact** interactions only before assigning to `chat@DB`. What would you change to also capture secreted or ECM-receptor signaling?
- [ ] `identifyOverExpressedGenes()` uses `selection.method = "meringue"` to identify spatially variable genes. How does this differ from a standard Wilcoxon DE test, and why is spatial variability a more appropriate criterion here?
- [ ] `computeCommunProb()` is called with `distance.use = TRUE`, `scale.distance = 1`, `interaction.range = 250` µm, and `contact.range = 10` µm. What is the biological rationale for these two range parameters?
- [ ] Two filtering steps follow: `filterProbability()` removes non-significant links, and `filterCommunication(min.links = 10, min.cells.sr = 10)` removes interactions supported by too few cells. Why are both filters needed?
- [ ] `computeAvgCommunProb()` uses `nboot = 100` permutations to compute cell-type-level average communication. The significant interactions table is written to `CosMx/results/SpatialCellChat/cell_cell_communication.csv`.

### Tasks — Visualization in the Jupyter Notebook

- [ ] Load the saved SpatialCellChat object. Check `dim(chat@net$prob)` (cell types × cell types × interactions) and `dim(chat@data.signaling)` to understand the object structure.
- [ ] Run `aggregateNet(chat)` to summarize the network and `computeCommunProbPathway(chat)` to group interactions by signaling pathway. What pathways are in `chat@netP$pathways`?
- [ ] Use `spatialDimPlot()` to visualize cell-type assignments on tissue coordinates. Does the spatial layout match the annotation from notebook 2?
- [ ] Use `spatialFeaturePlot()` for a pathway of interest (e.g., `pathway.show = "NOTCH"`). Does the signal localize to a specific tissue region or cell-type boundary?
- [ ] Load `cell_cell_communication.csv` and compare significant interactions to the CellPhoneDB results from notebook 3. Which appear in both? Which are unique to SpatialCellChat? What might a discrepancy indicate?

### Reflect

> - SpatialCellChat requires physical proximity between sender and receiver. What signaling categories would it systematically miss compared to CellPhoneDB?
> - The analysis uses only the Cell-Cell Contact database subset. How would results change if you included Secreted Signaling?
> - `interaction.range = 250` µm approximates paracrine signaling. How would you determine an appropriate range for a different tissue or hypothesis?
> - An interaction is significant in CellPhoneDB (notebook 3) but absent from SpatialCellChat. List the most plausible biological and technical explanations.

### 💡 Inspect Instead

```r
library(SpatialCellChat); library(dplyr)

chat <- readRDS("CosMx/results/SpatialCellChat/SpatialCellChat.rds")
chat <- aggregateNet(chat)

comm <- read.csv("CosMx/results/SpatialCellChat/cell_cell_communication.csv")
comm %>% arrange(desc(prob)) %>% head(20)

spatialDimPlot(chat, group.by = "final_annotation", point.size = 1)
```

---

## Connecting the Analyses

Once you have worked through all five notebooks, consider these integrative questions:

> - CellPhoneDB (notebook 3) and SpatialCellChat (notebook 4) both predict interactions between the same annotated cell types. What does it mean biologically when an interaction is found by both? When found by only one?
> - The annotation in notebook 2 is based on expression clusters. Does the spatial organization of those clusters in notebook 4's `spatialDimPlot` add confidence in the annotation, or raise any concerns?
> - The `Subset_dataset` notebook selects a spatial rectangle from the slide. How would you assess whether the interactions and cell-type compositions you find are representative of the whole tissue?
> - If you were to design a follow-up experiment based on the most compelling interaction found in notebook 4, what validation assay would you choose?
