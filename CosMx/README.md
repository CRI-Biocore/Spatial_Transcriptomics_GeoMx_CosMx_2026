# CosMx — Post-Class Self-Study Guide

> Work through this guide at your own pace after the workshop. Each notebook is independent — pre-computed results are in `CosMx/results/`, so you can open any notebook directly without having run the ones before it. The numbers suggest a logical progression, not a strict dependency.

---

## Dataset

For this section we use a dataset from FFPE human breast tissue, tested with the Whole-transcriptome CosMx panel as well as a 64 proteomics panel. The dataset was downloaded from [Bruker, multiomic data](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/cosmx-human-multiomic-breast-dataset/)
For this workshop, we focus only on the transcriptomics panel results. 


---

## Analysis workflow

Use the link to explore each step. The notebooks will show code, plots and notes to describe the process.

### Topics covered in the workshop

| Topic | Notebook / Resource |
|---|---|
| AnnData obj from flat files, QC and filtering | [1_cosmx_preprocessing_qc.ipynb](codes/1_cosmx_preprocessing_qc.ipynb) |
| Normalization, dimension reduction, clustering and annotation  | [2_cosmx_cell_annotation.ipynb](codes/2_cosmx_cell_annotation.ipynb) |
| Spatially-aware Cell-cell communication with SpatialCellChat — Visualizations | [4_cosmx_spatial_cell_com.ipynb](codes/4_cosmx_spatial_cell_com.ipynb) |


### Additional resources

| Topic | Notebook / Resource |
|---|---|
| Obj subsetting an MTX export | [Subset_dataset.ipynb](codes/Subset_dataset.ipynb) |
| Cell-cell communication with CellPhoneDB v5 | [3_cosmx_cell_comunication.ipynb](codes/3_cosmx_cell_comunication.ipynb) |
| Spatially-aware Cell-cell communication with SpatialCellChat — computation (heavy, HPC script) | [4_cosmx_spatial_cell_com.R](codes/4_cosmx_spatial_cell_com.R) |

---

## Notes

All the content for this section was created using the  `scverse_spatial` (python-based) and the `seurat_spatial` (R-based) conda environments. All the analysis require larger memory, try with at least **100GB**. If using Randi, don't forget to use the `tier2q` partition. 


---

## Step 1 — Preprocessing & QC

[1_cosmx_preprocessing_qc.ipynb](codes/1_cosmx_preprocessing_qc.ipynb)

### Context

CosMx flat files are loaded manually: the expression matrix (`*exprMat_file.csv.gz`) is read with `pyarrow` and immediately converted to a sparse `csr_matrix` to halve peak RAM usage, then combined with cell metadata (`*metadata_file.csv.gz`) into an `AnnData` object. Global pixel coordinates (`CenterX_global_px`, `CenterY_global_px`) are stored in `obsm["spatial"]` for `squidpy` compatibility. QC uses negative probe percentage as the primary noise metric alongside transcript counts, gene counts, and cell area.

### Tasks

- [ ] Inspect the multiple coordinates: What is the difference between local and global coordinates.
- [ ] Identify the negative probe columns (names contain `"NegPrb"`). `sc.pp.calculate_qc_metrics()` is called with `qc_vars=["negative_probe"]`, adding `pct_counts_negative_probe` to `adata.obs`. Why is this metric preferred over mitochondrial gene percentage for CosMx data?
- [ ] Review the FOV-level QC plots (average negative probe signal and average counts per FOV). What is the purpose of this step before applying cell-level filters?
  Which individual filter removes the most cells? What does `MAX_COUNTS = 7500` guard against?
- [ ] Compare the spatial scatter plots of `total_counts` and `pct_counts_negative_probe` before and after filtering. Is the noise spatially random, or concentrated in particular regions or FOVs?

### Reflect

> - Why is `pct_counts_negative_probe` a better noise metric here than mitochondrial gene percentage?
> - What biological or tissue-specific considerations should guide `MIN_AREA` and `MAX_AREA`? How would these change for a tissue with very large cells (e.g., neurons)?


---

## Step 2 — Clustering & Cell Type Annotation


[2_cosmx_cell_annotation.ipynb](codes/2_cosmx_cell_annotation.ipynb)

### Context

Starting from `BreastCancer_subset.h5ad` (31,999 cells), this notebook follows a standard scanpy workflow: normalization to the median library size, HVG selection (top 3,000 genes, `seurat_v3`), PCA, neighbor graph (15 neighbors, 17 PCs), UMAP, and Leiden clustering at resolutions 0.2, 0.4, and 0.6.

Annotation uses **CellTypist** with three pre-downloaded models from `data/celltypist_models/`:
- `Cells_Adult_Breast.pkl` — breast-specific cell types
- `Immune_All_Low.pkl` — broad immune categories
- `Immune_All_High.pkl` — fine-grained immune types

All three are run with `majority_voting=True` over `leiden_0.4` clusters. A consensus annotation selects the best model per cluster based on mean prediction confidence.

### Tasks

- [ ] Examine the elbow plot from `sc.pl.pca_variance_ratio()`. The notebook uses `n_pcs = 17` — where is the elbow, and does the choice seem well-supported?
- [ ] Look at the UMAP colored by `Mean.PanCK`, `Mean.CD68`, `Mean.CD45`, and `Mean.DAPI` (morphological markers from segmentation). Which populations do these highlight before any clustering label is applied?
- [ ] Compare spatial scatter plots for Leiden resolutions 0.2, 0.4, and 0.6. Which resolution produces the most spatially coherent clusters? The notebook selects `leiden_0.4` — do you agree?
- [ ] Review the `rank_genes_groups` dotplot for `leiden_0.4`. Do the top 5 markers per cluster give an intuition for cell identity before running CellTypist?
- [ ] Inspect `consensus_df` and `conf_summary`: most-frequent label and mean confidence per cluster per model. The final annotation uses `Immune_All_High` labels for most clusters but switches to `Cells_Adult_Breast` for clusters 0, 2, 3, and 7 (luminal epithelial). Why are two models needed?
- [ ] Review the luminal marker dotplot (`ESR1`, `PGR`, `FOXA1`, `KRT8`, `MKI67`, `TOP2A`, `KRT5`, `KRT14`, etc.) used to characterize epithelial clusters before finalizing annotation.
- [ ] Verify the `final_annotation` spatial scatter plot. Do luminal, immune, and stromal populations occupy distinct tissue regions?

### Reflect

> - The notebook uses three CellTypist models and builds a consensus. What are the risks of using only one model for a complex tissue like breast cancer?
> - The luminal marker dotplot suggests clusters 0, 2, 3, and 7 could be further subdivided. What biological distinction would you look for to justify splitting them?
> - Does the spatial organization of `final_annotation` support the expression-derived clustering, or do you see unexpected mismatches?


---

## Step 3 — Cell–Cell Communication (Expression-Based)


[3_cosmx_cell_comunication.ipynb](codes/3_cosmx_cell_comunication.ipynb)

### Context

This notebook uses **CellPhoneDB v5** (database at `data/cellphone_db_v5/cellphonedb.zip`, pre-downloaded) to infer ligand–receptor interactions from expression alone — without spatial constraints. CellPhoneDB requires two inputs: a lognorm expression `h5ad` (with the `lognorm` layer promoted to `X`) and a metadata TSV with columns `barcode_sample` and `cell_type`. Both the basic analysis method and the statistical permutation-based method are run. Results are visualized with `ktplotspy`.

This analysis intentionally ignores spatial coordinates — it establishes a baseline of putative interactions to compare against the proximity-constrained SpatialCellChat results in notebook 4.

### Tasks

- [ ] Explore the three result tables: `means`, `significant_means` (p < 0.05), and `interaction_scores`. How many interactions survive the significance filter? Which cell-type pairs have the most significant interactions?
- [ ] The dot plot (`kpy.plot_cpdb()`) is filtered to VEGF family genes for endothelial and luminal cell-type pairs. Do you see evidence of VEGF signaling between expected partners?
- [ ] Note the key limitation: CellPhoneDB does not know which cells are physically adjacent. High-scoring pairs here may never co-localize in the tissue. Keep this in mind when comparing to notebook 4.

### Reflect

> - The `threshold = 0.1` parameter excludes genes expressed in < 10% of a cell type's cells. How would changing this affect the number of estable interactions?
> - If an interaction scores highly here but is absent from notebook 4, what are the possible biological and technical explanations?
> - CellPhoneDB treats all cells of a given type as one pool regardless of tissue location. What biology might this miss in a heterogeneous tumor?


---

## Step 4 — Spatially-Aware Cell–Cell Communication

**Visualization:** [4_cosmx_spatial_cell_com.ipynb](codes/4_cosmx_spatial_cell_com.ipynb) 

**Compute:** [4_cosmx_spatial_cell_com.R](codes/4_cosmx_spatial_cell_com.R)

### Context and Workflow

**SpatialCellChat** infers cell–cell communication constrained by physical proximity. Because the compute steps are memory-intensive and slow, the workflow is intentionally split:


### Tasks

Work through the logic of `4_cosmx_spatial_cell_com.R` even if you do not rerun it:


- [ ] Understand the coordinate conversion: CosMx pixel coordinates are converted to micrometers using `conversion.factor = 0.12028` (120.28 nm/pixel). The tolerance factor `tol` is set to half the minimum centroid-to-centroid distance. According to the SpatialCellChat FAQ, why is an exact tolerance factor less critical here than it might seem?
- [ ] Compare significant interactions to the CellPhoneDB results from step 3. How similar does the amount of LR interactions change between cell types pairs.

### Reflect

> - SpatialCellChat requires physical proximity between sender and receiver. What signaling categories would it systematically miss compared to CellPhoneDB?
> - `interaction.range = 250` µm approximates paracrine signaling. How would you determine an appropriate range for a different tissue or hypothesis?
> - An interaction is significant in CellPhoneDB (notebook 3) but absent from SpatialCellChat. List the most plausible biological and technical explanations.


---

## Connecting the Analyses

Once you have worked through all the content, consider these integrative questions:

> - CellPhoneDB (notebook 3) and SpatialCellChat (notebook 4) both predict interactions between the same annotated cell types. What does it mean biologically when an interaction is found by both? When found by only one?
> - The annotation in notebook 2 is based on expression clusters. Does the spatial organization of those clusters in notebook 4's `spatialDimPlot` add confidence in the annotation, or raise any concerns?
> - If you were to design a follow-up experiment based on the most compelling interaction found in notebook 4, what validation assay would you choose?

----

## Relevant Links

- [CosMx Best Practices — Bruker](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/cosmx-smi-best-practices/)
- [SpatialCellChat Tutorial](https://htmlpreview.github.io/?https://github.com/jinworks/SpatialCellChat/blob/master/tutorial/SpatialCellChat_analysis_of_spatial_transcriptomics_data.html)
- [CellChat FAQ — Spatial Data](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/FAQ_on_applying_CellChat_to_spatial_transcriptomics_data.html)
- [CellPhoneDB Documentation](https://github.com/ventolab/CellphoneDB)
- [CellTypist](https://www.celltypist.org/)


