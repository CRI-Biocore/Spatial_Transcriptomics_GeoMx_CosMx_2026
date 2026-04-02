# Spatial transcriptomics - Bruker: Workshop material
### Created by Diana Vera Cruz

This repository contains notebooks, conda environments, and pre-computed results for the **Spatial Transcriptomics Workshop** focused on Bruker technologies. The workshop covers quality control, normalization, differential expression, cell-type deconvolution, and cell–cell communication analysis across two complementary spatial platforms.

Pre-computed results are provided for every step, so each notebook can be explored independently without re-running prior ones.

---

## Workshop description

This workshop introduces the fundamentals of Bruker Spatial Transcriptomics data analysis with a focus on quality control (QC), cell deconvolution, and cell–cell communication. Participants will learn the core principles of the platform, key experimental design considerations, and essential QC strategies to ensure reliable spatial transcriptomics data. The session will also highlight analytical approaches for resolving cell-type composition within spatial spots and exploring cellular interactions within tissue microenvironments, illustrated through real-world case examples.

**Topics covered**
 - Introduction to Bruker Spatial Transcriptomics technology
 - Study design considerations for spatial transcriptomics experiments
 - Quality control strategies and common pitfalls
 - Cell-type deconvolution in GeoMx
 - Cell–cell communication analysis in spatial contexts - CosMx

 **WORKSHOP SLIDES** (link to slides)

---

## Quick Reference

A summary of all available notebooks and resources. Click any notebook link to jump directly to that section, or open the platform-specific study guide for guided exercises.

### GeoMx

ADD QUICK NOTES ON DATASET HERE

| Topic | Notebook / Resource |
|---|---|
| Data loading — DCC, PKC, annotation Excel | `GeoMx/codes/1_geomx_setup_qc.ipynb` |
| Segment QC — reads, saturation, nuclei, area, NTC | `GeoMx/codes/1_geomx_setup_qc.ipynb` |
| Probe QC — low ratio and Grubbs outlier filtering | `GeoMx/codes/1_geomx_setup_qc.ipynb` |
| Limit of Quantification (LOQ) and gene filtering | `GeoMx/codes/1_geomx_setup_qc.ipynb` |
| Q3 normalization and PCA | `GeoMx/codes/2_geomx_limma_norm_DE.ipynb` |
| Differential expression — limma-voom + duplicateCorrelation | `GeoMx/codes/2_geomx_limma_norm_DE.ipynb` |
| Volcano plot and DEG heatmaps | `GeoMx/codes/2_geomx_limma_norm_DE.ipynb` |
| Background estimation from negative probes | `GeoMx/codes/3_geomx_SpatialDecon.ipynb` |
| Cell-type deconvolution — SpatialDecon + Kidney HCA reference | `GeoMx/codes/3_geomx_SpatialDecon.ipynb` |
| Cell abundance visualization — barplot and heatmap | `GeoMx/codes/3_geomx_SpatialDecon.ipynb` |

### CosMx

ADD QUICK NOTES ON DATASET HERE

| Topic | Notebook / Resource |
|---|---|
| Flat file ingestion — memory-efficient AnnData construction | `CosMx/codes/1_cosmx_preprocessing_qc.ipynb` |
| QC metrics — negative probes, counts, area, FOV-level checks | `CosMx/codes/1_cosmx_preprocessing_qc.ipynb` |
| Cell filtering and spatial visualization post-QC | `CosMx/codes/1_cosmx_preprocessing_qc.ipynb` |
| Spatial subsetting — FOV-based region selection | `CosMx/codes/Subset_dataset.ipynb` |
| MTX export for Seurat / SpatialCellChat | `CosMx/codes/Subset_dataset.ipynb` |
| Normalization, HVG selection, PCA, UMAP | `CosMx/codes/2_cosmx_cell_annotation.ipynb` |
| Leiden clustering at multiple resolutions | `CosMx/codes/2_cosmx_cell_annotation.ipynb` |
| CellTypist annotation — multi-model consensus | `CosMx/codes/2_cosmx_cell_annotation.ipynb` |
| Spatial scatter plots of cell types | `CosMx/codes/2_cosmx_cell_annotation.ipynb` |
| CellPhoneDB v5 — basic and statistical methods | `CosMx/codes/3_cosmx_cell_comunication.ipynb` |
| Ligand–receptor interaction scoring and visualization | `CosMx/codes/3_cosmx_cell_comunication.ipynb` |
| SpatialCellChat — compute (heavy, HPC script) | `CosMx/codes/4_cosmx_spatial_cell_com.R` |
| SpatialCellChat — object creation and visualization | `CosMx/codes/4_cosmx_spatial_cell_com.ipynb` |
| Spatially-aware communication — pathway and network plots | `CosMx/codes/4_cosmx_spatial_cell_com.ipynb` |

---

## Self-Study Guides

Step-by-step post-workshop guides with tasks, reflection questions, and inspect-instead code snippets for each notebook:

- [`GeoMx/README.md`](GeoMx/README.md) — GeoMx QC, DE, and SpatialDecon
- [`CosMx/README.md`](CosMx/README.md) — CosMx QC, annotation, CellPhoneDB, SpatialCellChat

---

## Conda Environments

Three environments are available on Randi. Activate with:

```bash
conda activate /gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2/conda_envs/{env_name}
```

| Environment | Language | Used for |
|---|---|---|
| `geomx_env` | R | GeoMx notebooks 1–3 |
| `scverse_spatial` | Python | CosMx notebooks 1–3 |
| `seurat_spatial` | R | CosMx notebook 4 (SpatialCellChat) |

To recreate any environment, you can find here the YAML files, they also contain notes on their installation.

```bash
conda env create -f conda_envs/{env_name}.yaml
```

Additional R packages for the GeoMx environment can be installed with:
```bash
Rscript conda_envs/geomx_Rpackages_install.R
```

---

## Launching Jupyter on the HPC (Randi)

All notebooks use paths **relative to the workshop root folder**. Always launch Jupyter from that directory.

```bash
# Submit a Jupyter job — replace {env_name} with the appropriate environment
sbatch --mem=50GB ./launch_jupyter_conda.sh \
  /gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2/conda_envs/{env_name}
```

Suggested memory per section:

| Section | Environment | Memory |
|---|---|---|
| GeoMx notebooks 1–3 | `geomx_env` | 10–50 GB |
| CosMx notebook 1 + Subset | `scverse_spatial` | 100–180 GB |
| CosMx notebooks 2–3 | `scverse_spatial` | 50–100 GB |
| CosMx notebook 4 (visualization) | `seurat_spatial` | 50 GB |
| CosMx notebook 4 (compute, `.R` script) | `seurat_spatial` | 180 GB, interactive node |

---

## Repository Structure

```
workshop/
├── data/                            # RANDI ONLY.
│   ├── GeoMx_Human_kidney/          # [See Datasets below]
│   ├── CosMx_Human_breast/          # [See Datasets below]
│   ├── cellphone_db_v5/             # CellPhoneDB v5.0.0 database zip (pre-downloaded)
│   └── celltypist_models/           # CellTypist models (pre-downloaded)
│       ├── Cells_Adult_Breast.pkl
│       ├── Immune_All_Low.pkl
│       └── Immune_All_High.pkl
├── conda_envs/
│   ├── geomx_env/                   # R: GeomxTools, limma, edgeR, SpatialDecon
│   ├── geomx_env.yaml
│   ├── geomx_Rpackages_install.R
│   ├── scverse_spatial/             # Python: scanpy, squidpy, celltypist, cellphonedb
│   ├── scverse_spatial_env.yaml
│   ├── seurat_spatial/              # R: Seurat, SpatialCellChat
│   └── seurat_spatial_env.yaml
├── GeoMx/
│   ├── codes/                       # R kernel notebooks
│   ├── results/                     # Pre-computed outputs
│   └── STUDY_GUIDE.md              # GeoMx self-study guide
├── CosMx/
│   ├── codes/                       # Python + R kernel notebooks
│   ├── results/                     # Pre-computed outputs
│   └── STUDY_GUIDE.md              # CosMx self-study guide
├── launch_jupyter_conda.sh
└── README.md
```

---

## Relevant Links

- [GeoMx Bioconductor Workflow](https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html)
- [SpatialDecon Bioconductor Vignette](https://bioconductor.org/packages/release/bioc/vignettes/SpatialDecon/inst/doc/SpatialDecon_vignette_NSCLC.html)
- [CosMx Best Practices — Bruker](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/cosmx-smi-best-practices/)
- [SpatialCellChat Tutorial](https://htmlpreview.github.io/?https://github.com/jinworks/SpatialCellChat/blob/master/tutorial/SpatialCellChat_analysis_of_spatial_transcriptomics_data.html)
- [CellChat FAQ — Spatial Data](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/FAQ_on_applying_CellChat_to_spatial_transcriptomics_data.html)
- [CellPhoneDB Documentation](https://github.com/ventolab/CellphoneDB)
- [CellTypist](https://www.celltypist.org/)