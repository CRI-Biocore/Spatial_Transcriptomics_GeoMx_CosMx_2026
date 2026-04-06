# Spatial transcriptomics - Bruker: Workshop material
### Created by Diana Vera Cruz

This repository contains notebooks, conda environments, and pre-computed results for the **Spatial Transcriptomics Workshop** focused on Bruker technologies. The workshop covers quality control, normalization, differential expression, cell-type deconvolution, and cell–cell communication analysis across two complementary spatial platforms.

---

## Workshop description

This workshop introduces the fundamentals of Bruker Spatial Transcriptomics data analysis with a focus on quality control (QC), cell deconvolution, and cell–cell communication. Participants will learn the core principles of the platform, key experimental design considerations, and essential QC strategies to ensure reliable spatial transcriptomics data. The session will also highlight analytical approaches for resolving cell-type composition within spatial spots and exploring cellular interactions within tissue microenvironments, illustrated through real-world case examples.

**Topics covered**
 - Introduction to Bruker Spatial Transcriptomics technology
 - Study design considerations for spatial transcriptomics experiments
 - Quality control strategies and common pitfalls
 - Cell-type deconvolution in GeoMx
 - Cell–cell communication analysis in spatial contexts - CosMx

 **WORKSHOP SLIDES** [Box link](https://uchicago.box.com/s/6kjmu2s7x3hng24zrgseehgrr2arnkz1)

 **Raw data** [Box link](https://uchicago.box.com/s/tp0j9x6bt8t6iqtpbqqbu6muaa4u9v4f)


---


## Self-Study Guides

Step-by-step post-workshop guides, you can also navigate and check the notebooks independently, to check the step by step process without running anything or you can also run the notebooks in Randi. 

- [`GeoMx`](GeoMx/README.md) — GeoMx QC, DE, and SpatialDecon
- [`CosMx`](CosMx/README.md) — CosMx QC, annotation, CellPhoneDB, SpatialCellChat

---

## Quick Reference

A summary of all available notebooks and resources. Click any notebook link to jump directly to that section, or open the platform-specific study guide for guided exercises.

### GeoMx

For this example we use a dataset from FFPE and FF kidney tissue sections from diabetic and healthy individuals. 


| Topic | Notebook / Resource |
|---|---|
| GeoMx obj from DCC, PKC data, QC and filtering with GeomxTools | `GeoMx/codes/1_geomx_setup_qc.ipynb` |
| Normalization, PCA and DE analysis with limma-voom | `GeoMx/codes/2_geomx_limma_norm_DE.ipynb` |
| Cell-type deconvolution — SpatialDecon + Kidney HCA reference | `GeoMx/codes/3_geomx_SpatialDecon.ipynb` |


### CosMx

For this section we use a dataset from FFPE human breast tissue, tested with the Whole-transcriptome CosMx panel. 


| Topic | Notebook / Resource |
|---|---|
| AnnData obj from flat files, QC and filtering | `CosMx/codes/1_cosmx_preprocessing_qc.ipynb` |
| Obj subsetting an MTX export | `CosMx/codes/Subset_dataset.ipynb` |
| Normalization, dimension reduction, clustering and annotation  | `CosMx/codes/2_cosmx_cell_annotation.ipynb` |
| Cell-cell communication with CellPhoneDB v5 | `CosMx/codes/3_cosmx_cell_comunication.ipynb` |
| Spatially-aware Cell-cell communication with SpatialCellChat — compute (heavy, HPC script) | `CosMx/codes/4_cosmx_spatial_cell_com.R` |
| Spatially-aware Cell-cell communication with SpatialCellChat — notebook with visualizations | `CosMx/codes/4_cosmx_spatial_cell_com.ipynb` |

---


## Running jupyter notebooks on the HPC (Randi)

The workshop folder, including datasets, codes and results, along with the conda environments is located at: 

`/gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2`

All the notebooks can be rerun in Randi using the help script *_run_jupyter_notebook.sh_* provided in the main project folder. This script loads the correct libraries, runs the script and store the output and a copy of the notebook in a folder of your choice.

`sbatch --mem={10G} ./run_jupyter_notebook.sh {conda_env_name} {notebook_workshop_relative_path}`

Please create a folder (`out_dir`) for your usage where you would like to store the notebook and the output files. After running the script, it will create a notebook in the output directory and also a `results` subfolder, if not present already, with all output files.

**Example, running 1_geomx_setup_qc.ipynb**
```bash
# Submit jupyter notebook job. 
## Set output directory of your choice. 
out_dir=$HOME/workshop 
mkdir -p $out_dir 

sbatch --mem=10G /gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2/run_jupyter_notebook.sh \
  geomx_env GeoMx/codes/1_geomx_setup_qc.ipynb $out_dir

```

If there is any error, you will notice that the job is no longer running, and if a notebook was saved in the output directory, it will have a red header with a note on which chunk caused the error.


**For experienced users only**

There is also an option for running the notebook directly in Randi, You need to use VSCode and log in Randi as remote, as well as have the jupyter extension. You can easily use the conda environments directly and copy the notebooks or create new ones. If you encounter any issue join us for office hours or raise an issue in the repository for some help. 

```bash
# Submit a Jupyter job — replace {env_name} with the appropriate environment
sbatch --mem=50GB ./launch_jupyter_conda.sh \
  /gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2/conda_envs/{env_name}
  ## Check err file, it will contain the URL to set the kernel. 
  ## Within 
```

The script *_launch_jupyter_conda.sh_* will load the selected conda environment and set the libraries, you can set up the memory and partition too. After submitting the job, you will observe two files, an err and out files named jupyter_nb_*.err/out.  
Wait a few seconds and check the content of the  jupyter_nb_*.err file, it will contain the URL that you can use to set the kernel.  


Suggested memory per section:

| Section | Environment | Memory |
|---|---|---|
| GeoMx notebooks 1–3 | `geomx_env` | 10 GB |
| CosMx notebook 1 + Subset | `scverse_spatial` | 100–180 GB |
| CosMx notebooks 2–3 | `scverse_spatial` | 50 GB |
| CosMx notebook 4 (visualization) | `seurat_spatial` | 30 GB |

| CosMx notebook 4 (compute, `.R` script) | `seurat_spatial` | 180 GB, use run_4_cosmx.sh script |

---


## Conda Environments


Three environments are available on Randi. 
The two helper scripts already load and add the correct path per conda environment when you ran a notebook, but you can also use them for your personal scripts using: 


```bash
conda activate /gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2/conda_envs/{env_name}
```

| Environment | Language | Used for |
|---|---|---|
| `geomx_env` | R | GeoMx notebooks 1–3 |
| `scverse_spatial` | Python | CosMx notebooks 1–3 |
| `seurat_spatial` | R | CosMx notebook 4 (SpatialCellChat) |


If you want to recreate the environments, you can find in this repository the YAML files, they also contain notes on their installation.

[conda environments YAMLs](conda_env/)


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