# Spatial transcriptomics - Bruker: Workshop material
### Created by Diana Vera Cruz

This repository contains jupyter notebooks and other code developed for the **Spatial Transcriptomics Workshop** focused on Bruker technologies. It is intended to be used as reference. 

---

## Workshop description

This workshop introduces the fundamentals of Bruker Spatial Transcriptomics data analysis with a focus on quality control (QC), cell deconvolution, and cell–cell communication. Participants will learn the core principles of the platform, key experimental design considerations, and essential QC strategies to ensure reliable spatial transcriptomics data. The session will also highlight analytical approaches for resolving cell-type composition within spatial spots and exploring cellular interactions within tissue microenvironments, illustrated through real-world case examples.

**Topics covered**
 - Introduction to Bruker Spatial Transcriptomics technology
 - Study design considerations for spatial transcriptomics experiments
 - Quality control strategies and common pitfalls
 - Cell-type deconvolution in GeoMx
 - Cell–cell communication analysis in spatial contexts - CosMx

 ---

## Workshop resources


**Recording** [Recording](https://cri.uchicago.edu/wp-content/uploads/2026/04/DVCruz-Rcording-of-Spatial-Transcriptomics-Seminar_April2026.mp4)

 **Workshop slides** [Box link](https://uchicago.box.com/s/6kjmu2s7x3hng24zrgseehgrr2arnkz1)

 **Raw data** [Box link](https://uchicago.box.com/s/tp0j9x6bt8t6iqtpbqqbu6muaa4u9v4f)
  Also within workshop folder in Randi.

 **Workshop folder in Randi** *_`/gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2`_*

To access content from this and other CRI seminars visit: 
[CRI Seminar series](https://cri.uchicago.edu/seminar-series/)

---


## Repository navigation

The repository is separated by technology for clarity. Independent jupyter notebooks related to each of the topics covered in the workshops are provided for you to navigate and inspect, they include the step by step process without needing to run anything. The workshop slides contain the name of the notebook you might want to check out so you can follow along in the same order as we covered then.

### Technologies

- [`GeoMx`](GeoMx/) — GeoMx QC, DE, and SpatialDecon
- [`CosMx`](CosMx/) — CosMx QC, annotation, CellPhoneDB, SpatialCellChat

---

## Quick Reference

Links to the main topics covered and the notebooks and other content related to them. A detailed description and additional material used per topic is included per tecnology, check their READMEs.

### GeoMx

For this example we use a dataset from FFPE and FF kidney tissue sections from diabetic and healthy individuals. 

| Topic | Notebook / Resource |
|---|---|
| GeoMx obj from DCC, PKC data, QC and filtering with GeomxTools | [1_geomx_setup_qc.ipynb](GeoMx/codes/1_geomx_setup_qc.ipynb) |
| Normalization, PCA and DE analysis with limma-voom | [2_geomx_limma_norm_DE.ipynb](GeoMx/codes/2_geomx_limma_norm_DE.ipynb) |
| Cell-type deconvolution — SpatialDecon + Kidney HCA reference | [3_geomx_SpatialDecon.ipynb](GeoMx/codes/3_geomx_SpatialDecon.ipynb) |


### CosMx

For this section we use a dataset from FFPE human breast tissue, tested with the Whole-transcriptome CosMx panel. 

| Topic | Notebook / Resource |
|---|---|
| AnnData obj from flat files, QC and filtering | [1_cosmx_preprocessing_qc.ipynb](CosMx/codes/1_cosmx_preprocessing_qc.ipynb) |
| Normalization, dimension reduction, clustering and annotation  | [2_cosmx_cell_annotation.ipynb](CosMx/codes/2_cosmx_cell_annotation.ipynb) |
| Cell-cell communication with CellPhoneDB v5 | [3_cosmx_cell_comunication.ipynb](CosMx/codes/3_cosmx_cell_comunication.ipynb) |
| Spatially-aware Cell-cell communication with SpatialCellChat | [4_cosmx_spatial_cell_com.ipynb](CosMx/codes/4_cosmx_spatial_cell_com.ipynb) |

---


## Accessing the workshop content in Randi

The workshop folder, including datasets, codes and results, along with the conda environments used are located at: 

`/gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2`


### Workshop folder structure

```
workshop/
├── data/                            
│   ├── GeoMx_Human_kidney/          # [See Datasets below]
│   ├── CosMx_Human_breast/          # [See Datasets below]
│   ├── cellphone_db_v5/             # CellPhoneDB v5.0.0 database zip (pre-downloaded)
│   └── celltypist_models/           # CellTypist models (pre-downloaded)
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
│   └── README.md              # GeoMx README, self-study description
├── CosMx/
│   ├── codes/                       # Python + R kernel notebooks
│   ├── results/                     # Pre-computed outputs
│   └── README.md              # CosMx README, self-study description
├── launch_jupyter_conda.sh
└── README.md
```

### Conda Environments

Three environments are available on Randi for your use. 

```bash
conda activate /gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2/conda_envs/{env_name}
```

| Environment | Language | Used for |
|---|---|---|
| `geomx_env` | R | GeoMx notebooks 1–3 |
| `scverse_spatial` | Python | CosMx notebooks 1–3 |
| `seurat_spatial` | R | CosMx notebook 4 (SpatialCellChat) |


If you want to recreate the environments in a different location (outside of Randi), you can find in this repository the YAML files, they also contain notes on their installation.

[conda environments YAMLs](conda_env/)


### For experienced users only

There is also an option for running the notebook directly in Randi, You need to use VSCode and log in Randi as remote, as well as have the jupyter extension. You can easily use the conda environments directly and copy the notebooks or create new ones. If you encounter any issue join us for office hours or raise an issue in the repository for some help. 

We provide a script to submit a job launching jupyter with the desired conda environment, the following example is a command ran within the workshop folder using the geomx_env conda enviroment.

```bash
# Submit a Jupyter job 
sbatch --mem=50GB ./launch_jupyter_conda.sh \
  /gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2/conda_envs/geomx_env
  ## Check err file, it will contain the URL to set the kernel. 
```

The script *_launch_jupyter_conda.sh_* will load the selected conda environment and set the libraries, you can set up the memory and partition too. After submitting the job, you will observe two files, an err and out files named jupyter_nb_*.err/out.  
Wait a few seconds and check the content of the  jupyter_nb_*.err file, it will contain the URL that you can use to set the kernel.  

** DO NOT FORGET TO CANCEL THE JOB ONCE YOU FINISHED**

`scancel {JOB_ID}`


Suggested memory per section:

| Section | Environment | Memory |
|---|---|---|
| GeoMx notebooks 1–3 | `geomx_env` | 10 GB |
| CosMx notebook 1 + Subset | `scverse_spatial` | 100–180 GB |
| CosMx notebooks 2–3 | `scverse_spatial` | 50 GB |
| CosMx notebook 4 (visualization) | `seurat_spatial` | 30 GB |
| CosMx notebook 4 (compute, `.R` script) | `seurat_spatial` | 180 GB, use run_4_cosmx.sh script |

---

