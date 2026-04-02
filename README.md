# Spatial transcriptomics - Bruker: Workshop material

This folder includes datasets, scripts and conda environments you can use to recapitulate the various analysis presented during the workshop. 

## Workshop description

This workshop introduces the fundamentals of Bruker Spatial Transcriptomics data analysis with a focus on quality control (QC), cell deconvolution, and cell–cell communication. Participants will learn the core principles of the platform, key experimental design considerations, and essential QC strategies to ensure reliable spatial transcriptomics data. The session will also highlight analytical approaches for resolving cell-type composition within spatial spots and exploring cellular interactions within tissue microenvironments, illustrated through real-world case examples.

**Topics covered**
 - Introduction to Bruker Spatial Transcriptomics technology
 - Study design considerations for spatial transcriptomics experiments
 - Quality control strategies and common pitfalls
 - Cell-type deconvolution in GeoMx
 - Cell–cell communication analysis in spatial contexts - CosMx


## Conda environments

A set of 3 conda environments are avaiable for your use in Randi. 

To use them, load conda and use the command: 
`conda activate /gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2/conda_env/{env_name}`

    * **scverse_spatial**: (/gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2/conda_env/scverse_spatial) Python-based, includes all the scverse packages (scanpy, squidpy, etc) and cellphoneDB.
    * **seurat_spatial**: R-based, 
    * geomx_env: R-based, GeomxTools. 

## Jupyter notebooks

For your convenience, the codes sections include jupyter notebooks, I also included a script to launch jupyter with a given conda environment. In all the notebooks, the job was launched from the main workshop folder, and all the paths inside them are relative to this folder too. 

`sbatch --mem=50GB ./launch_jupyter_conda.sh /gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2/conda_env/{env_name}`


## 1. GeoMx

    * Dataset pre-processing and QC. (1_cosmx_preprocessing_qc.ipynb)
    * DE analysis with limma-voom.
    * Cell-type deconvolution.


You can also refer to the codes 

## 2. CosMx

    * Dataset pre-processing and QC. 
    * Annotation
    * Cell-cell communication.


### Relevant links: 

[CosMx Best practices](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/cosmx-smi-best-practices/)
