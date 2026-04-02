#!/bin/bash
#SBATCH --job-name=spatialCellChat
#SBATCH --partition=tier2q
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --mem=200GB

module load miniconda3/23.1.0

CONDA_ENV=/gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2/conda_envs/seurat_spatial
source activate $CONDA_ENV

echo "Using conda environment: $CONDA_ENV"

export PATH=$CONDA_ENV/bin:$PATH

Rscript /gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2/CosMx/codes/4_cosmx_spatial_cell_com.R


