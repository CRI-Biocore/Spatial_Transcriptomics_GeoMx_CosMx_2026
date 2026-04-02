#!/bin/bash
#SBATCH --job-name=jupyter_nb
#SBATCH --output=jupyter_nb-%j.out
#SBATCH --error=jupyter_nb-%j.err
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50GB

## USAGE: sbatch launch_jupyter.sh /path/to/your/conda_env
## sbatch --mem=10GB launch_jupyter_conda.sh /path/to/your/env
module load miniconda3/23.1.0

## Conda environment ready for python. 
## If no argument is provided, defaults to sc_v5.
CONDA_ENV=${1:-/gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2/conda_envs/scverse_spatial}
source activate $CONDA_ENV

print "Using conda environment: $CONDA_ENV"

## Make sure the $PATH is set. 
export PATH=$CONDA_ENV/bin:$PATH

echo "PATH: $PATH"

cd $SLURM_SUBMIT_DIR

HOST_IP=`/sbin/ip route get 8.8.8.8 | awk '{print $7;exit}'`
PORT_NUM=$(shuf -i15001-30000 -n1)


echo "Jupyter URL: http://$HOST_IP:$PORT_NUM"

jupyter lab  --no-browser --ip=$HOST_IP --port=$PORT_NUM

