#!/bin/bash
#SBATCH --job-name=snakeHead
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000
#SBATCH --time=48:00:00
#SBATCH --output=snakeHead.out
#SBATCH --partition=all


module purge
#echo $PATH
source /home/sukem128/miniconda3/etc/profile.d/conda.sh
conda info
conda --version
conda activate snakemake
snakemake --profile CAUCluster --keep-going --local-cores 1
