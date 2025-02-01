#!/bin/env bash

#SBATCH --job-name=FosKDE
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=15G
#SBATCH --contiguous
#SBATCH --time=5:00:00
#SBATCH --output=logs/%u_%x_%A_%a.out
#SBATCH --mail-user=cz15@princeton.edu
#SBATCH --mail-type=END

module load anacondapy/2020.11
conda activate neuro

python make_heatmap_fullres_weighted_modkernel.py step2 ${condition_name} ${chunk_size}
