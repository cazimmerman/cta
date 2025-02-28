#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 8                      # number of cores
#SBATCH -t 60                # time (minutes)
#SBATCH -o logs/cfos_heatmap_step4_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/cfos_heatmap_step4_%j.err        # STDERR #add _%a to see each array job

module load anacondapy/2020.11
conda activate neuro

python make_heatmap_fullres_weighted_modkernel.py step4 ${condition_name} ${chunk_size} 
