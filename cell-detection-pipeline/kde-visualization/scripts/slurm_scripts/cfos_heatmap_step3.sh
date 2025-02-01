#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 1                      # number of cores
#SBATCH -t 10                # time (minutes)
#SBATCH -o logs/cfos_heatmap_step3_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/cfos_heatmap_step3_%j.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 100000 #100gbs 

module load anacondapy/2020.11
conda activate neuro

python make_heatmap_fullres.py step3 ${condition1_name} ${chunk_size} 
