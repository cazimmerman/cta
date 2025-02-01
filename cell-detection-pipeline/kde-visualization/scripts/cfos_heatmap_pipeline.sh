#!/bin/env bash

condition_name=csds_susceptible_c0
chunk_size=50000 # how many kernel evals will happen per core
# n_chunks=
# max_array_job=`echo "${n_chunks}-1" | bc`
# OUT0=$(sbatch --parsable \
# 	--export=ALL,condition_name=${condition_name},chunk_size=${chunk_size} \
# 	slurm_scripts/cfos_heatmap_step0.sh)
# echo $OUT0

# OUT1=$(sbatch --parsable --dependency=afterok:${OUT0} \
# 	--export=ALL,condition_name=${condition_name},chunk_size=${chunk_size} \
# 	slurm_scripts/cfos_heatmap_step1.sh)
# echo $OUT1

OUT2=$(sbatch --parsable --array=0-480 \
	--export=ALL,condition_name=${condition_name},chunk_size=${chunk_size} \
	--exclude=./bad_nodenames.txt slurm_scripts/cfos_heatmap_step2.sh)
echo $OUT2

# OUT3=$(sbatch --parsable   \
# 	--export=ALL,condition_name=${condition_name},chunk_size=${chunk_size} \
# 	slurm_scripts/cfos_heatmap_step3.sh)
# echo $OUT3

# OUT4=$(sbatch --parsable   \
# 	--export=ALL,condition_name=${condition_name},chunk_size=${chunk_size} \
# 	slurm_scripts/cfos_heatmap_step4.sh)
# echo $OUT4
