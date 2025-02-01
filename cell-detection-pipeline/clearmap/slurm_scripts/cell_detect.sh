#!/bin/env bash

#SBATCH --job-name=celldetect
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=15G
#SBATCH --contiguous
#SBATCH --time=4:00:00
#SBATCH --output=clearmap/logs/%u_%x_%A_%a.out
#SBATCH --mail-user=YOUR_EMAIL@princeton.edu
#SBATCH --mail-type=END

echo "directory: `pwd`"
echo "user:      `whoami`"
echo "host:      `hostname`"
cat /proc/$$/status | grep Cpus_allowed_list
echo ""

authdir=./authfiles
export XDG_RUNTIME_DIR=$authdir

module load anacondapy/2020.11
conda activate ClearMap_spock

tsleep=`echo "$SLURM_ARRAY_TASK_ID*0.2+1" | bc`
echo "Sleeping for $tsleep seconds"
sleep $tsleep

request_and_sample=`echo $sample_dir | cut -d "/" -f6,7`
cell_block_fname=${output_rootpath}/${request_and_sample}/${imaging_request}/rawdata/resolution_3.6x/cells_blocks/cells_block${SLURM_ARRAY_TASK_ID}.p
EXIT=0
while [[ ! -f "$cell_block_fname" ]]
do
	echo "In while loop, running cell detection on block ${SLURM_ARRAY_TASK_ID}"
	python /jukebox/YOUR_DIR/lightsheet-pipeline/clearmap/clearmap_cell_detect.py ${sample_dir} ${imaging_request} ${blocks_per_job} ${output_rootpath} ${clearmap_params_file}
	EXIT=$?
done
echo "Exit status is $EXIT"
exit 0
