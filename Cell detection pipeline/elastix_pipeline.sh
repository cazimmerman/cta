#!/bin/env bash

# ./elastix_pipeline.sh $request_name $imaging_request $sample_name
# ./elastix_pipeline.sh zimmerman_01 1 111

username=cz15
output_rootpath="/jukebox/YOUR_DIR/clearmap-output"
atlas='Allen' # 'Princeton' or 'Allen'

# Makes sure that you supplied at least two command line arguments
if [ "$#" -lt 2 ]
then
	echo "Error: Incorrect number of command line arguments"
	echo ""
	echo "Usage: clearmap_pipeline.sh request_name imaging_request"
	echo "e.g.: clearmap_pipeline.sh zimmerman_01 imaging_request_1"
	exit 1
fi

request_name=$1
imaging_request="imaging_request_$2"
request_dir="/jukebox/YOUR_DATA/$1"

# Check to see if sample_name argument was provided
if [ "$#" -eq 3 ]
then
	sample_name=$3
	echo "Running script for Request name: $1, Imaging request: $2, Sample name: $3"
	declare -a arr=("/jukebox/YOUR_DATA/${request_name}/${request_name}-${sample_name}")
else
	echo "Running script for Request name: $1, Imaging request: $2, all samples"
	declare -a arr=(${request_dir}/*)
fi
echo ""
# echo $arr
echo "Sample directories that will be run are:"
for d in "${arr[@]}"
do
	echo $d
done
echo ""
# Loop through all sample dirs in the request folder and run the full pipeline for each

blocks_per_job=1 # number of cell blocks to process in a single array job. Currently only works with 1 right now
# due to slurm error. TODO - allow parallel processing in each array job

# for sample_dir in "${request_dir}"/*
for sample_dir in "${arr[@]}"
do

	echo "Working on sample: $sample_dir"
	echo ""
	# Step 0: Link over raw 488 and 642 files to destination directory.
	echo "Step 0: Generating symlink files:"
	module load anacondapy/2020.11
	conda activate ClearMap_spock
	python /jukebox/YOUR_DIR/lightsheet-pipeline/elastix/spim_symlink.py ${sample_dir} ${imaging_request} ${output_rootpath} 2>&1 | tee
	conda deactivate
	module unload anacondapy/2020.11
	echo ""

	## Downsizing, both channels one per array job, can start without dependency
	OUT0=$(sbatch --parsable --export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},output_rootpath=${output_rootpath},\
atlas=${atlas} --array=0-1 elastix/slurm_scripts/spim_downsize.sh)
	echo "Step 1: Downsizing:"
	echo $OUT0

	# Inverse registration, both transformations one per array job
	OUT1=$(sbatch --parsable --dependency=afterok:${OUT0##* } --array=0-1 \
	--export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},output_rootpath=${output_rootpath},\
atlas=${atlas} elastix/slurm_scripts/spim_inverse_register.sh)


	echo "Step 2: Inverse registration:"
	echo $OUT1

done
