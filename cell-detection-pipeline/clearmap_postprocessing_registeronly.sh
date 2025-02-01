#!/bin/env bash

# ./clearmap_postprocessing_registeronly.sh $request_name $imaging_request $sample_name
# ./clearmap_postprocessing_registeronly.sh zimmerman_01 1 111

username=cz15
output_rootpath="/jukebox/YOUR_DIR/clearmap-output"
clearmap_params_file='/jukebox/YOUR_DIR/lightsheet-pipeline/clearmap/cell_detection_parameter.p'
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

	## Register cells to atlas space and make CSV data frame of counts in each region
	# Dependent on merge block step and inverse registration step
	OUT4=$(sbatch --parsable --export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},output_rootpath=${output_rootpath},\
atlas=${atlas} clearmap/slurm_scripts/postprocessing_registeronly.sh)
	echo "Step 4: Register cells to atlas"
	echo $OUT4

done
