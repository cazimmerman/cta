#!/bin/env bash

# ./clearmap_pipeline.sh $request_name $imaging_request $sample_name
# ./clearmap_pipeline.sh zimmerman_01 1 111

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
	# Creates a JSON file (the blockfile) which just contains the number of blocks to do cell detection
	echo "Step 0: Initializing ClearMap workspace:"
	module load anacondapy/2020.11
	conda activate ClearMap_spock
	python /jukebox/YOUR_DIR/lightsheet-pipeline/clearmap/clearmap_preprocessing.py ${sample_dir} ${imaging_request} ${output_rootpath} 2>&1 | tee
	conda deactivate
	module unload anacondapy/2020.11
	sample_name=$(basename ${sample_dir})
	blockfile="${output_rootpath}/$1/${sample_name}/${imaging_request}/rawdata/resolution_3.6x/block_processing_info.json"
	echo $blockfile

	# Read JSON file to figure out how many blocks there are to process
	if [ -a $blockfile ]
	then
		# echo $blockfile
		n_blocks=`cat ${blockfile} | cut -d " " -f 2 | cut -d "}" -f 1`
		echo "Have ${n_blocks} blocks for clearmap to process"
	else
		echo "Block file not found. Skipping this sample"
		echo ""
		continue
	fi
	max_array_job=`echo "(${n_blocks}+${blocks_per_job}-1)/${blocks_per_job}-1" | bc` # bash equivalent of ceil()
	# echo "Max array job for cell detection: ${max_array_job}"
	echo "Submitting batch jobs:"

	## Create stitched.npy memmap volume file
	OUT1=$(sbatch --parsable --export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},output_rootpath=${output_rootpath} \
	clearmap/slurm_scripts/create_memmap_vol.sh)
	echo "Step 1: Memmmap volume step"
	echo $OUT1

	# # ## Run the individual blocks in array jobs
	OUT2=$(sbatch --parsable --dependency=afterok:${OUT1##* } \
	--exclude=./bad_nodenames.txt \
	--export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},blocks_per_job=${blocks_per_job},\
output_rootpath=${output_rootpath},clearmap_params_file=${clearmap_params_file} \
	--array=0-${max_array_job} clearmap/slurm_scripts/cell_detect.sh)
	echo "Step 2: Cell detection on blocks:"
	echo $OUT2

	# ## Merge the blocks
	OUT3=$(sbatch --parsable --dependency=afterok:${OUT1##* }:${OUT2##* } \
	--export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},output_rootpath=${output_rootpath},\
clearmap_params_file=${clearmap_params_file} \
	clearmap/slurm_scripts/merge_blocks.sh)
	echo "Step 3: Merge cell detection blocks:"
	echo $OUT3

	## Register cells to atlas space and make CSV data frame of counts in each region
	# Dependent on merge block step and inverse registration step
	OUT4=$(sbatch --parsable --dependency=afterok:${OUT1##* }:${OUT3##* } \
	--export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},output_rootpath=${output_rootpath},\
atlas=${atlas} clearmap/slurm_scripts/postprocessing_registeronly.sh)
 echo "Step 4: Register cells to atlas"
 echo $OUT4





done
