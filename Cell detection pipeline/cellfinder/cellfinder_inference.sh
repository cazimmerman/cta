#!/bin/env bash

#SBATCH --job-name=cellfinder_inference
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=250G
#SBATCH --time=4-00:00:00
#SBATCH --output=logs/%u_%x_%j.out
#SBATCH --mail-user=YOUR_EMAIL@princeton.edu
#SBATCH --mail-type=END

echo "directory: `pwd`"
echo "user:      `whoami`"
echo "host:      `hostname`"
cat /proc/$$/status | grep Cpus_allowed_list
echo ""

request_name=$1
imaging_request="imaging_request_$2"
sample_name=$3
signal_dir="/jukebox/YOUR_DATA/$request_name/$request_name-$sample_name/$imaging_request/rawdata/resolution_3.6x/Ex_647_Em_690_corrected/"
background_dir="/jukebox/YOUR_DATA/$request_name/$request_name-$sample_name/$imaging_request/rawdata/resolution_3.6x/Ex_488_Em_525_corrected/"
output_dir="/jukebox/YOUR_DIR/clearmap-output/cellfinder/$request_name/$request_name-$sample_name/"
cell_candidates="/jukebox/YOUR_DIR/clearmap-output/$request_name/$request_name-$sample_name/$imaging_request/rawdata/resolution_3.6x/cells_raw_filtered.xml"
model="/jukebox/YOUR_DIR/lightsheet-pipeline/cellfinder/model.h5"
echo "Request name:     $request_name"
echo "Sample name:      $sample_name"
echo "Imaging session:  $imaging_request"
echo "Output directory: $output_dir"
echo "Cellfinder model: $model"
echo ""

if [ -f "$cell_candidates" ]; then
    mkdir -p "$output_dir/points/"
    cp -n $cell_candidates "$output_dir/points/cells.xml"
    module load anacondapy/2020.11
    conda activate cellfinder
    cellfinder -s $signal_dir -b $background_dir -o $output_dir -v 2 2 2 --orientation "asl" --trained-model $model --batch-size 128 --no-register --no-detection --no-analyse --no-figures
else
    echo "Cell candidates file does not exist!"
    echo "$cell_candidates"
fi
