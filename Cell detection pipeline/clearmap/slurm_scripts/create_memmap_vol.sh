#!/bin/env bash

#SBATCH --job-name=memmap
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=10G
#SBATCH --contiguous
#SBATCH --time=5:00:00
#SBATCH --output=clearmap/logs/%u_%x_%A.out
#SBATCH --mail-user=YOUR_EMAIL@princeton.edu
#SBATCH --mail-type=END

echo "directory: `pwd`"
echo "user:      `whoami`"
echo "host:      `hostname`"
cat /proc/$$/status | grep Cpus_allowed_list
echo ""

module load anacondapy/2020.11
conda activate ClearMap_spock

request_and_sample=`echo $sample_dir | cut -d "/" -f6,7`
stitched_file_fname=${output_rootpath}/${request_and_sample}/${imaging_request}/rawdata/resolution_3.6x/stitched.npy
python /jukebox/YOUR_DIR/lightsheet-pipeline/clearmap/clearmap_create_memmap_vol.py ${sample_dir} ${imaging_request} ${output_rootpath}
if [[ ! -f "$stitched_file_fname" ]]
then
	echo "Stitched file not found after code ran. Error."
	exit 1
fi
exit 0
