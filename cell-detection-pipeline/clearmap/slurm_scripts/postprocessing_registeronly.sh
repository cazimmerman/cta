#!/bin/env bash

#SBATCH --job-name=postprocessingregister
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=10G
#SBATCH --contiguous
#SBATCH --time=8:00:00
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

xvfb-run -d python /jukebox/YOUR_DIR/lightsheet-pipeline/clearmap/clearmap_postprocessing_register.py ${sample_dir} ${imaging_request} ${output_rootpath} ${atlas}
