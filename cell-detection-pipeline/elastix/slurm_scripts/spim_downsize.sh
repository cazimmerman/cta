#!/bin/env bash

#SBATCH --job-name=downsize
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=20G
#SBATCH --contiguous
#SBATCH --time=5:00:00
#SBATCH --output=elastix/logs/%u_%x_%A_%a.out
#SBATCH --mail-user=YOUR_EMAIL@princeton.edu
#SBATCH --mail-type=END

echo "directory: `pwd`"
echo "user:      `whoami`"
echo "host:      `hostname`"
cat /proc/$$/status | grep Cpus_allowed_list
echo ""

module load anacondapy/2020.11
conda activate lightsheet

xvfb-run -d python /jukebox/YOUR_DIR/lightsheet-pipeline/elastix/spim_downsize.py ${sample_dir} ${imaging_request} ${output_rootpath} ${atlas}
