#!/bin/env bash

#SBATCH --job-name=pystripe
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=10G
#SBATCH --contiguous
#SBATCH --time=10:00:00
#SBATCH --output=pystripe/logs/%u_%x_%A.out
#SBATCH --mail-user=YOUR_EMAIL@princeton.edu
#SBATCH --mail-type=END

echo "directory: `pwd`"
echo "user:      `whoami`"
echo "host:      `hostname`"
cat /proc/$$/status | grep Cpus_allowed_list
echo ""

module load anacondapy/2020.11
conda activate lightsheet

echo "Input directory (path to stitched images):" $input_dir
echo "Path to flat.tiff file:" $flat
echo "Output directory (does not need to exist):" $corrected_dir

pystripe -i $input_dir -f $flat -o $corrected_dir -s1 256 -s2 512 -w 'db3'
