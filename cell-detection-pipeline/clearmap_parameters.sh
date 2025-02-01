#!/bin/env bash

# ./clearmap_parameters.sh

module load anacondapy/2020.11
conda activate ClearMap_spock
xvfb-run -d python /jukebox/YOUR_DIR/lightsheet-pipeline/clearmap/clearmap_parameters.py
