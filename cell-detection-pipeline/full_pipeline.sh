#!/bin/env bash

# ./pystripe_pipeline.sh $request_name $imaging_request $sample_name
# ./pystripe_pipeline.sh zimmerman_01 1 111

request_name=$1
imaging_request="imaging_request_$2"
sample_name=$3
fpath="${request_name}/${request_name}-${sample_name}/${imaging_request}"
module load anacondapy/2020.11
jobOne = $(`python /jukebox/YOUR_DIR/lightsheet-pipeline/pystripe/pystripe_batch.py ${fpath}`)
echo ${jobOne}