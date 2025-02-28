#!/bin/env bash

# ./cellfinder_pipeline.sh $request_name $imaging_request $sample_name
# ./cellfinder_pipeline.sh zimmerman_01 1 111

request_name=$1
imaging_request=$2
sample_name=$3

sbatch cellfinder/cellfinder_inference_inparts.sh $request_name $imaging_request $sample_name 1
sbatch cellfinder/cellfinder_inference_inparts.sh $request_name $imaging_request $sample_name 2
sbatch cellfinder/cellfinder_inference_inparts.sh $request_name $imaging_request $sample_name 3
sbatch cellfinder/cellfinder_inference_inparts.sh $request_name $imaging_request $sample_name 4
