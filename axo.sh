#!/bin/bash
# usage: bash run_axo_directly.sh <ctsfile> <output_dir>

ctsfile=$1
output_dir=$2

mkdir -p $output_dir

source /mnt/disk7t/xwj/soft/miniforge3/etc/profile.d/conda.sh
conda activate py3.12

python parse_axo_directly.py $ctsfile $output_dir