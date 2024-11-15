#!/bin/bash

set -euo pipefail

# This script intends to simulate multiple samples.
# It copies the barcode01/reads.fq within barcode01 to create 15 new barcodes.

# ARGS:
if [[ $# != 2 ]]; then
    >&2 echo "ERROR: two arguments required: output path and input file."
    exit 1
fi

OUTPUT_PATH=$1  # Output folder. This folder is created from scratch.
FASTQ_TO_PROPAGATE=$2  # FASTQ file to be used as reference

# SET PATHS
mkdir $OUTPUT_PATH/test_data
PATH_SAMPLES=$OUTPUT_PATH/test_data

# Recreate many barcodes 1x barcode01 to show the order of samples in report
for i in {01..15}
do
    mkdir $PATH_SAMPLES/barcode$i
    echo barcode$i
    cp $FASTQ_TO_PROPAGATE $PATH_SAMPLES/barcode$i/
done


