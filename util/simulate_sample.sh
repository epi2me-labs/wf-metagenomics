#!/bin/bash

set -euo pipefail

# This script intends to simulate a real time run,
# where new fq files are input while the wf is being executed.
# It copies the barcode01/reads.fq within barcode01 to create another 2 fq files within this directory (renaming read names)
# Besides it creates 12 new barcodes with 1 fq and also creates 3 fq in barcode10.

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
# Copy the initial barcode file to be used within the real time simulation
mkdir $PATH_SAMPLES/barcode01
cp $FASTQ_TO_PROPAGATE $PATH_SAMPLES/barcode01/

sleep 15s # to let nextflow start
# simulate two more fq files. Barcode01 should have results x3
for i in {1..2}
do
    sed "s/seq/seq_$i/g" $FASTQ_TO_PROPAGATE > $PATH_SAMPLES/barcode01/reads_$i.fq
    echo reads_$i.fq
    sleep 2s
done
# Recreate many barcodes 1x barcode01 to increase the possibilities of a bad timing
for i in {02..12}
do
    mkdir $PATH_SAMPLES/barcode$i
    echo barcode$i
    cp $FASTQ_TO_PROPAGATE $PATH_SAMPLES/barcode$i/
    sleep 3s
done
# duplicate fq into another barcode more to make sure everything works as expected
for i in {1..2}
do
    sed "s/seq/seq_$i/g" $FASTQ_TO_PROPAGATE > $PATH_SAMPLES/barcode10/reads_$i.fq
    echo reads_$i.fq
    sleep 2s
done

