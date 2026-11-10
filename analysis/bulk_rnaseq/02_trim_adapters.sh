#!/bin/bash

STORAGE_DIR="/mnt/vstor/SOM_PATH_DKB50/members/rxr456"
TRIM_DIR="$STORAGE_DIR"/trim
BAM_DIR="$STORAGE_DIR"/mapped
RAW_DIR="$STORAGE_DIR"/raw_processed
SORT_BAM_DIR="$STORAGE_DIR/sortindSAM"
QUANT_DIR="$STORAGE_DIR/quant"
ADAPTER_DIR="/home/rxr456/adapters.fasta"

mkdir -p $TRIM_DIR
echo "Start trimming the adapter..."
for SAMPLE_DIR in "$RAW_DIR"/D*; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIR")
    R1=$(ls "$SAMPLE_DIR"/*_1.fq.gz)
    R2=$(ls "$SAMPLE_DIR"/*_2.fq.gz)

    O1="$TRIM_DIR/$SAMPLE_NAME""_trimmed_1.fastq"
    O2="$TRIM_DIR/$SAMPLE_NAME""_trimmed_2.fastq"

    echo "Running Cutadapt on $SAMPLE_NAME..."
    echo "$O1"
    echo "$O2"

    cutadapt -a file:"$ADAPTER_DIR" -A file:"$ADAPTER_DIR" -o "$O1" -p "$O2" "$R1" "$R2" -j 36 2> $TRIM_DIR/${SAMPLE_NAME}_cutadapt.txt
    echo "Cutadapt on $SAMPLE_NAME done"
done
echo "Cutadapt quantification complete for all samples!"