#!/bin/bash

# Set directories
INPUT_DIR="/home/user/SR/BCSIR-dengu_2025/data"
OUTPUT_DIR="/home/user/SR/BCSIR-dengu_2025/MixCR/results"
THREADS=8

# Create output directory
mkdir -p $OUTPUT_DIR

# Get all R1 fastq files
for R1_FILE in $INPUT_DIR/*_R1.fastq.gz; do
    # Get base filename (without _R1.fastq.gz)
    BASE_NAME=$(basename "$R1_FILE" _R1.fastq.gz)
    
    # Construct R2 filename
    R2_FILE="$INPUT_DIR/${BASE_NAME}_R2.fastq.gz"
    
    # Check if R2 file exists
    if [[ ! -f "$R2_FILE" ]]; then
        echo "Warning: R2 file not found for $BASE_NAME: $R2_FILE"
        continue
    fi
    
    # Output prefix for this sample
    SAMPLE_OUTPUT="$OUTPUT_DIR/$BASE_NAME"
    
    echo "Processing sample: $BASE_NAME"
    echo "R1: $R1_FILE"
    echo "R2: $R2_FILE"
    echo "Output: $SAMPLE_OUTPUT"
    echo "----------------------------------------"
    
    # Run MixCR with correct syntax for version 4.6.0
    mixcr analyze rna-seq \
        --species hsa \
        --threads $THREADS \
        --force-overwrite \
        "$R1_FILE" "$R2_FILE" \
        "$SAMPLE_OUTPUT"
    
    echo "Completed analysis for: $BASE_NAME"
    echo "========================================"
done

echo "All samples processed!"
