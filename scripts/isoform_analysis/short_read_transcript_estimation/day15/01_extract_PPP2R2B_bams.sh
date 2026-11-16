#!/bin/bash

set -euo pipefail

# Configuration Paths
DIR="/mnt/faruq2/lab_users/vivek/analysis/RNA-Seq/SCA12/batch2/dragen"
OUT="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read/ppp2r2b_bams"
THREADS=20
PPP2R2B_LOC="chr5:146400000-147200000"

# Sample Arrays
CASES=(
    "1108_NEU1.bam" "1108_NEU2.bam" "1108_NEU3.bam" "1108_NSCP10.bam"
    "ATC2555_NEU3.bam" "ATC2555_NEU5.bam" "ATC2555_NEU6.bam" 
    "ATC2555_NSCP14.bam" "ATC2555_NSCP8.bam" "150_N2.bam" 
    "150_NEU_5JAN1.bam" "150_NEU_5JAN2.bam" "150_NSC_P14.bam" "150_NSC_P8.bam"
)
CONTROLS=(
    "GOC29_N3.bam" "GOC29_N4.bam"
)

# Merge arrays
ALL_SAMPLES=("${CASES[@]}" "${CONTROLS[@]}")

# Ensure output directory exists
mkdir -p "$OUT"

# Process samples
for sample in "${ALL_SAMPLES[@]}"; do
    input_bam="$DIR/$sample"
    output_bam="$OUT/$sample"

    echo "Processing: $sample..."

    # Extract specific region using multi-threading
    if ! samtools view -@ "$THREADS" -b "$input_bam" "$PPP2R2B_LOC" > "$output_bam"; then
        echo "Error: Failed to subset $sample" >&2
        continue
    fi

    # Index the subsetted BAM file using multi-threading
    samtools index -@ "$THREADS" "$output_bam"
done

echo "All samples processed successfully!"
