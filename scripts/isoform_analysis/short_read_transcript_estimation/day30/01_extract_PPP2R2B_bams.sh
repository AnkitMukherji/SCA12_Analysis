#!/bin/bash

set -euo pipefail

# Configuration Paths
DIR="/mnt/faruq2/lab_users/shahrumi/RNA-Seq/SCA12_DIV30_neurons_rna_seq/RNA_Fastq_Dragen_Mapping"
OUT="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read/day30/ppp2r2b_bams"
THREADS=20
PPP2R2B_LOC="chr5:146400000-147200000"

# Sample Arrays
CASES=(
    "102A_S720_L007.bam" "102B_S721_L007.bam" "104A_S761_L008.bam" "104B_S762_L008.bam"
    "105A_S763_L008.bam" "105C_S764_L008.bam" "107A_S765_L008.bam" "107B_S766_L008.bam"
    "109A_S767_L008.bam" "109B_S768_L008.bam" "112A_S769_L008.bam" "112B_S770_L008.bam"
    "113A_S771_L008.bam" "113B_S772_L008.bam" "114A_S773_L008.bam" "114B_S774_L008.bam"
    "115A_S775_L008.bam" "115B_S776_L008.bam" "116A_S777_L008.bam" "116B_S778_L008.bam"
    "117A_S779_L008.bam" "117B_S780_L008.bam"
)
CONTROLS=(
    "24A_S713_L007.bam" "24B_S714_L007.bam" "25A_S590_L006.bam" "25B_S715_L007.bam"
    "28A_S716_L007.bam" "28B_S717_L007.bam" "29A_S718_L007.bam" "29B_S719_L007.bam"
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
