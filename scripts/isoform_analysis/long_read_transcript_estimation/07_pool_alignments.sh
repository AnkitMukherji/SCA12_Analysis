#!/usr/bin/env bash
set -euo pipefail

DIR="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_long_read"

mkdir -p "$DIR/flair/pooled"

shopt -s nullglob
bams=("$DIR"/alignment/*.bam)

if (( ${#bams[@]} != 6 )); then
    echo "ERROR: Expected 6 BAMs, found ${#bams[@]}" >&2
    exit 1
fi

samtools merge \
    -@ 16 \
    -f \
    "$DIR/flair/pooled/pooled_reads.bam" \
    "${bams[@]}"

samtools index \
    -@ 16 \
    "$DIR/flair/pooled/pooled_reads.bam"

samtools quickcheck -v \
    "$DIR/flair/pooled/pooled_reads.bam"

echo "Pooled BAM created successfully."
