#!/usr/bin/env bash
set -euo pipefail

REF="/mnt/faruq2/lab_data/reference/human_ref/hg38/dragen_fasta/hg38.fa"
IN="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_long_read"

mkdir -p "$IN/alignment/logs" "$IN/alignment/qc"

shopt -s nullglob
fastqs=("$IN"/primer_qc/trimmed_fastqs/*.primer_trimmed.fastq.gz)

if (( ${#fastqs[@]} == 0 )); then
    echo "ERROR: No trimmed FASTQ files found" >&2
    exit 1
fi

for fq in "${fastqs[@]}"; do

    sample=$(basename "$fq" .primer_trimmed.fastq.gz)

    echo "Aligning $sample"

    minimap2 \
        -ax splice \
        -s 40 \
        -G 500k \
        --MD \
        --secondary=no \
        -t 24 \
        "$REF" \
        "$fq" \
        2> "$IN/alignment/logs/${sample}.minimap2.log" \
    | samtools sort \
        -@ 8 \
        -o "$IN/alignment/${sample}.bam" \
        -

    samtools index "$IN/alignment/${sample}.bam"

    samtools flagstat \
        "$IN/alignment/${sample}.bam" \
        > "$IN/alignment/qc/${sample}.flagstat.txt"

    samtools quickcheck -v \
        "$IN/alignment/${sample}.bam"

    echo "Completed $sample"

done
