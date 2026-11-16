#!/usr/bin/env bash
set -euo pipefail

INPUT="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_long_read/fastqs"
OUT="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_long_read/primer_qc/logs"

for fq in "$INPUT"/*.fastq.gz;
do
    
    sample=$(basename "$fq" .fastq.gz)
    echo "Processing: $sample"

    cutadapt \
        -a "UPM_RC=ACTCTGCGTTGATACCACTGCTTGCCCTATAGTGAGTCGTATTAG" \
        -a "GSP_RC=AAGCTTGGCGTAATC" \
        --action=none \
        -e 0.1 \
        -O 15 \
        --info-file "$OUT/${sample}_primer_matches.tsv" \
        -o /dev/null \
        "$fq"
done
