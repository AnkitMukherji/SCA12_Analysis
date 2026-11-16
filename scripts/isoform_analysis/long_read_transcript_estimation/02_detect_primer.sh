#!/usr/bin/env bash
set -euo pipefail

INPUT="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_long_read/fastqs"
OUT="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_long_read/primer_qc"

mkdir -p "$OUT/logs"

UPM="CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT"
UPM_RC="ACTCTGCGTTGATACCACTGCTTGCCCTATAGTGAGTCGTATTAG"

GSP="GATTACGCCAAGCTT"
GSP_RC="AAGCTTGGCGTAATC"

for fq in "$INPUT"/*.fastq.gz; do

    sample=$(basename "$fq" .fastq.gz)

    echo "Processing: $sample"

    # Full-length UPM at 5' end
    cutadapt \
        -g "^${UPM}" \
        --action=none \
        -e 0.10 \
        -O 45 \
        -o /dev/null \
        "$fq" \
        > "$OUT/logs/${sample}_UPM_forward.txt"

    # Reverse-complement UPM at 3' end
    cutadapt \
        -a "${UPM_RC}" \
        --action=none \
        -e 0.10 \
        -O 45 \
        -o /dev/null \
        "$fq" \
        > "$OUT/logs/${sample}_UPM_reverse.txt"

    # GSP cloning overhang at 5' end
    cutadapt \
        -g "^${GSP}" \
        --action=none \
        -e 0.10 \
        -O 15 \
        -o /dev/null \
        "$fq" \
        > "$OUT/logs/${sample}_GSP_forward.txt"

    # Reverse-complement GSP overhang at 3' end
    cutadapt \
        -a "${GSP_RC}" \
        --action=none \
        -e 0.10 \
        -O 15 \
        -o /dev/null \
        "$fq" \
        > "$OUT/logs/${sample}_GSP_reverse.txt"

done
