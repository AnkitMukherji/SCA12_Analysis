#!/usr/bin/env bash
set -euo pipefail

BASE="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_long_read"

INPUT="${BASE}/fastqs"
OUTPUT="${BASE}/primer_qc/trimmed_fastqs"
LOGS="${BASE}/primer_qc/trimming_logs"

mkdir -p "$OUTPUT" "$LOGS"

UPM="CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT"

OVERHANG="GATTACGCCAAGCTT"
GSP_GENE="GGCGGTTTTCCATGTTGAGATCCCAGACTTTG"

revcomp() {
    printf '%s' "$1" | rev | tr 'ACGT' 'TGCA'
}

UPM_RC=$(revcomp "$UPM")
OVERHANG_RC=$(revcomp "$OVERHANG")

shopt -s nullglob

for fq in "$INPUT"/*.fastq.gz; do

    sample=$(basename "$fq" .fastq.gz)

    echo "Processing: $sample"

    # remove UPM in either orientation
    cutadapt \
        -g "^${UPM}" \
        -a "${UPM_RC}$" \
        -e 0.10 \
        -O 40 \
        -j 8 \
        --json "$LOGS/${sample}_UPM.json" \
        -o "$OUTPUT/${sample}.UPM_trimmed.fastq.gz" \
        "$fq" \
        > "$LOGS/${sample}_UPM.log"

    # remove ONLY the 15-nt cloning overhang
    cutadapt \
        -g "^${OVERHANG}" \
        -a "${OVERHANG_RC}$" \
        -e 0.10 \
        -O 15 \
        -j 8 \
        --json "$LOGS/${sample}_overhang.json" \
        -o "$OUTPUT/${sample}.primer_trimmed.fastq.gz" \
        "$OUTPUT/${sample}.UPM_trimmed.fastq.gz" \
        > "$LOGS/${sample}_overhang.log"

    gzip -t "$OUTPUT/${sample}.primer_trimmed.fastq.gz"

    echo "Completed: $sample"

done
