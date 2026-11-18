#!/usr/bin/env bash
set -euo pipefail

TOOL="/mnt/faruq2/lab_data/tools/SQANTI3"
DIR="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_long_read/sqanti3"

# Reference
GENOME="$DIR/input_data/references/hg38.chr5.fa"

# FLAIR reconstructed isoforms
ISO="$DIR/input_data/flair_isoforms/flair_isoforms.chr5.gtf"

# Reference GTF
REF="$DIR/input_data/references/refseq.chr5.gtf"

# Short-read junctions
SJ_DIR="$DIR/input_data/short_read_junctions"

# Short-read BAMs
SR_BAM="$DIR/input_data/bams"

# RefTSS annotations from CAGE
CAGE="$DIR/input_data/orthogonal_data/CAGE.chr5.bed"

OUT="$DIR/refseq"

mkdir -p "$OUT"

"$TOOL/sqanti3_qc.py" \
    --isoforms "$ISO" \
    --refGTF "$REF" \
    --refFasta "$GENOME" \
    --coverage "$SJ_DIR" \
    --SR_bam "$SR_BAM" \
    --ratio_TSS_metric max \
    --isoform_hits \
    --CAGE_peak "$CAGE" \
    --report both \
    --cpus 24 \
    --chunks 4 \
    --output PPP2R2B_SQANTI3 \
    --dir "$OUT"
