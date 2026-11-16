#!/usr/bin/env bash
set -euo pipefail

DIR="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_long_read"
GENOME="/mnt/faruq2/lab_data/reference/human_ref/hg38/dragen_fasta/hg38.fa"
GTF="$DIR/reference/gencode.v50.annotation.primer_capped.gtf"
REFSEQ="$DIR/reference/GCF_000001405.40_GRCh38.p14_genomic_ucsc_nomenclature_empty_transcripts_removed.primer_capped.gtf"
BAM="$DIR/flair/pooled/pooled_reads.bam"
JUNCTIONS="$DIR/long_read_junctions/pooled_reads.junctions.bed"
OUTDIR="$DIR/flair/primary"

mkdir -p "$OUTDIR"

# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------

for file in "$GENOME" "$GTF" "$REFSEQ" "$BAM" "$JUNCTIONS"; do
    if [[ ! -s "$file" ]]; then
        echo "ERROR: Missing or empty file: $file" >&2
        exit 1
    fi
done

samtools quickcheck -v "$BAM"

# ------------------------------------------------------------
# Validate junction BED
# ------------------------------------------------------------

awk '
BEGIN { FS="\t" }
NF < 6 {
    print "ERROR: Junction BED has fewer than 6 columns at line " NR > "/dev/stderr"
    exit 1
}
$5 !~ /^[0-9]+$/ {
    print "ERROR: Non-numeric junction support at line " NR > "/dev/stderr"
    exit 1
}
$6 != "+" && $6 != "-" {
    print "ERROR: Invalid strand at line " NR > "/dev/stderr"
    exit 1
}
' "$JUNCTIONS"

# ------------------------------------------------------------
# Run against both annotations
# ------------------------------------------------------------

names=("GENCODE" "RefSeq")
gtfs=("$GTF" "$REFSEQ")

for i in "${!names[@]}"; do

    name="${names[$i]}"
    annotation="${gtfs[$i]}"

    OUT="$OUTDIR/$name"
    mkdir -p "$OUT"

    echo "================================================"
    echo "Running FLAIR v3.0.1: $name"
    echo "================================================"

    flair transcriptome \
        -b "$BAM" \
        -g "$GENOME" \
        -f "$annotation" \
        --junction_bed "$JUNCTIONS" \
        --junction_support 2 \
        --ss_window 10 \
        --support 3 \
        --stringent \
        --check_splice \
        --end_window 25 \
        --max_ends 10 \
        --filter comprehensive \
        --keep_intermediate \
        -t 24 \
        -o "$OUT/pooled_reads_transcriptome" \
        2>&1 | tee "$OUT/flair_transcriptome.log"

    echo "Completed: $name"

done

echo "All FLAIR reconstructions completed."
