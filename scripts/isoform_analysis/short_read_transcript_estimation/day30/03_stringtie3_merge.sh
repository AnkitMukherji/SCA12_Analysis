#!/bin/bash

# Configuration Paths
DIR="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read/day30/assemblies"
GENCODE="/mnt/faruq2/lab_data/reference/gencode/gencode50/gencode.v50.annotation.gtf"
REFSEQ="/mnt/faruq2/lab_data/reference/ncbi_refseq_gtf/GCF_000001405.40_GRCh38.p14_genomic_ucsc_nomenclature_empty_transcripts_removed.gtf"
THREADS=20

# Create lists of GTF files to merge
find "$DIR" -name "transcripts_gencode.gtf" -type f > "$DIR/gencode_assembly_list.txt"
find "$DIR" -name "transcripts_refseq.gtf" -type f > "$DIR/refseq_assembly_list.txt"

# Merge Gencode transcripts using the input list file (-A or xargs)
stringtie --merge \
  -G "$GENCODE" \
  -p "$THREADS" \
  -o "$DIR/gencode_transcripts_merged.gtf" \
  $(cat "$DIR/gencode_assembly_list.txt")

# Merge RefSeq transcripts using the input list file
stringtie --merge \
  -G "$REFSEQ" \
  -p "$THREADS" \
  -o "$DIR/refseq_transcripts_merged.gtf" \
  $(cat "$DIR/refseq_assembly_list.txt")
