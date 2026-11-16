#!/bin/bash

DIR="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read"
GENCODE="/mnt/faruq2/lab_data/reference/gencode/gencode50/gencode.v50.annotation.gtf"
REFSEQ="/mnt/faruq2/lab_data/reference/ncbi_refseq_gtf/GCF_000001405.40_GRCh38.p14_genomic_ucsc_nomenclature_empty_transcripts_removed.gtf"

gffcompare \
  -r "$GENCODE" \
  -o $DIR/annotation/assembled_vs_gencode \
  $DIR/annotation/PPP2R2B_assembled_gencode.gtf

gffcompare \
  -r "$REFSEQ" \
  -o $DIR/annotation/assembled_vs_refseq \
  $DIR/annotation/PPP2R2B_assembled_refseq.gtf
