#!/bin/bash

DIR="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read"
REGION="chr5:146400000-147200000"

gffread \
  -r "$REGION" \
  -T \
  $DIR/assemblies/gencode_transcripts_merged.gtf \
  -o $DIR/annotation/PPP2R2B_assembled_gencode.gtf

gffread \
  -r "$REGION" \
  -T \
  $DIR/assemblies/refseq_transcripts_merged.gtf \
  -o $DIR/annotation/PPP2R2B_assembled_refseq.gtf
