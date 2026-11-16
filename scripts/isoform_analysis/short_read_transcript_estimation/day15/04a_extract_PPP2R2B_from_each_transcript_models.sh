#!/bin/bash

DIR="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read/day15/assemblies"
REFSEQ_LIST="$DIR/refseq_assembly_list.txt"
GENCODE_LIST="$DIR/gencode_assembly_list.txt"
REGION="chr5:146400000-147200000"

for i in $(cat "$REFSEQ_LIST" "$GENCODE_LIST");
do
    dir=$(dirname "$i")
    
    if [ "$(basename "$i")" == "transcripts_gencode.gtf" ];
    then
        gffread \
            -r "$REGION" \
            -T \
            "$i" \
            -o "$dir/PPP2R2B_transcripts_gencode.gtf"
    else
        gffread \
            -r "$REGION" \
            -T \
            "$i" \
            -o "$dir/PPP2R2B_transcripts_refseq.gtf"
    fi
done
