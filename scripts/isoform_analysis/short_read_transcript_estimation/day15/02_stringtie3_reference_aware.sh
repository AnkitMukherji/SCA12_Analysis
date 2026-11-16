#!/bin/bash

# Configuration Paths
DIR="/mnt/faruq2/lab_users/vivek/analysis/RNA-Seq/SCA12/batch2/dragen"
OUT="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read"
GENCODE="/mnt/faruq2/lab_data/reference/gencode/gencode50/gencode.v50.annotation.gtf"
REFSEQ="/mnt/faruq2/lab_data/reference/ncbi_refseq_gtf/GCF_000001405.40_GRCh38.p14_genomic_ucsc_nomenclature_empty_transcripts_removed.gtf"
THREADS=20

# Standard Bash arrays for cases and controls
CASES=("1108_NEU1.bam" "1108_NEU2.bam" "1108_NEU3.bam" "1108_NSCP10.bam" "ATC2555_NEU3.bam" "ATC2555_NEU5.bam" "ATC2555_NEU6.bam" "ATC2555_NSCP14.bam" "ATC2555_NSCP8.bam" "150_N2.bam" "150_NEU_5JAN1.bam" "150_NEU_5JAN2.bam" "150_NSC_P14.bam" "150_NSC_P8.bam")
CONTROLS=("GOC29_N3.bam" "GOC29_N4.bam")

# Merge arrays together to iterate through all samples
ALL_SAMPLES=("${CASES[@]}" "${CONTROLS[@]}")

# Loop through each BAM file
for SAMPLE in "${ALL_SAMPLES[@]}";
do
    sample_name=$(basename $SAMPLE .bam)
    mkdir -p "$OUT/assemblies/$sample_name"
    
    # Loop through references using variable expansion
    for REF in "$GENCODE" "$REFSEQ"; do
        if [ "$REF" == "$GENCODE" ]; then
            ref_name="gencode"
        else
            ref_name="refseq"
        fi
        
        # Run StringTie assembly
        stringtie \
            "$DIR/$SAMPLE" \
            -G "$REF" \
            -N \
            -p "$THREADS" \
            -o "$OUT/assemblies/$sample_name/transcripts_${ref_name}.gtf"
    done
done
