###### Creating splice junctions from long-reads using intron-prospector ######

WORK="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_long_read" 
GENOME="/mnt/faruq2/lab_data/reference/human_ref/hg38/dragen_fasta/hg38.fa" 

intron-prospector \
    --genome-fasta="${GENOME}" \
    --intron-bed6="${WORK}/long_read_junctions/pooled_reads.junctions.bed" \
    -C 0.0 \
    "${WORK}/flair/pooled/pooled_reads.bam"
