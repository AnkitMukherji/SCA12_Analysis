DIR="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_long_read/reference"
REF="/mnt/faruq2/lab_data/reference/human_ref/hg38/dragen_fasta/hg38.fa"

samtools faidx $REF chr5 > "$DIR/chr5.fa"

blastn \
    -task blastn-short \
    -query "$DIR/PPP2R2B_GSP.fa" \
    -subject "$DIR/chr5.fa" \
    -strand both \
    -dust no \
    -outfmt \
      "6 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore" \
    > "$DIR/PPP2R2B_GSP_hits.tsv"
