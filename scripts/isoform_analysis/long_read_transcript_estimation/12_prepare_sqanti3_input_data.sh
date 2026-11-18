#!/bin/bash
set -eou pipefail

DIR="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_long_read"

mkdir -p \
    "$DIR/sqanti3/input_data/references" \
    "$DIR/sqanti3/input_data/flair_isoforms" \
    "$DIR/sqanti3/input_data/bams" \
    "$DIR/sqanti3/input_data/short_read_junctions" \
    "$DIR/sqanti3/input_data/orthogonal_data"

# Extract chromosome 5 from human ref
REF="/mnt/faruq2/lab_data/reference/human_ref/hg38/dragen_fasta/hg38.fa"
OUTREF="$DIR/sqanti3/input_data/references/hg38.chr5.fa"

samtools faidx "$REF" chr5 > "$OUTREF"
samtools faidx "$OUTREF"

# Extract chromosome 5 from gtf
GENCODE="$DIR/reference/gencode.v50.annotation.primer_capped.gtf"
REFSEQ="$DIR/reference/GCF_000001405.40_GRCh38.p14_genomic_ucsc_nomenclature_empty_transcripts_removed.primer_capped.gtf"

awk 'BEGIN {FS=OFS="\t"} /^#/ {print; next} $1=="chr5" {print}' "$GENCODE" > "$DIR/sqanti3/input_data/references/gencode.chr5.gtf"
awk 'BEGIN {FS=OFS="\t"} /^#/ {print; next} $1=="chr5" {print}' "$REFSEQ" > "$DIR/sqanti3/input_data/references/refseq.chr5.gtf"

# Extract chromosome 5 from FLAIR isoform gtf
ISOGTF="$DIR/flair/primary/RefSeq/pooled_reads_transcriptome.isoforms.gtf"

awk 'BEGIN {FS=OFS="\t"} /^#/ {print; next} $1=="chr5" {print}' "$ISOGTF" > "$DIR/sqanti3/input_data/flair_isoforms/flair_isoforms.chr5.gtf"

# Extract chromosome 5 from dragen bams of short read rna-seq data
BATCH1_DIR="/mnt/faruq2/lab_users/vivek/analysis/RNA-Seq/SCA12/batch1/dragen"
BATCH2_DIR="/mnt/faruq2/lab_users/vivek/analysis/RNA-Seq/SCA12/batch2/dragen"

declare -A SAMPLE_DIRS=(
    # Batch 1 sample
    [NEU-GOC-29-SCA12]="$BATCH1_DIR"

    # Batch 2 samples
    [150_N2]="$BATCH2_DIR"
    [150_NEU_5JAN1]="$BATCH2_DIR"
    [150_NEU_5JAN2]="$BATCH2_DIR"
    [150_NSC_P14]="$BATCH2_DIR"
    [150_NSC_P8]="$BATCH2_DIR"
    [ATC2555_NEU3]="$BATCH2_DIR"
    [ATC2555_NEU5]="$BATCH2_DIR"
    [ATC2555_NEU6]="$BATCH2_DIR"
    [ATC2555_NSCP14]="$BATCH2_DIR"
    [ATC2555_NSCP8]="$BATCH2_DIR"
    [GOC29_N3]="$BATCH2_DIR"
    [GOC29_N4]="$BATCH2_DIR"
)

THREADS=20

for sample in "${!SAMPLE_DIRS[@]}";
do

    BAM="${SAMPLE_DIRS[$sample]}/${sample}.bam"
    FINAL="$DIR/sqanti3/input_data/bams/${sample}.chr5.bam"

    echo "Processing: $sample"

    [[ -s "$BAM" ]] || {
        echo "Missing BAM: $BAM" >&2
        exit 1
    }

    samtools view -@ "$THREADS" -b -o "$FINAL" "$BAM" chr5
    samtools quickcheck -v "$FINAL"
    samtools index -@ "$THREADS" "$FINAL"

    echo -n "Reads retained: "
    samtools view -c "$FINAL"

    echo "Completed: $sample"
    echo

done

# Extract chromosome 5 from STAR splice-junction files
for sample in "${!SAMPLE_DIRS[@]}"; 
do

    INPUT="${SAMPLE_DIRS[$sample]}/${sample}.SJ.out.tab"
    OUTPUT="$DIR/sqanti3/input_data/short_read_junctions/${sample}.SJ.out.tab"
    
    echo "Processing junctions: $sample"
    [[ -f "$INPUT" ]] || { echo "Missing junction file: $INPUT" >&2; exit 1; }
    
    awk '
    BEGIN {FS=OFS="\t"}
    $1=="chr5" && NF == 9 && $4 ~ /^[012]$/ {
        print $0;
        retained++
        next
    }
    { removed++ }
    END {
        print "Retained:", retained+0 > "/dev/stderr"
        print "Removed:", removed+0 > "/dev/stderr"
    }' "$INPUT" > "$OUTPUT"
done

# Extract chromosome 5 from CAGE peaks
extract_chr5_bed() {

    local INPUT="$1"
    local OUTPUT="$2"

    awk '
    BEGIN {FS=OFS="\t"}

    /^#/ {next}
    /^track/ {next}
    /^browser/ {next}

    $1=="chr5" {print}' "$INPUT" > "$OUTPUT"

    echo "Output: $OUTPUT"
    wc -l "$OUTPUT"
}

CAGE="/mnt/faruq2/lab_data/tools/SQANTI3/human.refTSS_v3.1.hg38.bed"
extract_chr5_bed "$CAGE" "$DIR/sqanti3/input_data/orthogonal_data/CAGE.chr5.bed"
