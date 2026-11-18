#!/usr/bin/env bash
set -euo pipefail

DIR="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_long_read"
GENOME="/mnt/faruq2/lab_data/reference/human_ref/hg38/dragen_fasta/hg38.fa"
GTF="$DIR/reference/gencode.v50.annotation.primer_capped.gtf"
REFSEQ="$DIR/reference/GCF_000001405.40_GRCh38.p14_genomic_ucsc_nomenclature_empty_transcripts_removed.primer_capped.gtf"
BAM="$DIR/flair/pooled/pooled_reads.bam"
JUNCTIONS="$DIR/short_read_junctions"
MERGED_JUNCTIONS="$JUNCTIONS/merged"
OUTDIR="$DIR/flair/primary"

mkdir -p "$OUTDIR" "$MERGED_JUNCTIONS"

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
# Prepare short-read RNA-seq splice junctions
# ------------------------------------------------------------

# =============================================================================
# Function: merge STAR SJ.out.tab files
#
# Junctions are grouped by:
#   chr, intron_start, intron_end, strand, motif
#
# For each junction:
#   column 6 = annotated if annotated in ANY replicate
#   column 7 = SUM unique-mapping read support
#   column 8 = SUM multi-mapping read support
#   column 9 = MAX overhang observed across replicates
# =============================================================================

combine_junctions() {

    OUTPUT=$1
    shift

    echo "Creating merged junction file:"
    echo "  ${OUTPUT}"

    VALID_FILES=()

    for SJ in "$@"
    do
        if [[ -s "${SJ}" ]]; then
            echo "  + ${SJ}"
            VALID_FILES+=("${SJ}")
        else
            echo "WARNING: Missing/empty junction file:"
            echo "         ${SJ}"
        fi
    done

    if (( ${#VALID_FILES[@]} == 0 )); then
        echo "ERROR: No valid SJ.out.tab files found."
        rm -f "${OUTPUT}"
        return 1
    fi

    awk 'BEGIN {OFS="\t"}
    {
        key=$1 OFS $2 OFS $3 OFS $4 OFS $5

        if (!(key in annotated))
            annotated[key]=0

        # annotated status: 1 if annotated in any replicate
        if ($6 > annotated[key])
            annotated[key]=$6

        # sum read support across replicates
        unique[key]+=$7
        multi[key]+=$8

        # retain maximum observed overhang
        if ($9 > overhang[key])
            overhang[key]=$9
    }
    END {
        for (key in unique)
            print key, annotated[key], unique[key], multi[key], overhang[key]
    }' "${VALID_FILES[@]}" \
    | sort -k1,1 -k2,2n -k3,3n \
    > "${OUTPUT}"

    echo "  Junctions: $(wc -l < "${OUTPUT}")"
    echo
}

# =============================================================================
# GOC29_NEURON
# =============================================================================

combine_junctions \
    "${MERGED_JUNCTIONS}/GOC29_NEURON.SJ.out.tab" \
    "${JUNCTIONS}/NEU-GOC-29-SCA12.SJ.out.tab" \
    "${JUNCTIONS}/GOC29_N3.SJ.out.tab" \
    "${JUNCTIONS}/GOC29_N4.SJ.out.tab"


# =============================================================================
# 150_NEURON
# =============================================================================

combine_junctions \
    "${MERGED_JUNCTIONS}/150_NEURON.SJ.out.tab" \
    "${JUNCTIONS}/150_N2.SJ.out.tab" \
    "${JUNCTIONS}/150_NEU_5JAN1.SJ.out.tab" \
    "${JUNCTIONS}/150_NEU_5JAN2.SJ.out.tab"


# =============================================================================
# 150_NSC
# =============================================================================

combine_junctions \
    "${MERGED_JUNCTIONS}/150_NSC.SJ.out.tab" \
    "${JUNCTIONS}/150_NSC_P14.SJ.out.tab" \
    "${JUNCTIONS}/150_NSC_P8.SJ.out.tab"


# =============================================================================
# AT2555_NEURON
# =============================================================================

combine_junctions \
    "${MERGED_JUNCTIONS}/AT2555_NEURON.SJ.out.tab" \
    "${JUNCTIONS}/ATC2555_NEU3.SJ.out.tab" \
    "${JUNCTIONS}/ATC2555_NEU5.SJ.out.tab" \
    "${JUNCTIONS}/ATC2555_NEU6.SJ.out.tab"


# =============================================================================
# AT2555_NSC
# =============================================================================

combine_junctions \
    "${MERGED_JUNCTIONS}/AT2555_NSC.SJ.out.tab" \
    "${JUNCTIONS}/ATC2555_NSCP14.SJ.out.tab" \
    "${JUNCTIONS}/ATC2555_NSCP8.SJ.out.tab"
    
cat \
  "${MERGED_JUNCTIONS}/GOC29_NEURON.SJ.out.tab" \
  "${MERGED_JUNCTIONS}/150_NEURON.SJ.out.tab" \
  "${MERGED_JUNCTIONS}/150_NSC.SJ.out.tab" \
  "${MERGED_JUNCTIONS}/AT2555_NEURON.SJ.out.tab" \
  "${MERGED_JUNCTIONS}/AT2555_NSC.SJ.out.tab" \
  > "${MERGED_JUNCTIONS}/all_samples.SJ.out.tab"

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
        --junction_tab "${MERGED_JUNCTIONS}/all_samples.SJ.out.tab" \
        --junction_support 2 \
        --ss_window 15 \
        --support 5 \
        --stringent \
        --check_splice \
        --end_window 75 \
        --max_ends 1 \
        --filter default \
        -t 24 \
        -o "$OUT/pooled_reads_transcriptome" \
        2>&1 | tee "$OUT/flair_transcriptome.log"

    echo "Completed: $name"

done

echo "All FLAIR reconstructions completed."
