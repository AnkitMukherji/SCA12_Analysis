BAM_DIR="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_long_read/alignment"
PPP2R2B_REGION="chr5:146300000-147300000"

# Extract PPP2R2B locus
mkdir -p ${BAM_DIR}/ppp2r2b_bam
for BAM in ${BAM_DIR}/*.bam
do

    SAMPLE=$(basename ${BAM} .bam)

    samtools view \
        -bh \
        ${BAM} \
        "${PPP2R2B_REGION}" \
        > ${BAM_DIR}/ppp2r2b_bam/${SAMPLE}.PPP2R2B.bam

    samtools index \
        ${BAM_DIR}/ppp2r2b_bam/${SAMPLE}.PPP2R2B.bam

done

# Count PPP2R2B-supporting reads
mkdir -p ${BAM_DIR}/qc
for BAM in ${BAM_DIR}/ppp2r2b_bam/*.bam
do

    SAMPLE=$(basename ${BAM} .PPP2R2B.bam)

    COUNT=$(samtools view -c -F 2308 ${BAM}) # -F 2308 removes unmapped, secondary and supplementary reads for a clean primary read count

    echo -e "${SAMPLE}\t${COUNT}"

done > ${BAM_DIR}/qc/PPP2R2B_mapped_reads.tsv

# Extract the actual reads for PPP2R2B
mkdir -p ${BAM_DIR}/fastq
for BAM in ${BAM_DIR}/ppp2r2b_bam/*.bam
do

    SAMPLE=$(basename ${BAM} .PPP2R2B.bam)

    samtools fastq \
        -F 2308 \
        ${BAM} \
        | gzip \
        > ${BAM_DIR}/fastq/${SAMPLE}.PPP2R2B.fastq.gz

done
