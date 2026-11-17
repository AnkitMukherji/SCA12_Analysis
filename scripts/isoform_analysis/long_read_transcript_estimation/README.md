# Targeted 5′ RACE long-read workflow

The numbered steps start at demultiplexed FASTQs and a barcode assignment table. Supply the matching genome and reference annotations, confirm software and study paths, and run stages in order from prepared inputs. The two scripts numbered 09 address primer location and reference boundaries. SQANTI3 needs external short-read alignments, junctions and CAGE peaks.

## Script and notebook guide

| File | Purpose |
| --- | --- |
| [`01_merge_barcodes.py`](01_merge_barcodes.py) | Groups compressed reads from barcode folders into a sample FASTQ and tabulates contributing reads. |
| [`02_detect_primer.sh`](02_detect_primer.sh) | Screens both orientations for expected primers or overhangs and writes primer-screening outputs. |
| [`03_primer_matches.sh`](03_primer_matches.sh) | Checks for primer matches with Cutadapt to assess the expected primer sequence in the data. |
| [`04_trim_5race_primers.sh`](04_trim_5race_primers.sh) | Trims the universal primer and cloning overhang and writes trimmed FASTQs and per-run logs. |
| [`05_run_alignment.sh`](05_run_alignment.sh) | Aligns trimmed long reads with minimap2 and writes BAM and alignment QC summaries. |
| [`06_assess_alignment.sh`](06_assess_alignment.sh) | Extracts locus alignments, summarizes eligible primary reads and exports corresponding FASTQ reads. |
| [`07_pool_alignments.sh`](07_pool_alignments.sh) | Merges locus BAMs, indexes the pooled BAM and checks its integrity. |
| [`08_long_read_junctions.sh`](08_long_read_junctions.sh) | Extracts supported introns from the pooled BAM with intron-prospector. |
| [`09_blast_gsp.sh`](09_blast_gsp.sh) | Uses BLAST on the chr5 genome sequence to locate the gene-specific primer. |
| [`09_cap_gtf.py`](09_cap_gtf.py) | Clips target reference transcript boundaries for primer-aware comparison and records the adjustments. |
| [`10_flair_transcriptome.sh`](10_flair_transcriptome.sh) | Runs FLAIR v3.0.1 for each capped annotation and writes reconstructed isoform models and logs. |
| [`11_star_junction_prob.sh`](11_star_junction_prob.sh) | Prepares the short-read STAR junction files used as orthogonal splice support. |
| [`12_prepare_sqanti3_input_data.sh`](12_prepare_sqanti3_input_data.sh) | Stages chr5 references, reconstructed models, short-read inputs and supporting evidence for SQANTI3. |
| [`13_sqanti3_qc.sh`](13_sqanti3_qc.sh) | Runs SQANTI3 on the RefSeq-based reconstructed transcript models and generates QC output. |

[Repository overview](../../../README.md)
