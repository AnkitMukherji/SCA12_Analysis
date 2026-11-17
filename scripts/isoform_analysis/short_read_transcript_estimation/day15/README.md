# Day 15 short-read PPP2R2B assembly

Supply aligned and indexed study BAMs, compatible GENCODE and RefSeq GTFs, the two PPP2R2B reference catalogs and STAR junction tables. Confirm configured BAM inputs before assembly. Run the numbered stages in order, with 04a and 04b collecting models at different levels. GTF and junction coordinates must refer to the same genome build.

## Script and notebook guide

| File | Purpose |
| --- | --- |
| [`01_extract_PPP2R2B_bams.sh`](01_extract_PPP2R2B_bams.sh) | Extracts locus BAMs and indexes the outputs with SAMtools. |
| [`02_stringtie3_reference_aware.sh`](02_stringtie3_reference_aware.sh) | Runs reference-guided StringTie assembly separately with GENCODE and RefSeq. |
| [`03_stringtie3_merge.sh`](03_stringtie3_merge.sh) | Builds per-reference input lists and merges the assemblies. |
| [`04a_extract_PPP2R2B_from_each_transcript_models.sh`](04a_extract_PPP2R2B_from_each_transcript_models.sh) | Extracts PPP2R2B transcript models from individual assemblies. |
| [`04b_extract_PPP2R2B_from_merged_transcript_models.sh`](04b_extract_PPP2R2B_from_merged_transcript_models.sh) | Extracts PPP2R2B models from the merged catalog for each reference. |
| [`05_compare_ppp2r2b_transcript_models_with_reference.sh`](05_compare_ppp2r2b_transcript_models_with_reference.sh) | Uses GffCompare to compare assembled models against each reference annotation. |
| [`06_classify_isoforms.py`](06_classify_isoforms.py) | Classifies exon and intron chains relative to both reference catalogs. |
| [`07_extract_junctions.py`](07_extract_junctions.py) | Reads STAR splice-junction tables and summarizes junction support for the locus. |

[Repository overview](../../../../README.md)
