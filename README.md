# SCA12 analysis: clinical phenotypes, transcript models, and organoid scRNA-seq

This repository contains research scripts for spinocerebellar ataxia type 12 (SCA12). The clinical scripts describe phenotype and repeat-size patterns. The bulk RNA workflows reconstruct transcript models at the **PPP2R2B** locus and compare their exon chains, splice junctions, and transcript ends with reference annotations. A separate organoid workflow performs exploratory single-cell preprocessing and marker-based annotation.

[Repository map](#repository-map) · [Clinical metadata](#clinical-metadata) · [Reference catalogs](#reference-catalogs) · [Short-read reconstruction](#short-read-reconstruction) · [Long-read reconstruction](#long-read-reconstruction) · [Organoid single-cell workflow](#organoid-single-cell-workflow) · [Reproduction and interpretation](#reproduction-and-interpretation)

## Folder guides

- [scripts](scripts/README.md)

## Repository map

| Area | Entry point | Required inputs and key outputs |
| --- | --- | --- |
| Clinical summaries | [Metadata plotting and reporting](scripts/metadata_viz/) | Study workbooks → figures, Sankey HTML/PNG, optional Word report |
| Reference comparison | [Reference catalog scripts](scripts/isoform_analysis/compare_references/) | GRCh38 GTFs → PPP2R2B GTFs, transcript catalogs and comparisons |
| Short-read reconstruction | [Day 15](scripts/isoform_analysis/short_read_transcript_estimation/day15/) and [Day 30](scripts/isoform_analysis/short_read_transcript_estimation/day30/) | Aligned RNA-seq BAMs, reference annotations, and STAR splice-junction tables → assembled transcript models, reference classifications and junction evidence |
| Long-read reconstruction | [5′ RACE workflow](scripts/isoform_analysis/long_read_transcript_estimation/) | Barcode assignment table, long-read FASTQs, genome/annotations and short-read support → aligned reads, reconstructed isoforms and QC reports |
| Organoid single-cell analysis | [Preprocessing and manual annotation](scripts/organoid_single_cell_analysis/) | Filtered count matrix, barcodes and features → AnnData objects, cluster marker rankings, annotation evidence and plots |

## Clinical metadata

Plots cover demographics, allele sizes, onset-age associations, disease-duration associations, concordance among clinical scales, and phenotype frequencies. The first workbook supplies demographics, repeat sizes and clinical scores. The second supplies categorical phenotypes for the proportions plots and Sankey. Run the relevant plotting modules before the report generator if you want to include their figures. See the [metadata guide](scripts/metadata_viz/README.md) for workbook names, fields, output names and the known path and report-image mismatches. Some paths are resolved from the scripts themselves and require review even when commands are launched at the repository root.

Report the denominator for each analysis after filtering missing values. The repeat-size/onset-age analysis estimates a cross-sectional association. Measuring intergenerational anticipation requires family comparisons, and measuring progression within individuals requires longitudinal observations. The Sankey shows co-occurring categories.

## Reference catalogs

The reference scripts compare **GENCODE v50** and **NCBI RefSeq GRCh38.p14** annotations. The catalog builder extracts PPP2R2B features from a GTF, writes a target-only GTF and a `.catalog.tsv`, and records exon, intron and CDS coordinates. It accounts for strand when identifying the first exon. The remaining scripts compare exon and intron chains and summarize first exons.

The catalog builder takes three positional arguments. Use a source label and output prefix appropriate to the chosen GTF:

```bash
python scripts/isoform_analysis/compare_references/reference_catalog.py input.gtf GENCODE gencode_ppp2r2b
```

Build both catalogs, then check the file paths in the comparison scripts before running them with pandas. Keep chromosome naming consistent and preserve GTF's 1-based inclusive coordinates.

## Short-read reconstruction

The Day 15 and Day 30 directories each contain the following stages:

1. Extract and index BAMs for the PPP2R2B locus using SAMtools.
2. Use StringTie for reference-guided per-sample assembly against GENCODE and RefSeq. Check the configured BAM inputs for each stage.
3. Merge assemblies separately for each annotation, then extract PPP2R2B models from individual and merged GTFs using the `04a` and `04b` scripts.
4. Compare the locus models with the corresponding references using GffCompare.
5. Classify assembled models using both reference catalogs. Categories distinguish exact exon chains, known intron chains with different termini, novel junctions, new combinations of known junctions, and unresolved single-exon models.
6. Use the junction extractor to summarize supporting reads from STAR `SJ.out.tab` files. Read the resulting junction support alongside the transcript classifications.

## Long-read reconstruction

| Stage | What the scripts do |
| --- | --- |
| 01–04 | Merge barcode FASTQs by sample, screen for primer matches, and trim the universal primer and cloning overhang with Cutadapt. Review primer orientation and trimming logs. |
| 05–08 | Align with minimap2, check and subset alignments to the locus with SAMtools, pool alignments, and extract supporting introns with intron-prospector. |
| 09 | Use BLAST to locate the gene-specific primer and the GTF utility to adjust reference transcript boundaries for primer-aware comparison. |
| 10 | Run FLAIR v3.0.1 separately against the two capped annotation sets. The script checks required inputs, BAM integrity, and junction BED fields before reconstructing transcript models. |
| 11–13 | Prepare STAR junction support, chromosome-level reference and short-read inputs, and run SQANTI3 QC on the RefSeq-based reconstructed models with CAGE peaks as orthogonal evidence for transcription starts. |

The long-read scripts depend on Bash, Python, Cutadapt, minimap2, SAMtools, intron-prospector, BLAST+, FLAIR and SQANTI3. The short-read scripts additionally use StringTie, GffCompare and STAR junction outputs. Review each script's paths, sample lists, reference files, outputs and tool options in sequence. The scripts are documented entry points for their original environment, with adaptation required for a new installation.

## Organoid single-cell workflow

Use the [environment specification](scripts/organoid_single_cell_analysis/environment.yaml) and [dependency record](scripts/organoid_single_cell_analysis/requirements.lock.txt). The two scripts use the current working directory as the project root. Supply the study count matrix, barcode and feature files in the layout expected by the preprocessing script. From the repository root in the prepared environment:

```bash
python scripts/organoid_single_cell_analysis/01_preprocessing.py
python scripts/organoid_single_cell_analysis/02_manual_annotation.py
```

Preprocessing reads the sparse count matrix and applies cell and gene filters. It records Scrublet doublet flags, normalizes and log-transforms counts, selects variable genes, and constructs a PCA-based neighbor graph for Leiden clustering and UMAP. Wilcoxon-ranked cluster markers are saved alongside the preprocessed AnnData object. Manual annotation reuses those rankings, compares them with identity panels, and records candidate labels, matched-marker evidence and unresolved clusters. Candidate identities are cluster-level suggestions. Program markers such as proliferation or stress help describe cell state.

The current checkout provides preprocessing and manual-annotation scripts. The older guide also describes a local automated-reference workflow and result inventories. Use the repository map above to locate the files available here.

## Reproduction and interpretation

Study inputs and generated result directories are external to this checkout. The root requirements file includes machine-specific package references. Use the subworkflow's recorded environment where supplied and check versions and mounted paths before execution. The commands above describe code entry points. This documentation review did not rerun the analyses.

Compare putative transcript structures against **both** annotations and inspect splice support, transcript ends, and assay-specific biases. Short reads support individual junctions more directly than complete exon chains. Targeted long reads can connect distant junctions, while the primer and alignment strategy affect which ends are observable. A transcript structure's difference from an annotation identifies a candidate for follow-up. Promoter activity, repeat-dependent effects and allele-specific expression require additional supporting evidence.

Organoid cluster markers support provisional labels. Cell proportions describe captured cells, and inference across conditions requires suitable biological replicates. The registry and organoid scripts operate on separate inputs.
