# PPP2R2B reference catalogs

Use compatible GRCh38 GTFs for GENCODE and RefSeq. Run reference_catalog.py separately for both sources, passing an input GTF, source label and output prefix. Then review the catalog filenames expected by the two pandas scripts.

## Script and notebook guide

| File | Purpose |
| --- | --- |
| [`compare_annotations.py`](compare_annotations.py) | Compares annotated transcript pairs by exact exon chain, shared intron chain and partial intron overlap. |
| [`compare_first_exons.py`](compare_first_exons.py) | Summarizes first exons from the reference catalogs for start-site comparisons. |
| [`reference_catalog.py`](reference_catalog.py) | Extracts locus-specific GTF records and writes a transcript TSV with coordinates, exon and intron chains, first exon and coding segments. |

[Repository overview](../../../README.md)
