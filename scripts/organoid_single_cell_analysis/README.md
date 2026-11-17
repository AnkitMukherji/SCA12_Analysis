# Organoid single-cell RNA-seq scripts

This workflow preprocesses a filtered organoid scRNA-seq count matrix and suggests cluster identities using ranked marker genes. Prepare the environment from [environment.yaml](environment.yaml), or inspect [requirements.lock.txt](requirements.lock.txt) for the recorded Python dependency versions. The input matrix, features and barcodes are study inputs. Check the folder paths in the scripts before execution.

| Order | Script | Input and output |
| --- | --- | --- |
| 1 | [Preprocessing](01_preprocessing.py) | Count matrix, barcodes and features → preprocessed H5AD with counts, embedding, clusters and differential-expression rankings. |
| 2 | [Manual annotation](02_manual_annotation.py) | Preprocessed H5AD → annotated H5AD, cluster marker evidence, uncertainty and visual summaries. |

From the repository root in the prepared environment:

```bash
python scripts/organoid_single_cell_analysis/01_preprocessing.py
python scripts/organoid_single_cell_analysis/02_manual_annotation.py
```

## 01_preprocessing.py: counts to clusters and DE

### Load and quality control

The matrix is transposed from genes × cells to cells × genes. Barcodes become observation IDs; gene IDs/symbols are retained and duplicate gene-symbol indices made unique. QC uses unnormalized counts. Cells must have >=300 genes, >=500 counts and <=20% mitochondrial counts. Ribosomal/mitochondrial tags do not remove genes.

Scrublet scores retained cells with expected doublet rate 0.06 and seed 29; predicted doublets stay present. Genes detected in fewer than three retained cells are removed. Filtered unnormalized counts are stored in `layers['counts']`.

### Expression, graph and clustering

`X` becomes log1p expression normalized to 10,000 counts per cell. The code flags 3,000 HVGs, computes 50 PCs from them, and builds a graph using 15 neighbours and 30 PCs. Other genes remain available for DE. No explicit scaling, regression or batch integration is performed.

UMAP visualizes the graph. Leiden is run separately at 0.3, 0.6 and 1.0; resolution 0.6 becomes `obs['cluster']`. UMAP does not generate cluster IDs.

### DE and saving

Wilcoxon `rank_genes_groups` compares each selected cluster with the rest using normalized/log expression, tie correction and detection fractions. The full ranking and first 30 entries per cluster are exported. The H5AD retains the ranking for annotation. A manifest records dimensions, settings, versions and hashes. Despite older comments about identity scoring, 01 does not assign cell types.

The script saves the preprocessed AnnData object, marker tables, QC plots and a manifest under the configured results and plots directories.

## 02_manual_annotation.py: ranked DE to broad identity suggestions

### Panels versus programs

Twelve broad identity `PANELS` cover neural progenitors, neurons, astroglial, oligodendroglial, myeloid/microglial, stromal/mesenchymal, endothelial, choroid-plexus-like, non-neural epithelial, endoderm-like, neural-crest-like and pluripotent-like cells.

Ten `PROGRAMS` cover proliferation, neurogenic transition, neuronal maturation, glutamatergic/GABAergic programs, forebrain/posterior/dorsal-midline patterning, stress and hypoxia/glycolysis. Programs are context and never compete as identities. There is no separate core-marker dictionary. The merged `marker_panel.json` is a convenience export; `marker_reference.csv` retains panel kind.

### Eligibility and overlap

The script extracts saved DE without rerunning it. It selects original ranks 1–200 with adjusted p<0.05, log2FC>1 and within-cluster detection>20%, and repeats at 100 and 300. Ranks are not renumbered after filtering.

For every cluster/panel it counts eligible matching genes, lists ranks and divides matches by available panel genes. Identity panels need >=70% availability. The best needs >=3 matches, >=30% panel overlap and >=25% relative separation from the runner-up.

### Uncertainty and propagation

The provisional label must agree at all three cutoffs. Clusters with <50 cells or >10% predicted doublets also remain unresolved. Passing cases receive moderate confidence; other cases receive low confidence and explicit reasons. Evidence index = panel fraction × relative margin, not a probability.

The cluster result is propagated to members, not predicted independently for each cell. No hand-written cluster-number-to-type dictionary overrides the rules. Expression/counts/clusters/UMAP are preserved. A cluster dendrogram is calculated for plotting only.

The separate specificity filter (log2FC>1, >20% inside, <20% outside) supplies a review table/dotplot. It has no p-value criterion and does not determine overlap eligibility; shared lineage markers can be valid even when expressed elsewhere.

Review `cluster_annotations.csv`, `panel_overlap.csv`, `matched_marker_evidence.csv` and `rank_cutoff_sensitivity.csv` together. A best candidate can remain unresolved. Panel/threshold changes should be transparent; obtaining a desired label does not validate it.

The manual stage saves the annotated AnnData object and corresponding evidence tables. Inspect the configured results and plots directories after execution.

## Interpretation

The annotation script assigns provisional labels to clusters. Review the marker panels, doublet flags, cluster sizes and cutoff sensitivity together. Captured cell proportions describe the input sample. Comparisons across conditions require biological replicates and an appropriate design.

[Repository overview](../../README.md)
