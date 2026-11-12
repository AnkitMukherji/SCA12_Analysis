# Sample29 analysis scripts

Sequential analysis of day-40 unguided cerebral-organoid scRNA-seq. Scripts use ordered blocks without custom function definitions. Documentation reflects current code on 12 November 2026.

## Execution order and environments

Run from `D:/sca12_analysis`: scripts use `Path.cwd()` as the root. In an editor, choose the matching interpreter and run imports, paths, loading and subsequent blocks in order; later blocks depend on earlier variables.

| Script | Environment | Main input | Saved H5AD under results/sample29 |
|---|---|---|---|
| `01_preprocessing.py` | `environments/sample29` | Compressed matrix/barcodes/features in data/single_cell_analysis/sample29 | `intermediate_data/sample29_preprocessed.h5ad` |
| `02_manual_annotation.py` | Same sample29 environment | Preprocessed object and its saved rank_genes_groups | `intermediate_data/sample29_manual_annotation.h5ad` |
| `03_automated_annotation.py` | `environments/sample29_reference` (Python 3.11) | Preprocessed object, references; optional manual object for comparison | `automated_annotation/sample29_automated_annotation.h5ad` |

```powershell
Set-Location D:\sca12_analysis
& .\environments\sample29\Scripts\python.exe .\scripts\sample29\01_preprocessing.py
& .\environments\sample29\Scripts\python.exe .\scripts\sample29\02_manual_annotation.py
& .\environments\sample29_reference\Scripts\python.exe -u .\scripts\sample29\03_automated_annotation.py
```

Run only the required stage: annotation reuses completed preprocessing. Named outputs are overwritten, but stale files are not automatically cleared. Figures are written to disk with a non-interactive Matplotlib backend.

`environment.yaml` specifies the sample29 Conda/Python environment with pinned pip dependencies; `requirements.lock.txt` lists its pip versions. They are not the reference-model environment specification. Reference versions/requirements are recorded under automated results. The root requirements file contains machine-specific build paths and is not a portable specification for these scripts.

## 01_preprocessing.py: counts to clusters and DE

### Load and quality control

The matrix is transposed from genes × cells to cells × genes. Barcodes become observation IDs; gene IDs/symbols are retained and duplicate gene-symbol indices made unique. QC uses unnormalized counts. Cells must have >=300 genes, >=500 counts and <=20% mitochondrial counts. Ribosomal/mitochondrial tags do not remove genes.

Scrublet scores retained cells with expected doublet rate 0.06 and seed 29; predicted doublets stay present. Genes detected in fewer than three retained cells are removed. Filtered unnormalized counts are stored in `layers['counts']`.

### Expression, graph and clustering

`X` becomes log1p expression normalized to 10,000 counts per cell. The code flags 3,000 HVGs, computes 50 PCs from them, and builds a graph using 15 neighbours and 30 PCs. Other genes remain available for DE. No explicit scaling, regression or batch integration is performed.

UMAP visualizes the graph. Leiden is run separately at 0.3, 0.6 and 1.0; resolution 0.6 becomes `obs['cluster']`. UMAP does not generate cluster IDs.

### DE and saving

Wilcoxon `rank_genes_groups` compares each selected cluster with the rest using normalized/log expression, tie correction and detection fractions. The full ranking and first 30 entries per cluster are exported. The H5AD retains the ranking for annotation. A manifest records dimensions, settings, versions and hashes. Despite older comments about identity scoring, 01 does not assign cell types.

The [preprocessing inventory](../../results/sample29/preprocessing/README.md) explains every output, its export stage and the AnnData fields.

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

The [manual output inventory](../../results/sample29/manual_annotation/README.md) explains every file, reason and plot. The manual H5AD is saved in intermediate_data.

## 03_automated_annotation.py: reference-based comparison

03 reads the original preprocessed object, aligns counts to HNOCA features, adapts the published scPoli model, and transfers hierarchy paths with shared-neighbour Jaccard-square weights. Braun developing-brain CellTypist is run separately on log1p expression normalized to 10,000 counts across all genes.

Raw predictions, uncertainty and supported parent/deeper labels are retained. Vote, margin, distance, doublet and hierarchy checks control acceptance. Previous manual labels supply comparison/coverage review, not training labels. `automated_review_label` is the main review output, distinct from raw and method-supported predictions.

HNOCA lacks explicit general non-neural epithelial/endoderm labels in its transferred taxonomy. Coverage conflicts are withheld in the main view. Regional names are hypotheses, not validated anatomy. See the existing [automated README](../../results/sample29/automated_annotation/README.md), [reference comparison](../../results/sample29/automated_annotation/reference_selection.md) and [analysis summary](../../results/sample29/automated_annotation/analysis_summary.md).

## Reproducibility

Seeds/thread limits control randomness/resources but cannot guarantee identical results across package/platform changes. Read saved settings alongside current code: edits do not update old outputs until execution. Changing 01 inputs/clustering requires regenerating dependent annotations. This documentation update does not rerun analyses or alter biological decisions.
