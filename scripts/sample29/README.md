# Sample 29 cerebral organoid analysis

## Set up the environment

With uv, from the repository root on Windows or macOS:

```sh
uv venv --python 3.13 .venv
uv pip install --python .venv -r scripts/sample29/requirements.lock.txt
```

If `.venv` already exists, skip its creation and run the installation command.
Close Python sessions using that environment before updating its packages.
On Windows PowerShell activate with `.\.venv\Scripts\Activate.ps1`.
On macOS activate with `source .venv/bin/activate`.
Select `.venv` as the interpreter when running the script sections in an editor.

From the repository root on Windows or macOS:

```sh
conda env create -f scripts/sample29/environment.yaml
conda activate sample29
```

Alternatively, create and activate a Python 3.13 virtual environment, then
install the same pinned dependencies with
`python -m pip install -r scripts/sample29/requirements.lock.txt`.
Create the environment separately on each machine. The analysis and reference
annotation have been tested on Windows. macOS has not been tested here.

## Run the analysis in sections

Restore the three filtered count files and `29.scRNA.umap.txt` under
`data/single_cell_analysis/test_data/`. Start a fresh Python session in the
Sample 29 environment. In VS Code use **Run Cell** on the numbered `# %%`
sections. In a console started from the repository root, paste complete
statements or blocks in order.

1. `workflow.py` filters cells, normalizes counts, clusters, ranks markers, and saves AnnData.
2. `review.py` applies the manual annotations in `annotation_review.json` and writes the clustering comparisons.
3. `atlas.py` downloads or reuses the public references, then runs the two annotation approaches below. With a completed manual review, start here.

Quality control keeps cells with at least 300 detected genes, at least 500
counts, and at most 20% mitochondrial reads by count. Scrublet uses an assumed
6% doublet rate. Flags are exported and cells are retained. Genes require
detection in three cells. The analysis uses total-count normalization to
10,000, log1p, 3,000 HVGs, 30 PCs for a 15-neighbor graph, and Leiden at
resolutions 0.3, 0.6, and 1.0. Resolution 0.6 is primary. The marker tests are
exploratory one-versus-rest Wilcoxon comparisons. See
[the manual marker sources](MARKER_SOURCES.md).

## Reference annotation

The [Single-cell Best Practices annotation chapter](https://www.sc-best-practices.org/cellular-structure/annotation/)
describes two approaches used here. The reference stage follows both and
reports the evidence separately.

- **Developing brain profile mapping.** Match each Sample 29 cluster mean to
  617 published [Braun et al.](https://doi.org/10.1126/science.adf1226) fetal brain
  profiles by Spearman correlation across 3,216 shared genes. The profile
  matrix was released by the [HNOCA authors](https://github.com/theislab/neural_organoid_atlas).
  The curated `ClusterID (PoolClean)` in Braun's Supplementary Table S2 links
  the profiles to classes. Cortex, other brain regions, and non-neural
  classes participate. This is cluster-level expression mapping, separate
  from HNOCA latent-space projection. The winning class, best competing
  class, score gap, and class agreement across 100 bootstrap resamples are
  reported.
- **Pan-fetal CellTypist classifier.** Apply the official
  [Pan-Fetal Human v2 model](https://www.celltypist.org/models) derived from
  [Suo et al.](https://doi.org/10.1126/science.abo0510) to each cell using
  6,668 shared genes. Rebuild log-normalized expression from the count layer
  after excluding duplicated gene symbols. CellTypist classifies cells
  before any cluster aggregation. The output keeps the two highest scoring
  labels and scores per cell, then reports cluster vote fractions and
  score margins. This model is useful for fetal non-neural lineages and
  has limited neural subtype coverage.

The downloaded reference assets remain under
`data/single_cell_analysis/references/`. Their URLs and hashes are in the
reference manifest. Both references are fetal tissue. Similarity, bootstrap
stability, vote share, and CellTypist scores measure different things and
are not calibrated organoid probabilities. Keep disagreements visible and
interpret broad programs alongside manual markers. Anatomical region names
describe the matched reference profiles.

Outputs are under `results/sample29/reference_annotation/`. Start with
`cluster_reference_comparison.csv` and `interpretation.md`. The marker audit,
cell votes, full per-cell assignments, and top brain reference profiles are
also saved. The former telencephalon-only analysis remains a historical
output under `results/sample29/atlas/`.

## Git and report

Track this directory, including both environment specifications. The Word
report builder `build_report.py` is optional and requires `python-docx`.
It builds the manual analysis report from existing tables and figures.
Keep local environments, downloaded models, reference matrices, and large
query files outside ordinary Git. The reference manifest and compact
comparison table are useful results to version if needed.
