"""Apply documented marker review, validate saved data, and build the analysis report."""

# %% 1. Imports and paths
import json
from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits

threadpool_limits(limits=4)
from sklearn.metrics import adjusted_rand_score

# A terminal session can start at the repository root or scripts/sample29.
working_directory = Path.cwd()
if "__file__" in globals() and Path(__file__).is_file():
    ROOT = Path(__file__).resolve().parents[2]
elif (working_directory / "scripts/sample29").is_dir():
    ROOT = working_directory
elif (
    working_directory.name == "sample29" and working_directory.parent.name == "scripts"
):
    ROOT = working_directory.parents[1]
else:
    raise RuntimeError("Start Python from the repository root or scripts/sample29.")
OUT = ROOT / "results/sample29"
FIG = ROOT / "plots/sample29"
sc.settings.figdir = FIG

# %% 2. Load and apply manual labels
adata = sc.read_h5ad(OUT / "sample29_annotated.h5ad")
annotations = pd.read_csv(
    OUT / "cluster_annotation_automatic.csv", dtype={"cluster": str}
)
annotations["reviewed_label"] = annotations.suggested_label
annotations["review_notes"] = (
    "Provisional marker-panel assignment; examine runner-up and marker expression."
)
overrides_path = ROOT / "scripts/sample29/annotation_review.json"
if overrides_path.exists():
    for cluster, entry in json.loads(overrides_path.read_text()).items():
        annotations.loc[annotations.cluster == cluster, "reviewed_label"] = entry[
            "label"
        ]
        annotations.loc[annotations.cluster == cluster, "review_notes"] = entry[
            "evidence"
        ]
        if "confidence" in entry:
            annotations.loc[annotations.cluster == cluster, "confidence"] = entry[
                "confidence"
            ]
annotations.to_csv(OUT / "cluster_annotation_reviewed.csv", index=False)
mapping = dict(zip(annotations.cluster, annotations.reviewed_label))
adata.obs["cell_type"] = pd.Categorical(adata.obs.cluster.astype(str).map(mapping))
adata.obs["cluster_label"] = pd.Categorical(
    [f"{c}: {mapping[str(c)]}" for c in adata.obs.cluster]
)
sc.pl.umap(adata, color="cluster_label", show=False, save="_reviewed_clusters.png")
plt.close("all")
sc.pl.umap(adata, color="cell_type", show=False, save="_reviewed_cell_types.png")
plt.close("all")
diagnostic = {
    "Epithelial / squamous": ["EPCAM", "ESRP1", "GRHL2", "KRT4", "KRT13", "SCEL"],
    "Choroid plexus": ["TTR", "OTX2", "KCNJ13", "SLC13A4"],
    "Regional neural": [
        "FOXG1",
        "LHX2",
        "EMX2",
        "PAX6",
        "HOXB9",
        "WNT1",
        "PAX3",
        "MSX1",
    ],
    "Endoderm": ["AFP", "APOA1", "HNF4A", "CDX2", "VIL1", "SI"],
}
sc.pl.dotplot(
    adata,
    {k: [g for g in v if g in adata.var_names] for k, v in diagnostic.items()},
    groupby="cluster",
    standard_scale="var",
    show=False,
    save="_review_diagnostic.png",
)
plt.close("all")
sc.pl.violin(
    adata,
    keys=[
        "Proliferation",
        "Forebrain identity",
        "Synaptic maturation",
        "Stress",
        "Hypoxia",
    ],
    groupby="cell_type",
    rotation=90,
    multi_panel=True,
    show=False,
    save="_development_by_cell_type.png",
)
plt.close("all")
neural = adata.obs.cell_type.str.contains("neuron", case=False, regex=False)
if neural.sum() > 0:
    sc.pl.umap(
        adata[neural],
        color=["DCX", "MAP2", "SNAP25", "SYP", "SLC17A7", "GAD1"],
        ncols=3,
        show=False,
        save="_neuronal_cells_development.png",
    )
    plt.close("all")
# %% 3. Summarize composition and compare clusterings
counts = (
    adata.obs.cell_type.value_counts()
    .rename_axis("cell_type")
    .reset_index(name="n_cells")
)
counts["percent"] = counts.n_cells / adata.n_obs * 100
counts.to_csv(OUT / "cell_type_composition.csv", index=False)
vendor = pd.read_csv(
    ROOT / "data/single_cell_analysis/test_data/29.scRNA.umap.txt",
    sep="\t",
    index_col=0,
)
common = adata.obs.index.intersection(vendor.index)
ari = (
    adjusted_rand_score(adata.obs.loc[common, "cluster"], vendor.loc[common, "cluster"])
    if len(common) > 0
    else None
)
if len(common):
    pd.crosstab(adata.obs.loc[common, "cluster"], vendor.loc[common, "cluster"]).to_csv(
        OUT / "vendor_cluster_comparison.csv"
    )
stability = []
for key in ["leiden_0.3", "leiden_1.0"]:
    stability.append(
        {
            "comparison": key,
            "adjusted_rand_index": adjusted_rand_score(
                adata.obs.cluster, adata.obs[key]
            ),
        }
    )
pd.DataFrame(stability).to_csv(OUT / "resolution_comparison.csv", index=False)
# Evaluate singlet-only clustering to identify doublet-sensitive primary clusters.
singlets = adata[~adata.obs.predicted_doublet].copy()
sc.pp.neighbors(singlets, n_neighbors=15, n_pcs=30, random_state=29)
sc.tl.leiden(
    singlets,
    resolution=0.6,
    key_added="singlet_cluster",
    random_state=29,
    flavor="igraph",
    n_iterations=2,
)
pd.crosstab(singlets.obs.cluster, singlets.obs.singlet_cluster).to_csv(
    OUT / "singlet_cluster_comparison.csv"
)
singlet_ari = adjusted_rand_score(singlets.obs.cluster, singlets.obs.singlet_cluster)
adata.obs.to_csv(OUT / "cell_metadata.csv")
adata.write_h5ad(OUT / "sample29_annotated.h5ad", compression="gzip")
assert adata.n_obs == len(pd.read_csv(OUT / "cell_metadata.csv"))
assert (
    adata.layers["counts"].shape == adata.shape
    and np.isfinite(adata.obsm["X_umap"]).all()
)
assert not adata.obs.cell_type.isna().any()
manifest = json.loads((OUT / "run_manifest.json").read_text())

# %% 4. Write the review report
report = f"""# Sample 29: single-cell neuronal organoid analysis

## Dataset and QC

Input: {manifest["input_shape"][0]:,} cells and {manifest["input_shape"][1]:,} features. Retained: {adata.n_obs:,} cells and {adata.n_vars:,} genes. Removed {manifest["input_shape"][0] - adata.n_obs:,} cells using >=300 genes, >=500 counts and <=20% mitochondrial counts. Counts are saved in the AnnData counts layer; X contains log1p expression after normalization to 10,000 counts per cell.

Scrublet flagged {manifest["doublets_flagged"]} possible doublets using a 6% expected rate. They remain in the primary analysis because loading metadata is unknown; flagged cells should be reviewed before downstream biological inference. Re-clustering the unflagged cells with the existing PCA gives ARI {singlet_ari:.3f} against their primary labels. This is a sensitivity check, not an independent annotation.

## Cluster annotations

Every cell has a cluster and provisional cell-type annotation. Reviewed assignments and evidence:

```csv
{annotations[["cluster", "n_cells", "reviewed_label", "confidence", "review_notes"]].to_csv(index=False)}```

Cell composition:

```csv
{counts.to_csv(index=False)}```

## Evidence and interpretation

The reviewed marker evidence indicates a heterogeneous sample with substantial epithelial and mesenchymal off-target populations, alongside neural progenitors and an immature neuronal cluster (cluster 11; 553 cells, 6.42% of retained cells). Squamous epithelial and hepatointestinal-like programs argue against treating all clusters as neuronal. Neural progenitors show differing regional programs, including forebrain-like and posterior-patterned signatures. FOXG1 is weak across the diagnostic panel, so the LHX2/EMX2-supported forebrain-like label should not be read as confirmed telencephalic identity. These observations do not establish a uniformly specified or mature neuronal organoid. Cluster 6 is strongly doublet-enriched and clusters 8 and 10 have additional size/state uncertainty.

See marker_expression_by_cluster.csv for gene detection fractions and mean log-normalized expression, cluster_markers_top30.csv for ranked markers, and cluster_identity_scores.csv for competing identities. Score margins and confidence labels are heuristic. Mixed developmental programs can occur within one lineage; cycling state is not a separate lineage. Shared markers such as VIM, SOX2, TTR and S100B cannot establish an identity alone.

Plots include progenitor/proliferation programs, forebrain identity, deep and upper cortical layer genes, synaptic maturation, stress and hypoxia. Positive neuronal and synaptic expression supports differentiation at the RNA level. It does not confirm normal development or physiological maturity. Assess expected identity against organoid age and induction protocol, then validate using spatial staining and electrophysiology. One organoid sample cannot establish reproducibility or donor-level effects. No trajectory is interpreted as developmental time.

## Clustering and reproducibility

Three Leiden resolutions are exported; primary resolution 0.6. Vendor clusters were compared only after independent clustering: {len(common)} matching cells, ARI {ari}. Resolution comparisons and vendor contingency tables are supplied. UMAP geometry alone does not define cell identity or lineage relationships.

Set up the pinned environment using scripts/sample29/README.md. Run workflow.py and review.py in order. Run atlas.py separately for reference annotation. See scripts/sample29/README.md for parameters, assumptions and references. run_manifest.json records input hashes, seed and versions; requirements.lock.txt records the complete package environment. Downloaded package cache and installed environment are under environments/. Existing project files and vendor analysis were preserved.

## Files

- sample29_annotated.h5ad: reviewed cell labels, raw counts layer, normalized expression, HVGs, PCA, neighbors, UMAP, clusters, scores and marker statistics.
- cell_metadata.csv and umap_coordinates.csv: cell-level exports.
- cluster_annotation_reviewed.csv and cell_type_composition.csv: annotations and composition.
- qc_all_cells.csv, qc_doublet_flags.csv: QC audit.
- plots/sample29/: QC, clustering, annotation, gene maps and marker dotplots.

## Sources

- [Scanpy preprocessing and clustering](https://scanpy.readthedocs.io/en/latest/tutorials/basics/clustering.html)
- [Scanpy Scrublet](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.scrublet.html)
- [Fair et al. 2020, cerebral organoid maturation](https://pmc.ncbi.nlm.nih.gov/articles/PMC7562943/)
- [Quadrato et al. 2017, organoid cell diversity](https://www.nature.com/articles/nature22047)
"""
(OUT / "analysis_report.md").write_text(report, encoding="utf-8")
print("Review and saved-object validation complete", flush=True)
