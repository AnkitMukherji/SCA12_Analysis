# Imports, paths and marker panels
import os

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
import json
import hashlib
import platform
import gzip
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.io
import scanpy as sc
import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits

threadpool_limits(limits=4)

working_directory = Path.cwd()

ROOT = working_directory
DATA = ROOT / "data/single_cell_analysis/sample29"
OUT = ROOT / "results/sample29"
OUT_QC = ROOT / "results/sample29/qc"
OUT_MAN_ANNOT = ROOT / "results/sample29/manual_annotation"
FIG_QC = ROOT / "plots/sample29/qc"
FIG_MAN_ANNOT = ROOT / "plots/sample29/manual_annotation"
OUT_QC.mkdir(parents=True, exist_ok=True)
OUT_MAN_ANNOT.mkdir(parents=True, exist_ok=True)
FIG_QC.mkdir(parents=True, exist_ok=True)
FIG_MAN_ANNOT.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = FIG_QC
sc.settings.figdir = FIG_MAN_ANNOT
sc.settings.n_jobs = 4
sc.set_figure_params(dpi=120, dpi_save=180, facecolor="white")
SEED = 29
PANELS = {
    "Radial glia / neural progenitors": ["SOX2", "PAX6", "NES", "VIM", "HES1", "HES5", "FABP7"],
    "Intermediate progenitors": ["EOMES", "NEUROG2", "INSM1", "HES6"],
    "Immature neurons": ["DCX", "TUBB3", "STMN2", "NEUROD1", "NEUROD2"],
    "Excitatory neurons": ["SLC17A6", "SLC17A7", "TBR1", "BCL11B", "SATB2", "NEUROD6"],
    "Inhibitory neurons": ["GAD1", "GAD2", "SLC32A1", "DLX1", "DLX2", "DLX5"],
    "Astroglia": ["GFAP", "AQP4", "SLC1A3", "ALDH1L1", "S100B"],
    "Oligodendrocyte lineage": ["OLIG1", "OLIG2", "SOX10", "PDGFRA", "MBP", "PLP1"],
    "Microglia": ["AIF1", "TYROBP", "P2RY12", "CX3CR1", "C1QA"],
    "Mesenchymal / fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"],
    "Mural / smooth muscle": ["ACTA2", "TAGLN", "MYL9", "RGS5", "PDGFRB"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "EMCN", "CLDN5"],
    "Choroid plexus / epithelial": ["TTR", "AQP1", "KRT18", "KRT8", "OTX2", "EPCAM"],
    "Endoderm / hepatic-like": ["AFP", "ALB", "APOA1", "APOA2", "TTR", "FGA"],
    "Pluripotent": ["POU5F1", "NANOG", "DPPA4", "LIN28A"],
    "Neural crest": ["SOX10", "TFAP2A", "FOXD3", "NGFR"],
}
PROGRAMS = {
    "Proliferation": ["MKI67", "TOP2A", "UBE2C", "CENPF", "PCNA", "MCM2"],
    "Forebrain identity": ["FOXG1", "LHX2", "EMX1", "EMX2", "PAX6"],
    "Deep layer cortical": ["TBR1", "BCL11B", "FEZF2", "TLE4"],
    "Upper layer cortical": ["SATB2", "CUX1", "CUX2", "POU3F2"],
    "Synaptic maturation": ["SYN1", "SYP", "SNAP25", "DLG4", "RBFOX3", "MAP2", "CAMK2A",],
    "Stress": ["FOS", "JUN", "HSPA1A", "DDIT3", "ATF4"],
    "Hypoxia": ["VEGFA", "BNIP3", "LDHA", "PGK1"],
}

# Load the count matrix
print("Reading matrix", flush=True)
features = pd.read_csv(DATA / "29.scRNA.filtered.features.tsv.gz", sep="\t", header=None)
barcodes = pd.read_csv(DATA / "29.scRNA.filtered.barcodes.tsv.gz", sep="\t", header=None)

with gzip.open(DATA / "29.scRNA.filtered.matrix.mtx.gz", "rb") as f:
    count_matrix = scipy.io.mmread(f).T.tocsr().astype(np.float32)

adata = ad.AnnData(
    count_matrix,
    obs=pd.DataFrame(index=barcodes[0].values),
    var=pd.DataFrame({"gene_ids": features[0].values, "gene_symbol": features[1].values},index=features[1].values,),
)
adata.var_names_make_unique()
input_shape = list(adata.shape)
print(adata)

# Filter cells
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPL", "RPS")) # Filter ribosomal genes
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, inplace=True)
name = "before"
fig, axs = plt.subplots(1, 3, figsize=(12, 3))
# `total_counts`: Total number of counts for a cell (depth / overall expression).
# `n_genes_by_counts`: The number of genes with at least 1 count in a cell. Calculated for all cells.
# `pct_counts_mt`: Proportion of total counts for a cell which are mitochondrial.
for ax, k in zip(axs, ["total_counts", "n_genes_by_counts", "pct_counts_mt"]):
    ax.hist(adata.obs[k], bins=60)
    ax.set_xlabel(k)
    ax.set_ylabel("Cells")
fig.tight_layout()
fig.savefig(FIG_QC / f"qc_{name}.png")
plt.close(fig)

# High library outliers flagged, not removed.
keep = (
    (adata.obs.n_genes_by_counts >= 300)
    & (adata.obs.total_counts >= 500)
    & (adata.obs.pct_counts_mt <= 20)
)
adata.obs["qc_pass"] = keep
adata.obs.to_csv(OUT_QC / "qc_all_cells.csv")
adata = adata[keep].copy()
name = "after"
fig, axs = plt.subplots(1, 3, figsize=(12, 3))
for ax, k in zip(axs, ["total_counts", "n_genes_by_counts", "pct_counts_mt"]):
    ax.hist(adata.obs[k], bins=60)
    ax.set_xlabel(k)
    ax.set_ylabel("Cells")
fig.tight_layout()
fig.savefig(FIG_QC / f"qc_{name}.png")
plt.close(fig)

print(f"QC retained {adata.n_obs} cells; running doublet detection", flush=True)
# Running doublet detection using Scrublet on raw (unnormalized) counts matrix
sc.pp.scrublet(adata, expected_doublet_rate=0.06, random_state=SEED, verbose=True)
adata.obs.to_csv(OUT_QC / "qc_doublet_flags.csv")
axes = sc.pl.scrublet_score_distribution(adata, show=False)
fig = axes[0].get_figure()
fig.savefig(FIG_QC)

# Filtering genes based on number of cells
# min_cells: Minimum number of cells expressed required for a gene to pass filtering.
sc.pp.filter_genes(adata, min_cells=3)
adata.layers["counts"] = adata.X.copy()

# Normalize expression and select features
# Normalize counts per cell. 
# Normalize each cell by total counts over all genes, so that every cell has the same total count after normalization.
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat")
adata.var.to_csv(OUT_QC / "gene_features_hvg.csv")

# Build the graph and cluster cells
sc.pp.pca(
    adata,
    n_comps=50,
    mask_var="highly_variable",
    svd_solver="arpack",
    random_state=SEED,
)

sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
plt.savefig(os.path.join(FIG_QC, "pca_variance_ratio_sample29.png"))
plt.close("all")

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=SEED)
sc.tl.umap(adata, random_state=SEED)
for r in [0.3, 0.6, 1.0]:
    sc.tl.leiden(
        adata,
        resolution=r,
        key_added=f"leiden_{r}",
        random_state=SEED,
        flavor="igraph",
        n_iterations=2,
    )
adata.obs["cluster"] = adata.obs["leiden_0.6"].copy()
sc.settings.figdir = FIG_QC
sc.pl.umap(
    adata,
    color=["leiden_0.3", "cluster", "leiden_1.0"],
    legend_loc="on data",
    show=False,
    save="_resolutions.png",
)
plt.close("all")
sc.pl.umap(
    adata,
    color=["n_genes_by_counts", "pct_counts_mt", "doublet_score"],
    show=False,
    save="_qc.png",
)
plt.close("all")

# Rank genes and score marker panels
print("Ranking markers and scoring cell identities", flush=True)
# `rank_genes_groups` performs differential gene expression (DGE) analysis to find marker genes that characterize specific groups of cells
sc.tl.rank_genes_groups(adata, "cluster", method="wilcoxon", use_raw=False, pts=True, tie_correct=True)
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv(OUT_QC / "cluster_markers_all.csv", index=False)
markers.groupby("group", observed=True).head(30).to_csv(OUT_QC / "cluster_markers_top30.csv", index=False)

panels = {k: [g for g in v if g in adata.var_names] for k, v in PANELS.items()}
(OUT_MAN_ANNOT / "marker_panels.json").write_text(
    json.dumps({"cell_types": PANELS, "programs": PROGRAMS}, indent=2)
)
for name, genes in {**panels, **PROGRAMS}.items(): # ** is the dictionary unpacking operator, merging the two into a single new dictionary
    genes = [g for g in genes if g in adata.var_names]
    if len(genes) >= 2:
        # Score is the average expression of a set of genes after subtraction by the average expression of a reference set of genes.
        sc.tl.score_genes(adata, genes, score_name=name, random_state=SEED, use_raw=False)

scores = adata.obs.groupby("cluster", observed=True)[list(panels)].mean()
scores.to_csv(OUT_MAN_ANNOT / "cluster_identity_scores.csv")

rows = []

for cluster in scores.index:
    ranked = scores.loc[cluster].sort_values(ascending=False)
    top = ranked.index[0]
    genes = panels[top]
    cluster_data = adata[adata.obs.cluster == cluster]
    fractions = np.asarray((cluster_data[:, genes].X > 0).mean(axis=0)).ravel()
    support = [g for g, f in zip(genes, fractions) if f >= 0.2]
    margin = float(ranked.iloc[0] - ranked.iloc[1])
    label = top if len(support) >= 2 and ranked.iloc[0] > 0.1 else "Unresolved"
    rows.append(
        {
            "cluster": cluster,
            "n_cells": cluster_data.n_obs,
            "suggested_label": label,
            "runner_up": ranked.index[1],
            "score": ranked.iloc[0],
            "margin": margin,
            "supporting_genes": ";".join(support),
            "confidence": "moderate" if margin >= 0.15 and len(support) >= 3 else "low",
            "doublet_fraction": float(cluster_data.obs.predicted_doublet.mean()),
        }
    )
annotations = pd.DataFrame(rows)
annotations.to_csv(OUT_MAN_ANNOT / "cluster_annotation_automatic.csv", index=False)

mapping = dict(zip(annotations.cluster, annotations.suggested_label))
adata.obs["cell_type"] = pd.Categorical(adata.obs.cluster.map(mapping))
adata.obs["cluster_label"] = pd.Categorical([f"{c}: {mapping[c]}" for c in adata.obs.cluster])
print(annotations.to_string(index=False))

# Plot marker expression
sc.settings.figdir = FIG_MAN_ANNOT

sc.pl.umap(adata, color=["cluster", "cell_type"], show=False, save="_annotated.png")
plt.close("all")

sc.pl.dotplot(
    adata,
    panels,
    groupby="cluster",
    standard_scale="var",
    show=False,
    save="_cell_identity.png",
)
plt.close("all")

sc.pl.dotplot(
    adata,
    {k: [g for g in v if g in adata.var_names] for k, v in PROGRAMS.items()},
    groupby="cluster",
    standard_scale="var",
    show=False,
    save="_development.png",
)
plt.close("all")

for name, genes in PROGRAMS.items():
    genes = [g for g in genes if g in adata.var_names]
    sc.pl.umap(
        adata,
        color=genes,
        ncols=3,
        show=False,
        save="_" + name.replace(" ", "_") + ".png",
    )
    plt.close("all")

# Export gene detection and abundance
stats = []
for cluster in adata.obs.cluster.cat.categories:
    cluster_data = adata[adata.obs.cluster == cluster]
    for g in sorted(set(sum(list(PANELS.values()) + list(PROGRAMS.values()), [])) & set(adata.var_names)):
        x = cluster_data[:, g].X
        stats.append(
            {
                "cluster": cluster,
                "gene": g,
                "mean_log1p": float(x.mean()),
                "fraction_expressing": float((x > 0).mean()),
            }
        )
pd.DataFrame(stats).to_csv(OUT_MAN_ANNOT / "marker_expression_by_cluster.csv", index=False)

# Save the analysis and run manifest
adata.obs.to_csv(OUT / "cell_metadata.csv")
pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names, columns=["UMAP1", "UMAP2"]).to_csv(OUT / "umap_coordinates.csv")

# HDF5 disallows slash characters in dataframe column keys.
adata.obs.columns = [c.replace("/", "_") for c in adata.obs.columns]
adata.write_h5ad(OUT / "sample29_annotated.h5ad", compression="gzip")
manifest = {
    "input_shape": input_shape,
    "retained_shape": list(adata.shape),
    "seed": SEED,
    "qc": {"min_genes": 300, "min_counts": 500, "max_mt_percent": 20},
    "doublets_flagged": int(adata.obs.predicted_doublet.sum()),
    "doublets_removed": False,
    "hvg": 3000,
    "pcs": 30,
    "neighbors": 15,
    "resolution": 0.6,
    "python": platform.python_version(),
    "scanpy": sc.__version__,
    "inputs": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.glob("*.gz")},
}
(OUT / "run_manifest.json").write_text(json.dumps(manifest, indent=2))
print("COMPLETE", flush=True)
