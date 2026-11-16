# Import packages
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
PREPROCESSED_DATA = ROOT / "results/sample29/intermediate_data"
OUT = ROOT / "results/sample29/preprocessing"
FIG = ROOT / "plots/sample29/preprocessing"

PREPROCESSED_DATA.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

sc.settings.figdir = FIG
sc.settings.n_jobs = 4
sc.set_figure_params(dpi=120, dpi_save=180, facecolor="white")
SEED = 29

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
fig.savefig(FIG / f"qc_{name}.png")
plt.close(fig)

# High library outliers flagged, not removed.
keep = (
    (adata.obs.n_genes_by_counts >= 300)
    & (adata.obs.total_counts >= 500)
    & (adata.obs.pct_counts_mt <= 20)
)
adata.obs["qc_pass"] = keep
adata.obs.to_csv(OUT / "qc_all_cells.csv")
adata = adata[keep].copy()
name = "after"
fig, axs = plt.subplots(1, 3, figsize=(12, 3))
for ax, k in zip(axs, ["total_counts", "n_genes_by_counts", "pct_counts_mt"]):
    ax.hist(adata.obs[k], bins=60)
    ax.set_xlabel(k)
    ax.set_ylabel("Cells")
fig.tight_layout()
fig.savefig(FIG / f"qc_{name}.png")
plt.close(fig)

print(f"QC retained {adata.n_obs} cells; running doublet detection", flush=True)
# Running doublet detection using Scrublet on raw (unnormalized) counts matrix
sc.pp.scrublet(adata, expected_doublet_rate=0.06, random_state=SEED, verbose=True)
adata.obs.to_csv(OUT / "qc_doublet_flags.csv")
sc.pl.scrublet_score_distribution(adata, show=False, save="_sample29.png")

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
adata.var.to_csv(OUT / "gene_features_hvg.csv")

# Build the graph and cluster cells
sc.pp.pca(
    adata,
    n_comps=50,
    mask_var="highly_variable",
    svd_solver="arpack",
    random_state=SEED,
)

sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
plt.savefig(os.path.join(FIG, "pca_variance_ratio_sample29.png"))
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
markers.to_csv(OUT / "cluster_markers_all.csv", index=False)
markers.groupby("group", observed=True).head(30).to_csv(OUT / "cluster_markers_top30.csv", index=False)

# Save the analysis and run manifest
adata.obs.to_csv(OUT / "cell_metadata.csv")
pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names, columns=["UMAP1", "UMAP2"]).to_csv(OUT / "umap_coordinates.csv")

# HDF5 disallows slash characters in dataframe column keys.
adata.obs.columns = [c.replace("/", "_") for c in adata.obs.columns]
adata.write_h5ad(PREPROCESSED_DATA / "sample29_preprocessed.h5ad", compression="gzip")
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
