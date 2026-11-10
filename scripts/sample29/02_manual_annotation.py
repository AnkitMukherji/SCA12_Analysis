# Import packages
import os

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
import json
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits

threadpool_limits(limits=4)

working_directory = Path.cwd()

ROOT = working_directory
PREPROCESSED_DATA = ROOT / "results/sample29/intermediate_data/sample29_preprocessed.h5ad"
ANNOTATED_DATA = ROOT / "results/sample29/intermediate_data"
OUT = ROOT / "results/sample29/manual_annotation"
FIG = ROOT / "plots/sample29/manual_annotation"

OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

sc.settings.figdir = FIG
sc.settings.n_jobs = 4
sc.set_figure_params(dpi=120, dpi_save=180, facecolor="white")
SEED = 29

adata = ad.read_h5ad(filename=PREPROCESSED_DATA)

# Marker panels and programs
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

panels = {k: [g for g in v if g in adata.var_names] for k, v in PANELS.items()}
(OUT / "marker_panels.json").write_text(
    json.dumps({"cell_types": PANELS, "programs": PROGRAMS}, indent=2)
)
for name, genes in {**panels, **PROGRAMS}.items(): # ** is the dictionary unpacking operator, merging the two into a single new dictionary
    genes = [g for g in genes if g in adata.var_names]
    if len(genes) >= 2:
        # Score is the average expression of a set of genes after subtraction by the average expression of a reference set of genes.
        sc.tl.score_genes(adata, genes, score_name=name, random_state=SEED, use_raw=False)

scores = adata.obs.groupby("cluster", observed=True)[list(panels)].mean()
scores.to_csv(OUT / "cluster_identity_scores.csv")

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
annotations.to_csv(OUT / "cluster_annotation_automatic.csv", index=False)

mapping = dict(zip(annotations.cluster, annotations.suggested_label))
adata.obs["cell_type"] = pd.Categorical(adata.obs.cluster.map(mapping))
adata.obs["cluster_label"] = pd.Categorical([f"{c}: {mapping[c]}" for c in adata.obs.cluster])
print(annotations.to_string(index=False))

# Plot marker expression
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
pd.DataFrame(stats).to_csv(OUT / "marker_expression_by_cluster.csv", index=False)

# Save the analysis
adata.obs.to_csv(OUT / "cell_metadata.csv")

# HDF5 disallows slash characters in dataframe column keys.
adata.obs.columns = [c.replace("/", "_") for c in adata.obs.columns]
adata.write_h5ad(ANNOTATED_DATA / "sample29_annotated.h5ad", compression="gzip")
