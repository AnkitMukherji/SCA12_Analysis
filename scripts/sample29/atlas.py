"""Compare manual labels with developing-brain and pan-fetal references."""

# %% 1. Imports, paths and reference versions
import hashlib
import json
import os
from importlib.metadata import version
from pathlib import Path
from urllib.request import urlretrieve

import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits

matplotlib.use("Agg")

working_directory = Path.cwd()
if (working_directory / "scripts/sample29").is_dir():
    root = working_directory
elif (
    working_directory.name == "sample29" and working_directory.parent.name == "scripts"
):
    root = working_directory.parents[1]
elif "__file__" in globals() and Path(__file__).is_file():
    root = Path(__file__).resolve().parents[2]
else:
    raise RuntimeError("Start Python from the repository root or scripts/sample29.")

results = root / "results/sample29"
output = results / "reference_annotation"
figures = root / "plots/sample29"
references = root / "data/single_cell_analysis/references"
braun_expression = references / "braun2023/HumanFetalBrainPool_cluster_expr_highvar.tsv"
braun_metadata = references / "braun2023/human_dev_science.adf1226_table_s2.xlsx"
fetal_model_path = references / "celltypist/Pan_Fetal_Human.pkl"

# Keep CellTypist's cache inside the project.
os.environ.setdefault("CELLTYPIST_FOLDER", str(fetal_model_path.parent))
from celltypist import annotate, models

reference_urls = {
    braun_expression: (
        "https://raw.githubusercontent.com/theislab/neural_organoid_atlas/"
        "2d4d16c7cfd0c900d64f5a3f56496d2e2c25c5f8/supplemental_files/"
        "HumanFetalBrainPool_cluster_expr_highvar.tsv"
    ),
    braun_metadata: (
        "https://www.dropbox.com/scl/fi/b8zk8124u0hemnolp33w3/"
        "human_dev_science.adf1226_table_s2.xlsx?dl=1&rlkey=li7dop8b8fvmpq0lqnbnvnwri"
    ),
    fetal_model_path: (
        "https://celltypist.cog.sanger.ac.uk/models/Pan_Fetal_Suo/v2/Pan_Fetal_Human.pkl"
    ),
}

# %% 2. Cache the public reference assets and load sample 29
for path, url in reference_urls.items():
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.is_file():
        temporary = path.with_suffix(path.suffix + ".part")
        urlretrieve(url, temporary)
        temporary.replace(path)

output.mkdir(parents=True, exist_ok=True)
figures.mkdir(parents=True, exist_ok=True)
adata = sc.read_h5ad(results / "sample29_annotated.h5ad")
manual = pd.read_csv(
    results / "cluster_annotation_reviewed.csv", dtype={"cluster": str}
)
if "cell_type" not in adata.obs or not manual.cluster.is_unique:
    raise ValueError("Run the manual review stage before reference annotation.")

# Duplicate symbols are ambiguous across atlas features.
unique_symbols = ~adata.var["gene_symbol"].duplicated(keep=False)
query = adata[:, unique_symbols].copy()
query.var_names = query.var["gene_symbol"].astype(str)
query.X = query.layers["counts"].copy()
query.uns.pop("log1p", None)
sc.pp.normalize_total(query, target_sum=10_000)
sc.pp.log1p(query)
clusters = query.obs["cluster"].astype(str)
cluster_ids = sorted(clusters.unique(), key=int)
observed_counts = clusters.value_counts().reindex(manual.cluster)
if not np.array_equal(observed_counts, manual.n_cells):
    raise ValueError("Cluster counts differ from the saved manual review.")

# %% 3. Align the published whole developing brain profiles by cluster ID
brain_profiles = pd.read_csv(braun_expression, sep="\t", index_col=0).T
brain_profiles.index = brain_profiles.index.str.removeprefix("X").astype(int)
brain_labels = pd.read_excel(braun_metadata, sheet_name="Sheet1")
brain_labels["reference_cluster"] = pd.to_numeric(
    brain_labels["ClusterID (PoolClean)"], errors="coerce"
)
brain_labels = brain_labels.dropna(subset=["reference_cluster"]).copy()
brain_labels["reference_cluster"] = brain_labels["reference_cluster"].astype(int)
brain_labels = brain_labels.set_index("reference_cluster")
if not brain_labels.index.is_unique or set(brain_profiles.index) != set(
    brain_labels.index
):
    raise ValueError("The published profile IDs and curated Braun IDs do not match.")
brain_labels = brain_labels.loc[brain_profiles.index]
brain_profiles = brain_profiles.loc[:, ~brain_profiles.columns.duplicated(keep=False)]
brain_genes = brain_profiles.columns.intersection(query.var_names)
if len(brain_genes) < 1_000:
    raise ValueError(f"Only {len(brain_genes)} genes overlap the brain reference.")

# Spearman correlation compares gene ranks across count and platform scales.
reference_ranks = rankdata(brain_profiles[brain_genes], axis=1)
reference_ranks -= reference_ranks.mean(axis=1, keepdims=True)
reference_ranks /= np.maximum(
    np.linalg.norm(reference_ranks, axis=1, keepdims=True), 1e-12
)
print(
    f"Braun reference: {len(brain_profiles)} profiles and {len(brain_genes)} shared genes."
)

# %% 4. Map cluster profiles and measure bootstrap stability
rng = np.random.default_rng(29)
brain_comparison = []
brain_hits = []
brain_classes = brain_labels["Class"].to_numpy()
with threadpool_limits(limits=4):
    for cluster in cluster_ids:
        values = query[clusters.eq(cluster), brain_genes].X.toarray()
        sample_size = len(values)
        means = [values.mean(axis=0)]
        for _ in range(100):
            positions = rng.integers(0, sample_size, size=sample_size)
            means.append(values[positions].mean(axis=0))

        query_ranks = rankdata(np.asarray(means), axis=1)
        query_ranks -= query_ranks.mean(axis=1, keepdims=True)
        query_ranks /= np.maximum(
            np.linalg.norm(query_ranks, axis=1, keepdims=True), 1e-12
        )
        similarities = query_ranks @ reference_ranks.T
        ranked = np.argsort(-similarities[0], kind="stable")
        best = ranked[0]
        brain_class = brain_classes[best]
        runner_up = ranked[brain_classes[ranked] != brain_class][0]
        wins = brain_classes[similarities[1:].argmax(axis=1)]
        best_profile = brain_labels.iloc[best]

        brain_comparison.append(
            {
                "cluster": cluster,
                "brain_class": brain_class,
                "brain_subclass": best_profile["Subclass"],
                "brain_reference_cluster": best_profile.name,
                "brain_reference_region": best_profile["CommonestRegion"],
                "brain_runner_up_class": brain_classes[runner_up],
                "brain_similarity": similarities[0, best],
                "brain_class_margin": similarities[0, best]
                - similarities[0, runner_up],
                "brain_class_bootstrap": np.mean(wins == brain_class),
            }
        )
        for rank, position in enumerate(ranked[:10], start=1):
            hit = brain_labels.iloc[position]
            brain_hits.append(
                {
                    "cluster": cluster,
                    "rank": rank,
                    "reference_cluster": hit.name,
                    "reference_class": hit["Class"],
                    "reference_subclass": hit["Subclass"],
                    "reference_region": hit["CommonestRegion"],
                    "similarity": similarities[0, position],
                }
            )

brain_comparison = pd.DataFrame(brain_comparison)
pd.DataFrame(brain_hits).to_csv(
    output / "brain_top_reference_profiles.csv", index=False
)

# %% 5. Apply the pan-fetal CellTypist classifier to individual cells
# CellTypist models use pickle. Load only the verified public release.
with fetal_model_path.open("rb") as stream:
    model_hash = hashlib.file_digest(stream, "sha256").hexdigest()
if model_hash != "e7f001ff8dbec35ea715eccdbde03f5ce8832743c532ff5a749811ce6b355794":
    raise ValueError(
        f"Unexpected SHA256 for the Pan-Fetal Human v2 model: {model_hash}"
    )
model = models.Model.load(fetal_model_path.as_posix())
model_genes = pd.Index(model.features)
shared_model_genes = model_genes.intersection(query.var_names)
if len(shared_model_genes) < 1_000:
    raise ValueError(f"Only {len(shared_model_genes)} genes overlap CellTypist.")
with threadpool_limits(limits=4):
    fetal_predictions = annotate(query, model=model, majority_voting=False)

decision_scores = fetal_predictions.decision_matrix
scores = decision_scores.to_numpy()
ranked_labels = np.argsort(-scores, axis=1, kind="stable")[:, :2]
positions = np.arange(len(scores))
best = ranked_labels[:, 0]
runner_up = ranked_labels[:, 1]
probabilities = fetal_predictions.probability_matrix.to_numpy()
cell_predictions = query.obs[["cluster", "cell_type", "predicted_doublet"]].copy()
cell_predictions["cluster"] = clusters
cell_predictions["fetal_label"] = decision_scores.columns[best]
cell_predictions["fetal_runner_up"] = decision_scores.columns[runner_up]
cell_predictions["fetal_score"] = probabilities[positions, best]
cell_predictions["fetal_score_margin"] = (
    probabilities[positions, best] - probabilities[positions, runner_up]
)
cell_predictions.to_csv(output / "celltypist_cells.csv", index_label="barcode")

# %% 6. Compare the two references with manual annotations
fetal_rows = []
for cluster in cluster_ids:
    cells = cell_predictions.loc[clusters.eq(cluster)]
    label_counts = cells["fetal_label"].value_counts()
    winning_label = label_counts.index[0]
    winner_cells = cells.loc[cells["fetal_label"].eq(winning_label)]
    fetal_rows.append(
        {
            "cluster": cluster,
            "fetal_top_label": winning_label,
            "fetal_top_fraction": len(winner_cells) / len(cells),
            "fetal_runner_up_label": label_counts.index[1]
            if len(label_counts) > 1
            else None,
            "fetal_median_score_for_top_label": winner_cells["fetal_score"].median(),
            "fetal_median_margin_for_top_label": winner_cells[
                "fetal_score_margin"
            ].median(),
            "fetal_fraction_top_label_and_score_ge_0_5": (
                winner_cells["fetal_score"].ge(0.5).sum() / len(cells)
            ),
        }
    )

comparison = manual[["cluster", "n_cells", "reviewed_label", "confidence"]].merge(
    brain_comparison, on="cluster", validate="one_to_one"
)
comparison = comparison.merge(
    pd.DataFrame(fetal_rows), on="cluster", validate="one_to_one"
)
comparison.to_csv(output / "cluster_reference_comparison.csv", index=False)
fetal_composition = pd.crosstab(cell_predictions.cluster, cell_predictions.fetal_label)
fetal_composition.to_csv(output / "celltypist_cluster_votes.csv")
print(
    comparison[
        ["cluster", "reviewed_label", "brain_class", "fetal_top_label"]
    ].to_string(index=False)
)

# %% 7. Audit manual markers against the top brain profile
marker_evidence = pd.read_csv(
    results / "marker_expression_by_cluster.csv", dtype={"cluster": str}
).merge(
    brain_comparison[["cluster", "brain_reference_cluster", "brain_class"]],
    on="cluster",
    validate="many_to_one",
)
marker_evidence["in_brain_reference"] = marker_evidence["gene"].isin(brain_genes)
marker_evidence["in_celltypist_model"] = marker_evidence["gene"].isin(
    shared_model_genes
)
reference_gene_median = brain_profiles.median()
marker_evidence["brain_profile_expression_minus_reference_median"] = [
    brain_profiles.at[cluster, gene] - reference_gene_median[gene]
    if available
    else np.nan
    for cluster, gene, available in zip(
        marker_evidence["brain_reference_cluster"],
        marker_evidence["gene"],
        marker_evidence["in_brain_reference"],
    )
]
marker_evidence.to_csv(output / "marker_reference_evidence.csv", index=False)

# %% 8. Plot evidence and record input versions
from matplotlib import pyplot as plt

ordered = comparison.set_index("cluster").loc[cluster_ids]
fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
y = np.arange(len(ordered))
axes[0].barh(y, ordered["brain_similarity"], color="#437b9e")
axes[1].barh(y, ordered["fetal_top_fraction"], color="#50896b")
axes[0].set_yticks(y, ordered.index)
axes[0].set_ylabel("Sample 29 cluster")
axes[0].set_xlabel("Braun nearest profile Spearman similarity")
axes[1].set_xlabel("Fraction of cells with the top CellTypist label")
axes[1].set_xlim(0, 1.2)
for row_number, (_, row) in enumerate(ordered.iterrows()):
    axes[0].text(
        row["brain_similarity"] + 0.01,
        row_number,
        row["brain_class"],
        va="center",
        fontsize=8,
    )
    axes[1].text(
        row["fetal_top_fraction"] + 0.01,
        row_number,
        row["fetal_top_label"],
        va="center",
        fontsize=8,
    )
axes[0].set_xlim(0, 1.1)
axes[0].invert_yaxis()
fig.tight_layout()
fig.savefig(
    figures / "reference_annotation_comparison.png", dpi=180, bbox_inches="tight"
)
plt.close(fig)

assets = {}
for path, url in reference_urls.items():
    with path.open("rb") as stream:
        assets[path.name] = {
            "url": url,
            "sha256": hashlib.file_digest(stream, "sha256").hexdigest(),
        }
provenance = {
    "brain_reference": "https://doi.org/10.1126/science.adf1226",
    "fetal_reference": "https://doi.org/10.1126/science.abo0510",
    "brain_reference_profiles": len(brain_profiles),
    "brain_shared_genes": len(brain_genes),
    "celltypist_model_classes": len(model.cell_types),
    "celltypist_shared_genes": len(shared_model_genes),
    "celltypist_version": version("celltypist"),
    "seed": 29,
    "brain_bootstraps": 100,
    "assets": assets,
}
(output / "reference_manifest.json").write_text(
    json.dumps(provenance, indent=2), encoding="utf-8"
)
