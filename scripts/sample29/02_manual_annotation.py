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
original_umap = adata.obsm["X_umap"].copy()
adata.obs["cluster"] = adata.obs["cluster"].astype(str)
cluster_order = sorted(adata.obs["cluster"].unique(), key=int)
adata.obs["cluster"] = pd.Categorical(adata.obs["cluster"], categories=cluster_order)

# Broad marker panels across developmental stages
PANELS = {
    "Neural progenitors": ["SOX2", "SOX1", "PAX6", "NES", "HES5", "FABP7", "PTPRZ1", "HES1"],
    "Neurons": ["DCX", "STMN2", "TUBB3", "NEUROD1", "ELAVL3", "RBFOX3", "SNAP25", "SYT1"],
    "Astroglial cells": ["GFAP", "AQP4", "ALDH1L1", "SLC1A3", "S100B"],
    "Oligodendroglial lineage": ["OLIG1", "OLIG2", "SOX10", "PDGFRA", "PLP1", "MBP", "MOG"],
    "Microglia / myeloid cells": ["AIF1", "TYROBP", "CSF1R", "C1QA", "P2RY12", "PTPRC"],
    "Mesenchymal / stromal cells": ["PRRX1", "TWIST1", "COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB", "TAGLN"],
    "Endothelial cells": ["PECAM1", "CDH5", "VWF", "EMCN", "CLDN5", "KDR", "ROBO4"],
    "Choroid plexus-like cells": ["TTR", "OTX2", "CLIC6", "KCNJ13", "SLC13A4", "AQP1"],
    "Non-neural epithelial cells": ["EPCAM", "CDH1", "KRT8", "KRT18", "TFAP2A", "GATA3", "KRT4", "KRT13"],
    "Endoderm-like cells": ["CDX2", "HNF4A", "CDH17", "VIL1", "AFP", "FABP1", "APOA2"],
    "Neural crest-like cells": ["SOX10", "FOXD3", "NGFR", "TFAP2A"],
    "Pluripotent-like cells": ["POU5F1", "NANOG", "DPPA4", "LIN28A"],
}
# Programs overlap stages and lineages.
PROGRAMS = {
    "Proliferation": ["MKI67", "TOP2A", "UBE2C", "PCNA", "MCM2"],
    "Neurogenic transition": ["EOMES", "NEUROG2", "INSM1", "HES6"],
    "Neuronal maturation": ["MAP2", "RBFOX3", "SYN1", "SYP", "SNAP25", "DLG4"],
    "Glutamatergic program": ["SLC17A6", "SLC17A7", "NEUROD2", "NEUROD6"],
    "GABAergic program": ["GAD1", "GAD2", "SLC32A1", "DLX1", "DLX2"],
    "Forebrain patterning": ["FOXG1", "EMX1", "EMX2", "LHX2"],
    "Posterior patterning": ["HOXB9", "HOXB8", "HOXA5", "PAX3"],
    "Dorsal midline patterning": ["LMX1A", "WNT3A", "WNT8B", "RSPO2"],
    "Heat shock / stress": ["HSPA6", "HSPA1A", "HSPH1", "DDIT3", "ATF4"],
    "Hypoxia / glycolysis": ["VEGFA", "BNIP3", "LDHA", "PGK1"],
}

marker_panel = {**PANELS, **PROGRAMS}
(OUT / "marker_panel.json").write_text(json.dumps(marker_panel, indent=2), encoding='utf-8')

TOP_N = 200
RANK_WINDOWS = [100, 200, 300]
MAX_ADJUSTED_P = 0.05
MIN_LOG2FC = 1.0
MIN_DETECTION = 0.20 # proportion of cells expressing a specific gene within the target cluster versus the comparison group
MIN_MATCHES = 3 # minimum number of genes from the panel that should match
MIN_PANEL_FRACTION = 0.30 # Proportion of genes in the panel present in a cluster / Proportion of genes in the panel available in the entire gene list
MIN_PANEL_AVAILABILITY = 0.70 # Proportion of genes in the panel available in the entire gene list
MIN_RELATIVE_MARGIN = 0.25
MIN_CLUSTER_CELLS = 50
MAX_DOUBLET_FRACTION = 0.10

DE_KEY = "rank_genes_groups"
params = adata.uns[DE_KEY]["params"]

de = sc.get.rank_genes_groups_df(adata, group=None, key=DE_KEY)
de = de.rename(columns={'group': 'cluster', 'names': 'gene'})
de['cluster'] = de.cluster.astype(str)
de['rank'] = de.groupby('cluster', sort=False).cumcount() + 1
de['eligible_positive_de'] = (de.pvals_adj < MAX_ADJUSTED_P) & (de.logfoldchanges > MIN_LOG2FC) & (de.pct_nz_group > MIN_DETECTION)
de['used_for_annotation'] = de.eligible_positive_de & (de['rank'] <= TOP_N)
de['cluster_specific'] = (de.logfoldchanges > MIN_LOG2FC) & (de.pct_nz_group > MIN_DETECTION) & (de.pct_nz_reference < 0.20)
de.to_csv(OUT / 'cluster_de_all.csv.gz', index=False, compression='gzip')
de.loc[de['rank'] <= TOP_N].to_csv(OUT / 'cluster_de_top200.csv', index=False)
de.loc[de.used_for_annotation].to_csv(OUT / 'cluster_de_annotation_input.csv', index=False)
de.loc[de.cluster_specific].to_csv(OUT / 'cluster_de_specific.csv', index=False)

reference_rows = []
for kind, collection in [('identity', PANELS), ('program', PROGRAMS)]:
    for panel, genes in collection.items():
        for gene in dict.fromkeys(genes):
            reference_rows.append({'kind': kind, 'panel': panel, 'gene': gene, 'available': gene in adata.var_names})
reference = pd.DataFrame(reference_rows)
reference.to_csv(OUT / 'marker_reference.csv', index=False)
evidence = de.merge(reference, on='gene', how='inner')
evidence.to_csv(OUT / 'marker_de_evidence.csv', index=False)
evidence.loc[evidence.used_for_annotation].to_csv(OUT / 'matched_marker_evidence.csv', index=False)

# Count overlaps for EVERY cluster and panel
# Score = matched eligible top-DE genes / panel genes present in the data.
# Both the count and fraction are required, so one marker cannot assign a label.
match_rows = []
for window in RANK_WINDOWS:
    for cluster in cluster_order:
        top = de.loc[(de.cluster == cluster) & de.eligible_positive_de & (de['rank'] <= window)]
        for kind, collection in [('identity', PANELS), ('program', PROGRAMS)]:
            for panel, genes in collection.items():
                available = [g for g in dict.fromkeys(genes) if g in adata.var_names]
                matches = top.loc[top.gene.isin(available)]
                fraction = len(matches) / len(available) if available else 0.0
                match_rows.append({'cluster': cluster, 'top_n': window, 'kind': kind,
                    'panel': panel, 'n_panel_genes': len(genes), 'n_available': len(available),
                    'panel_availability': len(available) / len(genes), 'n_matches': len(matches),
                    'panel_fraction': fraction, 'n_eligible_top_de': len(top),
                    'matched_genes': ';'.join(matches.gene),
                    'matched_original_ranks': ';'.join(matches['rank'].astype(str))})
matches = pd.DataFrame(match_rows)
matches.to_csv(OUT / 'panel_overlap_all_windows.csv', index=False)
matches.loc[matches.top_n == TOP_N].to_csv(OUT / 'panel_overlap.csv', index=False)

# Candidate selection and uncertainty, including sensitivity to rank cutoff
candidate_rows = []
for window in RANK_WINDOWS:
    for cluster in cluster_order:
        ranked = matches.loc[(matches.cluster == cluster) & (matches.top_n == window) & (matches.kind == 'identity') & (matches.panel_availability >= MIN_PANEL_AVAILABILITY)]
        ranked = ranked.sort_values(['panel_fraction', 'n_matches', 'panel'], ascending=[False, False, True])
        assert len(ranked) >= 2
        best = ranked.iloc[0]
        second = ranked.iloc[1]
        margin = (best.panel_fraction - second.panel_fraction) / best.panel_fraction if best.panel_fraction > 0 else 0.0
        reasons = []
        if best.n_matches < MIN_MATCHES:
            reasons.append('too_few_matching_genes')
        if best.panel_fraction < MIN_PANEL_FRACTION:
            reasons.append('insufficient_panel_coverage')
        if margin < MIN_RELATIVE_MARGIN:
            reasons.append('competing_panels_or_tie')
        candidate_rows.append({'cluster': cluster, 'top_n': window,
            'best_candidate': best.panel if best.n_matches > 0 else 'None',
            'runner_up': second.panel if second.n_matches > 0 else 'None',
            'n_matches': int(best.n_matches), 'n_available': int(best.n_available),
            'matched_genes': best.matched_genes, 'panel_fraction': best.panel_fraction,
            'runner_up_fraction': second.panel_fraction, 'relative_margin': margin,
            'evidence_index': best.panel_fraction * margin,
            'window_label': best.panel if not reasons else 'Unresolved',
            'overlap_review_reason': ';'.join(reasons)})
candidates = pd.DataFrame(candidate_rows)
candidates.to_csv(OUT / 'rank_cutoff_sensitivity.csv', index=False)

annotations = candidates.loc[candidates.top_n == TOP_N].set_index('cluster').copy()
annotations['n_cells'] = adata.obs.groupby('cluster', observed=True).size()
annotations['doublet_fraction'] = adata.obs.groupby('cluster', observed=True).predicted_doublet.mean()
annotations['stable_across_rank_windows'] = False
annotations['cell_type'] = 'Unresolved'
annotations['confidence'] = 'low'
annotations['review_reason'] = ''

for cluster in cluster_order:
    row = annotations.loc[cluster]
    windows = candidates.loc[candidates.cluster == cluster]
    stable = (windows.window_label == row.window_label).all()
    annotations.loc[cluster, 'stable_across_rank_windows'] = stable
    reasons = [r for r in row.overlap_review_reason.split(';') if r]
    if not stable:
        reasons.append('sensitive_to_top_gene_cutoff')
    if row.n_cells < MIN_CLUSTER_CELLS:
        reasons.append('small_cluster')
    if row.doublet_fraction > MAX_DOUBLET_FRACTION:
        reasons.append('doublet_enrichment')
    if not reasons:
        annotations.loc[cluster, 'cell_type'] = row.window_label
        annotations.loc[cluster, 'confidence'] = 'moderate'
    annotations.loc[cluster, 'review_reason'] = ';'.join(reasons) or 'coherent_overlap_requires_manual_review'
annotations.to_csv(OUT / 'cluster_annotations.csv')

# Propagate cluster inference to member cells
adata.obs['cell_type'] = pd.Categorical(adata.obs.cluster.astype(str).map(annotations.cell_type))
adata.obs['annotation_confidence'] = pd.Categorical(adata.obs.cluster.astype(str).map(annotations.confidence))
adata.obs['annotation_evidence_index'] = adata.obs.cluster.astype(str).map(annotations.evidence_index).astype(float)
adata.obs['annotation_review_reason'] = adata.obs.cluster.astype(str).map(annotations.review_reason)
composition = adata.obs.cell_type.value_counts().rename('n_cells').to_frame()
composition['percent_cells'] = 100 * composition.n_cells / adata.n_obs
composition.to_csv(OUT / 'cell_type_composition.csv')
adata.obs.to_csv(OUT / 'cell_metadata.csv')

# DE and marker plots
sc.tl.dendrogram(adata, groupby='cluster', use_rep='X_pca')
plot = sc.pl.rank_genes_groups_dotplot(adata, groupby='cluster', key=DE_KEY, n_genes=5, standard_scale='var', return_fig=True, show=False)
plot.savefig(FIG / 'de_top5_dotplot.png', dpi=180, bbox_inches='tight')
plt.close('all')

specific_top = {}
for cluster in cluster_order:
    genes = de.loc[(de.cluster == cluster) & de.cluster_specific].head(5).gene.tolist()
    if genes:
        specific_top[cluster] = genes

for filename, collection in [('de_specific_top5', specific_top), ('identity_markers', PANELS), ('developmental_programs', PROGRAMS)]:
    available = {name: [g for g in genes if g in adata.var_names] for name, genes in collection.items()}
    available = {name: genes for name, genes in available.items() if genes}
    plot = sc.pl.dotplot(adata, available, groupby='cluster', use_raw=False, standard_scale='var', return_fig=True, show=False)
    plot.savefig(FIG / (filename + '_dotplot.png'), dpi=180, bbox_inches='tight')
    plt.close('all')

for filename, color in [('umap_clusters', 'cluster'), ('umap_cell_types', 'cell_type'), ('umap_confidence', 'annotation_confidence')]:
    sc.pl.umap(adata, color=color, show=False, title=filename.replace('_', ' ') + ' (cluster-level)')
    plt.savefig(FIG / (filename + '.png'), dpi=180, bbox_inches='tight')
    plt.close('all')

# Overlap heatmap exposes competing identities, not only the winning panel.
overlap = matches.loc[(matches.top_n == TOP_N) & (matches.kind == 'identity')].pivot(index='cluster', columns='panel', values='panel_fraction').reindex(cluster_order)
fig, ax = plt.subplots(figsize=(12, 7))
im = ax.imshow(overlap.values, vmin=0, vmax=1, cmap='Blues', aspect='auto')
ax.set_xticks(range(len(overlap.columns)), overlap.columns, rotation=70, ha='right')
ax.set_yticks(range(len(overlap.index)), overlap.index)
ax.set_ylabel('Cluster')
ax.set_title('Fraction of panel markers matching eligible top ' + str(TOP_N) + ' DE genes')
fig.colorbar(im, ax=ax, label='Panel overlap fraction (not probability)')
fig.savefig(FIG / 'panel_overlap.png', dpi=180, bbox_inches='tight')
plt.close('all')

# Save method and explain all rules
settings = {'de_params': params, 'top_n': TOP_N, 'rank_windows': RANK_WINDOWS,
    'max_adjusted_p': MAX_ADJUSTED_P, 'min_log2fc': MIN_LOG2FC, 'min_detection': MIN_DETECTION,
    'min_matches': MIN_MATCHES, 'min_panel_fraction': MIN_PANEL_FRACTION,
    'min_panel_availability': MIN_PANEL_AVAILABILITY, 'min_relative_margin': MIN_RELATIVE_MARGIN,
    'min_cluster_cells': MIN_CLUSTER_CELLS, 'max_doublet_fraction': MAX_DOUBLET_FRACTION,
    'require_stability_across_windows': True}
(OUT / 'annotation_settings.json').write_text(json.dumps(settings, indent=2), encoding='utf-8')

adata.uns['annotation_method'] = {'method': 'ranked_DE_panel_overlap', 'status': 'rule_based_suggestions_for_manual_review', 'confidence': 'heuristic, not probability', 'settings_json': json.dumps(settings)}
adata.write_h5ad(ANNOTATED_DATA / 'sample29_manual_annotation.h5ad', compression='gzip')
