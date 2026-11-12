# %% Imports and paths - run from D:/sca12_analysis in environments/sample29_reference
# Sequential workflow; no custom function definitions. Original preprocessing is preserved.
import os
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
from pathlib import Path
ROOT = Path.cwd()
OUT = ROOT / "results/sample29/automated_annotation"
FIG = ROOT / "plots/sample29/automated_annotation"
REF = OUT / "references"
CACHE = OUT / "mapping_cache"
for directory in [OUT, FIG, REF, CACHE]:
    directory.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ["CELLTYPIST_FOLDER"] = str(REF / "celltypist_cache")
import json
import hashlib
import random
import platform
import importlib.metadata
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
ad.settings.allow_write_nullable_strings = True
# scArches 0.6.1 imports the legacy anndata.read alias; current H5AD needs anndata 0.12.
ad.read = ad.read_h5ad
import scanpy as sc
import torch
import celltypist
from celltypist import models
from scarches.models.scpoli import scPoli
from sklearn.neighbors import NearestNeighbors
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
threadpool_limits(limits=4)
torch.set_num_threads(4)
SEED = 29
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
sc.settings.n_jobs = 4
sc.set_figure_params(dpi=120, dpi_save=180, facecolor="white")

# %% Reference provenance: published organoid model plus separately trained primary-brain classifier
# He et al. Nature 2024 https://doi.org/10.1038/s41586-024-08172-8
# Minimal reference: https://doi.org/10.5281/zenodo.15004818 (CC-BY-4.0)
# Model: https://github.com/theislab/neural_organoid_atlas/tree/main/supplemental_files/scpoli_model_params
# Braun et al. Science 2023: https://doi.org/10.1126/science.adf1226
# Official CellTypist model is an separately trained first-trimester-brain classifier, not HNOCA.
INPUT = ROOT / "results/sample29/intermediate_data/sample29_preprocessed.h5ad"
REFERENCE = REF / "hnoca_minimal_for_mapping.h5ad"
MODEL_DIR = REF / "scpoli_HNOCA"
BRAIN_CLASSIFIER = REF / "Developing_Human_Brain.pkl"
assert INPUT.exists() and REFERENCE.exists() and BRAIN_CLASSIFIER.exists()
adata = ad.read_h5ad(INPUT)
original_umap = adata.obsm["X_umap"].copy()
original_clusters = adata.obs["cluster"].astype(str).copy()
assert "counts" in adata.layers and adata.var_names.is_unique
# Reference expression/metadata and feature alignment are inspected before model surgery.
reference = ad.read_h5ad(REFERENCE, backed="r")
model_attributes, reference_state, model_genes = scPoli._load_params(str(MODEL_DIR), map_location="cpu")
print("Reference shape:", reference.shape, flush=True)
print("Reference embeddings:", list(reference.obsm), flush=True)
print("Reference observation columns:", reference.obs.columns.tolist(), flush=True)
print("Model conditions:", model_attributes["condition_keys_"], flush=True)
print("Model label keys:", model_attributes["cell_type_keys_"], flush=True)


# %% Reproducible settings and download provenance
N_EPOCHS = 30
BATCH_SIZE = 256
K_NEIGHBORS = 100
MIN_VOTE = 0.70
MIN_MARGIN = 0.20
REFERENCE_DISTANCE_QUANTILE = 0.99
REFERENCE_CALIBRATION_CELLS = 5000
RETRAIN = False
reference_urls = {
    'hnoca_minimal_for_mapping.h5ad': 'https://zenodo.org/api/records/15004818/files/hnoca_minimal_for_mapping.h5ad/content',
    'scpoli_HNOCA/attr.pkl': 'https://raw.githubusercontent.com/theislab/neural_organoid_atlas/main/supplemental_files/scpoli_model_params/attr.pkl',
    'scpoli_HNOCA/model_params.pt': 'https://raw.githubusercontent.com/theislab/neural_organoid_atlas/main/supplemental_files/scpoli_model_params/model_params.pt',
    'scpoli_HNOCA/var_names.csv': 'https://raw.githubusercontent.com/theislab/neural_organoid_atlas/main/supplemental_files/scpoli_model_params/var_names.csv',
    'Developing_Human_Brain.pkl': 'https://celltypist.cog.sanger.ac.uk/models/Human_Devbrain_Braun/v1/Developing_Human_Brain.pkl',
}
# Files above are downloaded from the authors/official CellTypist repository.
# Download them to REF preserving relative filenames if reusing on another computer.
hashes = {}
for name, path in [('query', INPUT), ('reference', REFERENCE), ('weights', MODEL_DIR/'model_params.pt'), ('attributes', MODEL_DIR/'attr.pkl'), ('model_gene_order', MODEL_DIR/'var_names.csv'), ('classifier', BRAIN_CLASSIFIER)]:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        while True:
            block = handle.read(8*1024*1024)
            if not block:
                break
            digest.update(block)
    hashes[name] = digest.hexdigest()
settings = {'implementation_revision': '2026-11-11-v2', 'package_versions': {name:importlib.metadata.version(name) for name in ['scarches','torch','anndata','scanpy','scikit-learn']}, 'seed': SEED, 'n_epochs': N_EPOCHS, 'batch_size': BATCH_SIZE,
    'k_neighbors': K_NEIGHBORS, 'min_vote': MIN_VOTE, 'min_margin': MIN_MARGIN,
    'reference_distance_quantile': REFERENCE_DISTANCE_QUANTILE,
    'reference_calibration_cells': REFERENCE_CALIBRATION_CELLS, 'sha256': hashes,
    'reference_urls': reference_urls, 'method': 'scPoli frozen-expression query adaptation; Jaccard-square wkNN label transfer'}
cache_key = hashlib.sha256(json.dumps(settings, sort_keys=True).encode()).hexdigest()
(OUT/'run_settings.json').write_text(json.dumps(settings,indent=2),encoding='utf-8')
versions = {name: importlib.metadata.version(name) for name in ['numpy','pandas','scanpy','anndata','scarches','scvi-tools','torch','celltypist','scikit-learn']}
versions['python'] = platform.python_version()
(OUT/'environment_versions.json').write_text(json.dumps(versions,indent=2),encoding='utf-8')
(OUT/'reference_environment_requirements.txt').write_text('\n'.join(name+'=='+importlib.metadata.version(name) for name in sorted({dist.metadata['Name'] for dist in importlib.metadata.distributions() if dist.metadata['Name']})),encoding='utf-8')

# %% Align original counts to reference features; rescue aliases by unique Ensembl IDs
reference_var = reference.var.copy()
reference_var.index = reference.var_names.astype(str)
assert set(model_genes).issubset(reference_var.index)
query_ids = adata.var['gene_ids'].astype(str).str.split('.').str[0]
unique_query_ids = ~query_ids.duplicated(keep=False)
query_id_to_position = pd.Series(np.arange(adata.n_vars)[unique_query_ids], index=query_ids[unique_query_ids])
feature_rows = []
for gene in model_genes:
    position = -1
    match_type = 'missing_zero_padded'
    if gene in adata.var_names:
        position = adata.var_names.get_loc(gene)
        match_type = 'symbol'
    elif 'ensembl' in reference_var:
        ensembl = str(reference_var.loc[gene,'ensembl']).split('.')[0]
        if ensembl in query_id_to_position.index:
            position = int(query_id_to_position[ensembl])
            match_type = 'unique_ensembl_alias'
    feature_rows.append({'reference_gene':gene, 'query_position':position, 'match_type':match_type,
        'query_gene': str(adata.var_names[position]) if position>=0 else ''})
feature_map = pd.DataFrame(feature_rows)
# Never duplicate the same query gene into multiple model features.
duplicate_positions = feature_map.loc[feature_map.query_position>=0,'query_position'].duplicated(keep=False)
assert not duplicate_positions.any(), 'Ambiguous feature mapping requires review.'
feature_map.to_csv(OUT/'hnoca_feature_alignment.csv',index=False)
coverage = float((feature_map.query_position>=0).mean())
assert coverage>=0.80, 'Insufficient overlap with model features; do not force mapping.'
valid = np.flatnonzero(feature_map.query_position.to_numpy()>=0)
selector = sp.csr_matrix((np.ones(len(valid),dtype=np.float32), (feature_map.query_position.to_numpy()[valid],valid)),shape=(adata.n_vars,len(model_genes)))
query_counts = adata.layers['counts'].astype(np.float32) @ selector
query = ad.AnnData(X=query_counts.tocsr(),obs=adata.obs.copy(),var=pd.DataFrame(index=model_genes))
query.obs['batch'] = 'sample29_query'
for key in model_attributes['cell_type_keys_']:
    query.obs[key] = 'unknown'
print('HNOCA model feature coverage:',coverage,flush=True)

# %% Adapt query batch with the original expression network frozen
latent_file = CACHE/'query_latent.npz'
cache_manifest = CACHE/'manifest.json'
reuse = latent_file.exists() and cache_manifest.exists() and not RETRAIN
if reuse:
    reuse = json.loads(cache_manifest.read_text())['cache_key']==cache_key
if reuse:
    query_latent = np.load(latent_file)['latent']
    assert query_latent.shape==(adata.n_obs,10)
else:
    query_model = scPoli.load_query_data(query, str(MODEL_DIR), labeled_indices=[],
        unknown_ct_names=['unknown'], freeze=True, freeze_expression=True,
        remove_dropout=True, map_location='cpu')
    query_model.model.to('cpu')
    before = {name: value.detach().cpu().clone() for name,value in query_model.model.state_dict().items()}
    trainable = [name for name,value in query_model.model.named_parameters() if value.requires_grad]
    (OUT/'trainable_parameters.json').write_text(json.dumps(trainable,indent=2))
    query_model.train(n_epochs=N_EPOCHS,pretraining_epochs=N_EPOCHS,
        batch_size=BATCH_SIZE,lr=0.001,eta=10,weight_decay=0,
        prototype_training=False,unlabeled_prototype_training=False,
        use_early_stopping=False)
    # Old-condition embeddings/dispersion must remain fixed for the stored reference coordinates.
    changed = []
    with torch.no_grad():
        for name,value in query_model.model.state_dict().items():
            if name.startswith('embeddings.'):
                reference_rows = len(model_attributes['conditions_']['batch'])
                value[:reference_rows].copy_(before[name][:reference_rows])
            elif name=='theta':
                value[:,:len(model_attributes['conditions_']['batch'])].copy_(before[name][:,:len(model_attributes['conditions_']['batch'])])
            elif not torch.equal(value.cpu(),before[name]):
                changed.append(name)
    assert not changed, 'Reference expression parameters changed: '+str(changed)
    query_model.save(str(CACHE/'query_scpoli'),overwrite=True)
    # scArches 0.6.1 indexes with torch tensors; dense query avoids SciPy's newer sparse-index restriction.
    query.X = query.X.toarray() if sp.issparse(query.X) else query.X
    query_latent = query_model.get_latent(query,mean=True)
    np.savez_compressed(latent_file,latent=query_latent)
    cache_manifest.write_text(json.dumps({'cache_key':cache_key,'original_expression_weights_unchanged':True},indent=2))
adata.obsm['X_hnoca_scpoli'] = query_latent

# %% Reference-label hierarchy and unchanged reference coordinates
LABEL_KEYS = ['annot_level_1','annot_level_2','annot_level_3_rev2','annot_level_4_rev2']
assert set(LABEL_KEYS).issubset(reference.obs.columns)
reference_labels = reference.obs[LABEL_KEYS].astype(str).copy()
reference_labels.to_csv(OUT/'reference_label_hierarchy_cells.csv.gz',compression='gzip')
reference_labels.drop_duplicates().to_csv(OUT/'reference_label_hierarchy.csv',index=False)
reference_latent = np.asarray(reference.obsm['X_scpoli'],dtype=np.float32)
assert reference_latent.shape[1]==query_latent.shape[1]
# Confirm the downloaded model and stored coordinate system agree on actual reference counts.
check_positions = np.sort(np.random.default_rng(SEED).choice(reference.n_obs,128,replace=False))
check_ref = reference[check_positions,:].to_memory()
check_ref = check_ref[:,model_genes].copy()
check_ref.X = check_ref.layers['counts'].astype(np.float32)
check_ref.X = check_ref.X.toarray() if sp.issparse(check_ref.X) else check_ref.X
for key in model_attributes['cell_type_keys_']:
    check_ref.obs[key] = 'unknown'
check_model = scPoli.load(str(MODEL_DIR),adata=check_ref,map_location='cpu')
check_latent = check_model.get_latent(check_ref,mean=True)
reference_coordinate_error = float(np.max(np.abs(check_latent-reference_latent[check_positions])))
(OUT/'reference_coordinate_check.json').write_text(json.dumps({'n_cells':128,'max_absolute_error':reference_coordinate_error},indent=2))
assert np.allclose(check_latent,reference_latent[check_positions],atol=0.02,rtol=0.02), 'Stored reference coordinates do not match the downloaded model.'
del check_ref,check_model

# %% Exact kNN and published shared-neighbour Jaccard-square edge weighting
# Only reference neighbours reached by the query need their kNN lists computed.
neighbor_file = CACHE/('neighbors_'+cache_key+'.npz')
if reuse and neighbor_file.exists():
    saved_neighbors = np.load(neighbor_file)
    distances = saved_neighbors['distances']
    indices = saved_neighbors['indices']
    edge_weights = saved_neighbors['weights']
    distance_limits = saved_neighbors['distance_limits']
else:
    nn = NearestNeighbors(n_neighbors=K_NEIGHBORS,algorithm='kd_tree',leaf_size=40,n_jobs=4)
    nn.fit(reference_latent)
    distances,indices = nn.kneighbors(query_latent)
    neighbor_union,inverse = np.unique(indices,return_inverse=True)
    print('Reference cells needed for wkNN:',len(neighbor_union),flush=True)
    reference_neighbors = np.empty((len(neighbor_union),K_NEIGHBORS),dtype=np.int32)
    for start in range(0,len(neighbor_union),5000):
        end = min(start+5000,len(neighbor_union))
        reference_neighbors[start:end] = nn.kneighbors(reference_latent[neighbor_union[start:end]],return_distance=False)
    union_positions = inverse.reshape(indices.shape)
    edge_weights = np.zeros(indices.shape,dtype=np.float32)
    for cell in range(adata.n_obs):
        shared = np.isin(reference_neighbors[union_positions[cell]],indices[cell]).sum(axis=1)
        edge_weights[cell] = (shared/(2*K_NEIGHBORS-shared))**2
    # Reference-internal local density: leave self out when calibrating 100th-neighbour distances.
    calibration_positions = np.sort(np.random.default_rng(SEED).choice(reference.n_obs,min(REFERENCE_CALIBRATION_CELLS,reference.n_obs),replace=False))
    calibration_distances = nn.kneighbors(reference_latent[calibration_positions],n_neighbors=K_NEIGHBORS+1)[0][:,-1]
    calibration = pd.DataFrame({'broad_label':reference_labels.iloc[calibration_positions,0].to_numpy(),'distance_k':calibration_distances})
    global_limit = float(np.quantile(calibration_distances,REFERENCE_DISTANCE_QUANTILE))
    limits = calibration.groupby('broad_label').distance_k.quantile(REFERENCE_DISTANCE_QUANTILE)
    n_calibration = calibration.broad_label.value_counts()
    limits.loc[n_calibration[n_calibration<30].index] = global_limit
    nearest_broad = reference_labels.iloc[indices[:,0],0].to_numpy()
    distance_limits = pd.Series(nearest_broad).map(limits).fillna(global_limit).to_numpy()
    calibration.to_csv(OUT/'reference_distance_calibration.csv',index=False)
    limits.to_csv(OUT/'reference_distance_limits.csv')
    np.savez_compressed(neighbor_file,distances=distances.astype(np.float32),indices=indices.astype(np.int32),weights=edge_weights,distance_limits=distance_limits)
    del nn,reference_neighbors
weight_mass = edge_weights.sum(axis=1)
normalized_weights = edge_weights/np.maximum(weight_mass[:,None],1e-12)
distance_ratio = distances[:,-1]/np.maximum(distance_limits,1e-12)
adata.obs['hnoca_neighbor_weight_mass'] = weight_mass
adata.obs['hnoca_reference_distance_ratio'] = distance_ratio
adata.obs['hnoca_distance_flag'] = distance_ratio>1

# %% Transfer complete paths, retain raw labels, and accept only coherent confident branches
# Vote fractions are normalized graph weights, NOT calibrated probabilities.
accepted = np.full(adata.n_obs,'Unassigned',dtype=object)
accepted_depth = np.zeros(adata.n_obs,dtype=int)
accepted_path = np.full(adata.n_obs,'',dtype=object)
score_frames = []
rejection_notes = np.full(adata.n_obs,'',dtype=object)
distance_by_broad = pd.read_csv(OUT/'reference_distance_limits.csv',index_col=0).iloc[:,0]
global_distance_limit = float(pd.read_csv(OUT/'reference_distance_calibration.csv').distance_k.quantile(REFERENCE_DISTANCE_QUANTILE))
for depth,key in enumerate(LABEL_KEYS,start=1):
    paths = reference_labels[key].copy() if depth==1 else paths+' | '+reference_labels[key]
    codes,path_values = pd.factorize(paths,sort=True)
    neighbor_codes = codes[indices]
    votes = np.zeros((adata.n_obs,len(path_values)),dtype=np.float32)
    for cell in range(adata.n_obs):
        votes[cell] = np.bincount(neighbor_codes[cell],weights=normalized_weights[cell],minlength=len(path_values))
    order = np.argsort(votes,axis=1)
    best = order[:,-1]
    runner = order[:,-2] if votes.shape[1]>1 else best
    support = votes[np.arange(adata.n_obs),best]
    margin = support-votes[np.arange(adata.n_obs),runner]
    best_paths = np.asarray(path_values,dtype=object)[best]
    best_labels = pd.Series(best_paths).str.split(' | ',regex=False).str[-1].to_numpy()
    if depth==1:
        # Compare to the density of the voted broad class, not the single nearest cell's class.
        distance_limits = pd.Series(best_labels).map(distance_by_broad).fillna(global_distance_limit).to_numpy()
        distance_ratio = distances[:,-1]/np.maximum(distance_limits,1e-12)
        adata.obs['hnoca_reference_distance_ratio'] = distance_ratio
        adata.obs['hnoca_distance_flag'] = distance_ratio>1
    parent_paths = pd.Series(best_paths).str.rsplit(' | ',n=1).str[0].to_numpy()
    valid_name = ~pd.Series(best_labels).str.lower().isin(['nan','none','unknown','unassigned','']).to_numpy()
    coherent = np.ones(adata.n_obs,dtype=bool) if depth==1 else (parent_paths==accepted_path)&(accepted_depth==depth-1)
    keep = (support>=MIN_VOTE)&(margin>=MIN_MARGIN)&(weight_mass>0)&(distance_ratio<=1)&valid_name&coherent
    keep &= ~adata.obs.predicted_doublet.to_numpy(dtype=bool)
    for cell in range(adata.n_obs):
        if depth==1 or accepted_depth[cell]==depth-1:
            reasons = []
            if support[cell]<MIN_VOTE: reasons.append('low_vote')
            if margin[cell]<MIN_MARGIN: reasons.append('small_margin')
            if weight_mass[cell]<=0: reasons.append('zero_graph_weight')
            if distance_ratio[cell]>1: reasons.append('outside_reference_distance')
            if not valid_name[cell]: reasons.append('missing_reference_label')
            if not coherent[cell]: reasons.append('hierarchy_conflict')
            if adata.obs.predicted_doublet.iloc[cell]: reasons.append('doublet_flag')
            rejection_notes[cell] = ';'.join(reasons)
    accepted[keep] = best_labels[keep]
    accepted_path[keep] = best_paths[keep]
    accepted_depth[keep] = depth
    adata.obs['hnoca_raw_level_'+str(depth)] = pd.Categorical(best_labels)
    adata.obs['hnoca_vote_level_'+str(depth)] = support
    adata.obs['hnoca_margin_level_'+str(depth)] = margin
    votes_table = pd.DataFrame(votes,index=adata.obs_names,columns=path_values)
    votes_table.to_csv(OUT/('hnoca_votes_level_'+str(depth)+'.csv.gz'),compression='gzip')
adata.obs['hnoca_accepted_label'] = pd.Categorical(accepted)
adata.obs['hnoca_accepted_path'] = accepted_path
adata.obs['hnoca_accepted_depth'] = accepted_depth
adata.obs['hnoca_distinct_path_labels'] = [len(set(path.split(' | '))) if path else 0 for path in accepted_path]
adata.obs['hnoca_status'] = pd.Categorical(np.where(accepted_depth==0,'unassigned',np.where(accepted_depth==len(LABEL_KEYS),'deepest_reference_level_supported','broad_or_intermediate_only')))
adata.obs['hnoca_rejection_reasons'] = rejection_notes

# %% Wider-gene classifier: separately trained primary developing brain comparator
# Normalize ALL query genes to 10,000 before selecting classifier features.
classifier_query = ad.AnnData(X=adata.layers['counts'].copy(),obs=adata.obs[['cluster']].copy(),var=adata.var.copy())
sc.pp.normalize_total(classifier_query,target_sum=10000)
sc.pp.log1p(classifier_query)
brain_model = models.Model.load(BRAIN_CLASSIFIER.as_posix())
classifier_shared = np.intersect1d(classifier_query.var_names,brain_model.features)
pd.DataFrame({'model_gene':brain_model.features,'present':np.isin(brain_model.features,classifier_query.var_names)}).to_csv(OUT/'classifier_feature_coverage.csv',index=False)
classifier_result = celltypist.annotate(classifier_query,model=brain_model,majority_voting=False)
classifier_probabilities = classifier_result.probability_matrix
classifier_probabilities.to_csv(OUT/'braun_classifier_probabilities.csv.gz',compression='gzip')
classifier_order = np.argsort(classifier_probabilities.to_numpy(),axis=1)
classifier_best = classifier_order[:,-1]
classifier_second = classifier_order[:,-2]
classifier_max = classifier_probabilities.to_numpy()[np.arange(adata.n_obs),classifier_best]
classifier_margin = classifier_max-classifier_probabilities.to_numpy()[np.arange(adata.n_obs),classifier_second]
adata.obs['braun_raw_label'] = pd.Categorical(classifier_result.predicted_labels['predicted_labels'].reindex(adata.obs_names))
adata.obs['braun_max_score'] = classifier_max
adata.obs['braun_score_margin'] = classifier_margin
# Preserve raw scores rather than claiming that the classifier detects unseen populations.
adata.obs['braun_score_flag'] = (classifier_max<0.50)|(classifier_margin<0.20)
pd.crosstab(adata.obs.hnoca_raw_level_1,adata.obs.braun_raw_label).to_csv(OUT/'hnoca_braun_label_comparison.csv')
del classifier_query

# %% Compare to previous manual annotations without using them for training or forced labels
MANUAL_INPUT = ROOT/'results/sample29/intermediate_data/sample29_manual_annotation.h5ad'
adata.obs['reference_taxonomy_gap_flag'] = False
if MANUAL_INPUT.exists():
    manual = ad.read_h5ad(MANUAL_INPUT,backed='r')
    if 'cell_type' in manual.obs:
        adata.obs['previous_manual_label'] = pd.Categorical(manual.obs['cell_type'].reindex(adata.obs_names))
        # HNOCA EC means ENDOTHELIAL, not non-neural epithelial (paper Fig. 6).
        # These query labels have no explicit counterpart in this reference taxonomy.
        # This is a review flag only; manual labels do not train or overwrite predictions.
        adata.obs['reference_taxonomy_gap_flag'] = adata.obs.previous_manual_label.astype(str).isin(['Non-neural epithelial cells','Endoderm-like cells','Endoderm-like epithelial cells'])
        pd.crosstab(adata.obs.previous_manual_label,adata.obs.hnoca_raw_level_1).to_csv(OUT/'manual_hnoca_comparison.csv')
        pd.crosstab(adata.obs.previous_manual_label,adata.obs.braun_raw_label).to_csv(OUT/'manual_braun_comparison.csv')
    manual.file.close()
adata.obs['reference_review_status'] = pd.Categorical(np.where(adata.obs.reference_taxonomy_gap_flag,'taxonomy_gap_review_required',adata.obs.hnoca_status.astype(str)))
# Withhold a user-facing lineage assertion where the reference lacks a suspected query lineage.
# Raw/method-accepted predictions remain unchanged above; this is not manual relabelling.
adata.obs['automated_review_label'] = pd.Categorical(np.where(adata.obs.reference_taxonomy_gap_flag,'Review required: reference coverage',adata.obs.hnoca_accepted_label.astype(str)))

# %% Cluster summaries preserve within-cluster heterogeneity
summary_rows = []
for cluster in sorted(adata.obs.cluster.astype(str).unique(),key=int):
    group = adata.obs.loc[adata.obs.cluster.astype(str)==cluster]
    for column in ['hnoca_raw_level_1','hnoca_raw_level_2','hnoca_raw_level_3','hnoca_raw_level_4','hnoca_accepted_label','automated_review_label','braun_raw_label']:
        counts = group[column].value_counts()
        counts = counts[counts>0]
        for label,count in counts.items():
            summary_rows.append({'cluster':cluster,'annotation':column,'label':label,'n_cells':int(count),'fraction_of_cluster':float(count/len(group))})
summary = pd.DataFrame(summary_rows)
summary.to_csv(OUT/'cluster_label_distributions.csv',index=False)
cluster_overview = adata.obs.groupby('cluster',observed=True).agg(n_cells=('hnoca_accepted_depth','size'),mean_accepted_depth=('hnoca_accepted_depth','mean'),distance_flag_fraction=('hnoca_distance_flag','mean'),doublet_fraction=('predicted_doublet','mean'),classifier_score_flag_fraction=('braun_score_flag','mean'))
cluster_overview.to_csv(OUT/'cluster_uncertainty_summary.csv')
adata.obs.to_csv(OUT/'cell_annotations.csv')

# %% Plot original UMAP with reference labels and uncertainty; no UMAP recalculation
plot_specs = {
    'umap_original_clusters':['cluster'],
    'umap_automated_annotation':['automated_review_label'],
    'umap_hnoca_hierarchy':['hnoca_raw_level_1','hnoca_raw_level_2','hnoca_raw_level_3','hnoca_raw_level_4'],
    'umap_hnoca_accepted':['hnoca_accepted_label','hnoca_status'],
    'umap_reference_review':['reference_review_status','reference_taxonomy_gap_flag'],
    'umap_mapping_uncertainty':['hnoca_vote_level_1','hnoca_vote_level_4','hnoca_reference_distance_ratio','hnoca_accepted_depth'],
    'umap_braun_classifier':['braun_raw_label','braun_max_score','braun_score_margin','braun_score_flag'],
}
for name,columns in plot_specs.items():
    sc.pl.umap(adata,color=columns,ncols=2,show=False,wspace=0.7)
    plt.savefig(FIG/(name+'.png'),dpi=180,bbox_inches='tight')
    plt.savefig(FIG/(name+'.pdf'),bbox_inches='tight')
    plt.close('all')
composition = pd.crosstab(adata.obs.cluster,adata.obs.automated_review_label,normalize='index')
composition.plot.bar(stacked=True,figsize=(14,7),width=0.9)
plt.ylabel('Fraction of cells in cluster')
plt.legend(bbox_to_anchor=(1.02,1),loc='upper left',fontsize=7)
plt.savefig(FIG/'cluster_accepted_composition.png',dpi=180,bbox_inches='tight')
plt.savefig(FIG/'cluster_accepted_composition.pdf',bbox_inches='tight')
plt.close('all')

# Post-hoc biological check only: these genes do not train or select transferred labels.
diagnostic_genes = ['SOX2','PAX6','NES','FABP7','DCX','STMN2','FOXG1','EMX2','HOXB9','PAX3','GAD1','SLC17A6','TTR','CLIC6','PRRX1','COL1A1','DCN','EPCAM','CDH1','KRT13','CDX2','HNF4A','VIL1','PECAM1']
diagnostic_genes = [gene for gene in diagnostic_genes if gene in adata.var_names]
diagnostic_plot = sc.pl.dotplot(adata,diagnostic_genes,groupby='hnoca_accepted_label',use_raw=False,standard_scale='var',return_fig=True,show=False)
diagnostic_plot.savefig(FIG/'accepted_labels_marker_check.png',dpi=180,bbox_inches='tight')
diagnostic_plot.savefig(FIG/'accepted_labels_marker_check.pdf',bbox_inches='tight')
plt.close('all')

# %% Save annotated intermediate and methods
assert np.array_equal(original_umap,adata.obsm['X_umap'])
assert original_clusters.equals(adata.obs.cluster.astype(str))
adata.uns['automated_annotation'] = {'method':settings['method'],'reference':'HNOCA minimal Zenodo15004818','parameters_json':json.dumps(settings),'feature_coverage':coverage,'classifier_feature_coverage':len(classifier_shared)/len(brain_model.features),'confidence':'heuristic votes/distance; not calibrated accuracy'}
adata.write_h5ad(OUT/'sample29_automated_annotation.h5ad',compression='gzip')
(OUT/'README.md').write_text("""# Reference-based annotation of sample29

Run from D:/sca12_analysis:

    environments/sample29_reference/Scripts/python.exe -u scripts/sample29/03_automated_annotation.py

Use the Python 3.11 interpreter in environments/sample29_reference when executing sections interactively. The existing sample29 preprocessing environment is unchanged. The code uses sequential sections with no custom functions. reference_environment_requirements.txt and environment_versions.json record the installed environment. The anndata.read alias restores compatibility with scArches 0.6.1 while reading current H5AD encoding.

## Reference selection

HNOCA is the primary reference because it contains neural organoids across protocols and ages, including unguided organoids. The published minimal reference has 1,770,578 cells and 3,000 features. The saved authors' scPoli model is adapted to sample29 using raw counts and an unlabelled new batch, with the expression network frozen. Manual annotations are never training labels. Official URLs and SHA256 checksums are in run_settings.json; files are stored in references/. The HNOCA minimal file's published MD5 is 078675d6108e93cebc99676b6b0626aa. See reference_selection.md for the research comparison.

The separately trained Braun first-trimester human brain CellTypist model supplies a wider-gene classifier comparison. It is not an independent biological validation: HNOCA itself used developing-brain information, including Braun. Its primary-brain label vocabulary can misclassify organoid off-target cells. Do not combine its sigmoid scores numerically with HNOCA vote fractions or accept its highest score as proof that a cell belongs to a represented class.

Regional labels are reference-dependent hypotheses. Disagreement between HNOCA and Braun is retained rather than resolved by an unvalidated ensemble. A high vote for a regional NPC or neuron label does not establish anatomical origin. NPC means neural progenitor cell; IP intermediate progenitor; OPC oligodendrocyte progenitor; PSC pluripotent stem cell. NT-VGLUT denotes a glutamatergic-transporter-associated reference label, not functional neurotransmission measured in this sample.

The downloaded HNOCA taxonomy also lacks explicit non-neural epithelial and endoderm classes. EC means endothelial cell, MC mesenchymal cell, CP choroid plexus, and NC Derivatives neural crest derivatives. A reference_taxonomy_gap_flag marks cells previously assigned epithelial/endoderm identities by the query DE workflow. This is a review warning only, not a training label or an override of model predictions. HNOCA-supported nearest labels for those cells must not be treated as a validated replacement. The user-facing automated_review_label withholds a lineage assertion for these cells as 'Review required: reference coverage'; hnoca_accepted_label and all raw outputs remain unchanged for inspection. Manual labels themselves are uncertain, and the gap flag is not exhaustive because manual-unresolved cells can also be out of reference. hnoca_distinct_path_labels distinguishes genuinely more specific paths from labels repeated unchanged across four levels.

## Mapping and uncertainty

Reference feature alignment first uses gene symbols and then unambiguous Ensembl identifiers; absent features are zero-padded and reported. HNOCA uses original counts. CellTypist separately uses log1p(counts per 10,000) computed across all genes before model-feature selection.

The query is adapted for 30 epochs on CPU. Old reference embedding/dispersion rows are preserved and expression weights checked. A reference subset is re-encoded to check that the downloaded model matches the stored reference latent coordinates. Exact nearest neighbours use k=100. The published HNOCA graph formula is reproduced: weight = (shared neighbours/(2*k - shared neighbours))^2, retaining only query-to-reference nearest-neighbour edges. Reference neighbour lists are computed only for the union reached by the query. Explicit reference kneighbors includes self, matching this graph convention.

Transferred labels are complete reference hierarchy paths, not unrelated maxima independently combined across levels. A deeper label requires consistency with the accepted parent, at least 70% normalized weighted vote and a margin of at least 20 percentage points. Doublet-flagged cells and cells outside the reference-distance threshold remain unassigned. Otherwise the result falls back to its last supported parent. The deepest reference level is not necessarily a fine biological subtype; several labels repeat across levels.

Distance rejection compares the query's 100th-neighbour distance with the 99th percentile of reference leave-self-out distances for its voted broad class. Calibration uses a reproducible 5,000-cell reference sample; classes with fewer than 30 calibration cells use the overall threshold. This is a heuristic novelty flag, not a validated out-of-distribution detector. Integration may hide biological differences; high votes or classifier scores can occur for populations absent from the reference. Manual DE evidence and conflicting reference predictions still require review.

## Outputs to inspect

- cell_annotations.csv: raw labels at all levels, accepted label/path/depth, vote fractions, margins, reference-distance flags, rejection reasons, doublet flags and raw Braun predictions/scores.
- cluster_label_distributions.csv: all per-cell label proportions within original clusters, preserving mixed clusters rather than forcing one label onto every member.
- cluster_uncertainty_summary.csv: cluster-level uncertainty summaries.
- hnoca_votes_level_*.csv.gz: normalized graph votes for all reference paths; these are not calibrated probabilities.
- hnoca_feature_alignment.csv and classifier_feature_coverage.csv: feature coverage and missing features.
- reference_label_hierarchy.csv: actual reference hierarchy paths.
- reference_coordinate_check.json and run_settings.json: coordinate validation and reproducibility.
- sample29_automated_annotation.h5ad: annotated query with unchanged original counts, expression, clustering and UMAP, plus a separate mapped latent representation.

Plots are in plots/sample29/automated_annotation. All label/uncertainty plots use the original UMAP. A reference-space UMAP is not recomputed.

The mapping cache is checked against input/reference/model hashes, gene order, key package versions, settings and implementation revision. Set RETRAIN=True to force fresh mapping. Changing thresholds currently invalidates the cache conservatively. Original data and manual-annotation outputs are not overwritten. No query data are uploaded to an external service.

Sources: https://doi.org/10.1038/s41586-024-08172-8 ; https://doi.org/10.5281/zenodo.15004818 ; https://devsystemslab.github.io/HNOCA-tools/vignettes/atlas_mapping/ ; https://doi.org/10.1126/science.adf1226 ; https://www.sc-best-practices.org/cellular-structure/annotation/ . Acceptance thresholds and novelty rejection are explicit conservative additions, not published calibrated error guarantees.
""",encoding='utf-8')
reference.file.close()
print('Completed automated annotation:',OUT,flush=True)
