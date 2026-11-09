from pathlib import Path
import csv
from docx import Document
from docx.shared import Inches, Pt, RGBColor

R = Path(__file__).resolve().parents[2]
O = R / "results/sample29"
F = R / "plots/sample29"
d = Document()
sec = d.sections[0]
sec.top_margin = sec.bottom_margin = Inches(0.7)
style = d.styles["Normal"]
style.font.name = "Calibri"
style.font.size = Pt(10)
for n in ["Title", "Heading 1", "Heading 2"]:
    d.styles[n].font.color.rgb = RGBColor(0, 0, 0)


d.add_heading("Sample 29 single cell RNA sequencing report", 0)
d.add_paragraph(
    "Neuronal organoid composition and developmental interpretation | 28 October 2026"
)
d.add_paragraph(
    "The sample contains neural progenitors and immature neurons together with substantial epithelial and mesenchymal populations. The analysis supports heterogeneous differentiation, with regional variation among neural progenitors. It does not establish uniformly specified or functionally mature neuronal development."
)
d.add_heading("Data and processing", 1)
d.add_paragraph(
    "The filtered matrix contained 8,721 cells and 40,407 features. Quality control retained 8,611 cells (98.74%) and 32,091 genes detected in at least three retained cells. All 110 excluded cells exceeded 20% mitochondrial counts. The supplied vendor UMAP and markers were not used to determine the new clustering."
)
d.add_paragraph(
    "QC thresholds were at least 300 genes, at least 500 counts, and at most 20% mitochondrial counts. Counts were normalized to 10,000 per cell and log1p transformed. The workflow selected 3,000 highly variable genes, calculated 50 principal components, used 30 components and 15 neighbors for the graph, and calculated UMAP and Leiden clustering. Primary Leiden resolution was 0.6; resolutions 0.3 and 1.0 provide sensitivity comparisons. Seed 29 was used. No batch correction or cell cycle regression was applied."
)
d.add_paragraph(
    "Scrublet flagged 151 possible doublets (1.75%) using an assumed 6% expected rate. Flags were retained in the primary analysis because sample loading metadata was unavailable. The counts layer preserves untransformed counts; X stores normalized log expression."
)
d.add_heading("Principal findings", 1)
d.add_paragraph(
    "Cluster 11 comprises 553 immature neurons (6.42%). Neural progenitor related states occupy clusters 5, 7, 8, 10 and 12. Epithelial, squamous, mesenchymal and hepatointestinal-like programs indicate substantial differentiation outside the intended neuronal population. Cluster 4 has a more specific choroid plexus-like epithelial signature."
)
d.add_paragraph(
    "These assignments are marker-based interpretations, not reference-atlas predictions. Moderate and low confidence categories are qualitative assessments rather than probabilities."
)
d.add_page_break()
d.add_heading("Cellular structure", 1)
d.add_picture(str(F / "umap_reviewed_cell_types.png"), width=Inches(6.4))
d.add_paragraph(
    "Figure 1 Reviewed marker-supported identities. UMAP distances and apparent branches are not developmental trajectories.",
    style="Caption",
)
rows = list(csv.DictReader((O / "cluster_annotation_reviewed.csv").open()))
t = d.add_table(rows=1, cols=4)
t.style = "Light Shading Accent 1"
for c, s in zip(
    t.rows[0].cells, ["Cluster", "Cells", "Reviewed identity", "Confidence"]
):
    c.text = s
for r in rows:
    for c, s in zip(
        t.add_row().cells,
        [r["cluster"], r["n_cells"], r["reviewed_label"], r["confidence"]],
    ):
        c.text = s
d.add_page_break()
d.add_heading("Evidence supporting the annotations", 1)
for r in rows:
    d.add_heading("Cluster " + r["cluster"] + " " + r["reviewed_label"], 2)
    d.add_paragraph(r["review_notes"])
d.add_page_break()
d.add_heading("Diagnostic gene expression", 1)
d.add_picture(str(F / "dotplot__review_diagnostic.png"), width=Inches(6.4))
d.add_paragraph(
    "Figure 2 Diagnostic markers across clusters. Dot size is the fraction expressing a gene; color is mean expression scaled per gene. Relative color intensity cannot compare absolute expression between genes.",
    style="Caption",
)
d.add_paragraph(
    "EPCAM, ESRP1 and GRHL2 support epithelial identities; KRT4, KRT13 and SCEL distinguish squamous differentiation. TTR, OTX2, KCNJ13 and SLC13A4 support the choroid plexus-like assignment. AFP, APOA1, HNF4A, CDX2, VIL1 and SI support endodermal differentiation. Shared markers were interpreted in combinations."
)
d.add_paragraph(
    "LHX2 and EMX2 support a forebrain-like progenitor program, but FOXG1 is weak and does not establish telencephalic identity. HOXB9, WNT1, PAX3 and MSX1 support posterior or dorsal patterning in cluster 5. Protocol and organoid age are needed to decide whether this regional mixture is expected."
)
d.add_page_break()
d.add_heading("Neuronal development", 1)
d.add_picture(str(F / "umap_neuronal_cells_development.png"), width=Inches(6.4))
d.add_paragraph(
    "Figure 3 Developmental and neurotransmitter-related gene expression within the immature neuronal population. The coordinates are inherited from the full-sample UMAP.",
    style="Caption",
)
d.add_paragraph(
    "DCX, TUBB3, STMN2 and NEUROD1, together with CELF3, SCRT2, SRRM4 and ELAVL3, support immature neuronal identity. MAP2 and synaptic genes provide complementary RNA-level evidence. This panel does not establish a resolved excitatory or inhibitory subtype or functional synaptic maturity."
)
d.add_paragraph(
    "Validate neuronal identity and organization using staining for progenitor and neuronal proteins. Electrophysiology or calcium imaging is needed for functional assessment. Cortical layer marker expression alone cannot establish organized cortical lamination."
)
d.add_page_break()
d.add_heading("Developmental programs across clusters", 1)
d.add_picture(str(F / "dotplot__development.png"), width=Inches(6.4))
d.add_paragraph(
    "Figure 4 Proliferation, regional, layer, synaptic, stress and hypoxia marker panels. Colors are scaled within each gene.",
    style="Caption",
)
d.add_paragraph(
    "Stress-rich cluster 10 has uncertain lineage specificity. Cluster 8 contains only 29 cells and may represent a neural transition or mixed state. Cluster 6 is doublet-enriched at approximately 31% and should not be interpreted as a resolved biological lineage. No mature astrocyte, oligodendrocyte or microglial population was confidently resolved by this analysis; absence of an assigned cluster is not proof of biological absence."
)
d.add_page_break()
d.add_heading("Quality control and robustness", 1)
d.add_picture(str(F / "qc_before.png"), width=Inches(6.4))
d.add_paragraph("Figure 5 Input cell quality distributions.", style="Caption")
d.add_picture(str(F / "qc_after.png"), width=Inches(6.4))
d.add_paragraph("Figure 6 Distributions after filtering.", style="Caption")
d.add_paragraph(
    "Reclustering unflagged cells using the existing PCA produced adjusted Rand index 0.809 against primary assignments. This check tests graph sensitivity after doublet exclusion; it does not repeat feature selection or PCA on singlets. ARI against vendor clusters was 0.630 for 8,611 overlapping cells. Resolution 0.3 and 1.0 gave ARIs of 0.683 and 0.711 against resolution 0.6. These comparisons support broad structure while showing that exact boundaries depend on analysis choices."
)
d.add_page_break()
d.add_heading("Reproducibility and rerunning the analysis", 1)
d.add_paragraph(
    "Follow the cross-platform setup in scripts/sample29/README.md and run the numbered sections in workflow.py and review.py. Run atlas.py separately for the reference annotation."
)
d.add_paragraph(
    "Create a sample-specific Python 3.13 environment using environment.yaml with Conda or requirements.lock.txt with venv. Current versions include Python 3.13.15 and Scanpy 1.12.4. Installed packages and downloaded artifacts are saved under environments/sample29 and environments/package_cache. These directories are local installations, not portable source archives."
)
d.add_paragraph(
    "Reproduction requires the three original gzipped matrix, barcode and feature files at data/single_cell_analysis/test_data. review.py also requires 29.scRNA.umap.txt for the vendor comparison. Preserve all supplied input files in a controlled archive. run_manifest.json provides SHA256 hashes for the three count inputs; retain a retrieval location and access instructions alongside the archive. Recheck these hashes before rerunning. The code records hashes but does not automatically reject mismatched inputs."
)
d.add_paragraph(
    "annotation_review.json stores the human-reviewed cluster mapping and evidence. Cluster numbers can change with parameters, software or data. Reassess this file after any such change. Fixed seeds and pinned dependencies improve repeatability on the same platform; bit-identical results across operating systems or numerical libraries are not guaranteed. Biological reproducibility requires independent organoids and donors."
)
d.add_heading("Files to save in Git", 1)
d.add_paragraph(
    "Commit the analysis scripts, annotation review, marker sources, and pinned environment files in scripts/sample29. Commit .gitignore and the project README. The existing .gitignore already permits scripts/sample29, but ignores results and plots by default."
)
d.add_paragraph(
    "For a compact analysis record, additionally track results/sample29/run_manifest.json, cluster_annotation_reviewed.csv, cell_type_composition.csv, resolution_comparison.csv, analysis_report.md and this Word report. Optionally track selected final PNG figures. Add narrow exceptions for these exact files to .gitignore if choosing to version them; no ignore rules or commits were changed for this report."
)
d.add_paragraph(
    "Keep environments/, .venv/, package caches, render intermediates and logs out of Git. Store the original count files, sample29_annotated.h5ad and large generated tables in a controlled data archive or Git LFS if permitted. The annotated object is about 239 MB and exceeds ordinary GitHub file limits. Never rely on a local environment directory as the sole dependency record. Keep protected genomic data outside a public repository."
)
d.add_page_break()
d.add_heading("Limitations and supporting files", 1)
d.add_paragraph(
    "Only one sample was analyzed. Organoid age, donor, target region, differentiation protocol and replicate metadata were not provided. Marker tests are exploratory one-cluster-versus-rest Wilcoxon tests with tie correction and Benjamini-Hochberg adjustment; clustering and testing use the same cells. These tests do not measure donor-level differential expression. Ambient RNA was not corrected because an unfiltered droplet matrix was unavailable. Relative cluster proportions reflect capture and filtering as well as tissue composition."
)
d.add_paragraph(
    "The H5AD integrates counts, expression, cell labels, QC, embeddings, graph and marker statistics. CSV exports provide cell metadata, UMAP coordinates, all ranked markers, top 30 markers, per-cluster marker expression, automatic and reviewed annotations, composition and clustering comparisons. The plots directory contains additional individual gene maps, PCA variance, doublet-score distributions and resolution comparisons. These machine-readable exports accompany this report and retain detail beyond the displayed figures."
)
d.add_heading("References", 1)
for s in [
    "Scanpy preprocessing and clustering https://scanpy.readthedocs.io/en/latest/tutorials/basics/clustering.html",
    "Scanpy Scrublet https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.scrublet.html",
    "Fair et al 2020 Electrophysiological maturation of cerebral organoids https://pmc.ncbi.nlm.nih.gov/articles/PMC7562943/",
    "Quadrato et al 2017 Cell diversity and network dynamics in photosensitive human brain organoids https://www.nature.com/articles/nature22047",
]:
    d.add_paragraph(s)
d.save(O / "sample29_integrated_report.docx")
