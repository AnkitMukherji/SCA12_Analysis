# Sources for the manual marker panels

Sample 29 is a **cerebral organoid**. The analysis does not assume cortical
identity, anatomical origin within the cerebrum, or a particular induction
protocol. Layer-associated genes are exploratory patterning evidence.

The source papers and their marker tables support the **program-level**
interpretations below. The panels in `workflow.py` are curated combinations,
not verbatim lists from a single paper. The additional ranked genes in
`annotation_review.json` are **observations from sample 29**. They come from
`cluster_markers_all.csv` and do not all have external specificity evidence.

| Manual call or panel | Representative evidence | Primary reference and where to look |
| --- | --- | --- |
| Neural progenitors and radial glia (5, 7, 12) | SOX2, PAX6, NES, VIM, HES1/HES5, FABP7 | [Nowakowski et al. 2017](https://doi.org/10.1126/science.aap8809), fetal progenitor populations and supplementary WGCNA markers. [Fair et al. 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7562943/), Figure 4 and Figure S2 for PAX6 and SOX2 in organoids. |
| Intermediate progenitors and immature neurons (8, 11) | EOMES, NEUROG2, INSM1 for neurogenesis, DCX, TUBB3, STMN2 and NEUROD1 for immature neuronal profiles | [Nowakowski et al. 2017](https://doi.org/10.1126/science.aap8809), excitatory neurogenesis and intermediate progenitor populations. [Rosebrock et al. 2022](https://www.nature.com/articles/s41467-022-33364-z), neuronal differentiation and neuronal immunostaining. |
| Excitatory and inhibitory lineages and layer-associated programs | SLC17A7, TBR1, BCL11B, SATB2 and CUX2 versus GAD1/GAD2, DLX1/DLX2 and SLC32A1 | [Nowakowski et al. 2017](https://doi.org/10.1126/science.aap8809), excitatory and inhibitory trajectories. [Quadrato et al. 2017](https://doi.org/10.1038/nature22047), cerebral organoid cell diversity and projection neuron staining. |
| Regional patterning (5, 7, 12) | FOXG1, LHX2, EMX2 and PAX6 versus the posterior-associated HOXB9, WNT1, PAX3 and MSX1 program | [Braun et al. 2023](https://doi.org/10.1126/science.adf1226), cross-region developmental marker tables. [He et al. 2024](https://doi.org/10.1038/s41586-024-08172-8), Figure 2, organoid regional identity across protocols. |
| Astroglia, oligodendrocyte lineage and microglia screening | GFAP/AQP4/SLC1A3, OLIG2/SOX10/PDGFRA, and P2RY12/CX3CR1/AIF1 | [Braun et al. 2023](https://doi.org/10.1126/science.adf1226), class and cluster metadata. No clear clusters were assigned these labels in sample 29. |
| Epithelial and squamous epithelial (0–2) | EPCAM and epithelial keratins, with KRT4 and KRT13 supporting squamous differentiation | [Human oral mucosa spatial and single-cell analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC13530142/), KRT4/KRT13 expression across nonkeratinized epithelial layers. These markers do not establish oral tissue identity here. |
| Mesenchymal and vascular-associated (3, 6, 9) | COL1A1/COL1A2, DCN and LUM versus ACTA2/PDGFRB and PECAM1/CLDN5 | [Wälchli et al. 2024](https://doi.org/10.1038/s41586-024-07493-y), fetal and adult brain vascular and perivascular cell types. [Winkler et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8995178/), CLDN5, PDGFRB and DCN staining. |
| Choroid plexus-like epithelium (4) | TTR, OTX2, KCNJ13 and SLC13A4 with epithelial expression | [Pellegrini et al. 2020](https://doi.org/10.1126/science.aaz5626), choroid plexus organoid identity and supplementary scRNA-seq data, [GSE150903](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150903). TTR can also occur in hepatic-like populations. |
| Stress and hypoxia (10) | HSP-family genes and glycolytic programs, with PGK1 as one organoid stress-associated example | [Bhaduri et al. 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7433012/), comparison of primary tissue and organoid stress programs. HSPA6, HSP90AA1 and CRYAB were ranked in sample 29. |
| Endodermal or hepatointestinal-like (13) | AFP, ALB, APOA1, HNF4A alongside CDX2, VIL1 and SI | [Wesley et al. 2022](https://doi.org/10.1038/s41556-022-00989-7), hepatoblast developmental trajectories and Figure 4 marker validation. [Elmentaite et al. 2021](https://doi.org/10.1038/s41586-021-03852-1), intestinal epithelial developmental atlas. |
| Proliferation, synaptic maturation, pluripotency and neural crest screening | MKI67/TOP2A, SNAP25/SYP/MAP2, POU5F1/NANOG, SOX10/TFAP2A | [Fair et al. 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7562943/), proliferation and maturation assays. [Petropoulos et al. 2016](https://pubmed.ncbi.nlm.nih.gov/27662094/), embryo lineage expression. [Braun et al. 2023](https://doi.org/10.1126/science.adf1226), developmental classes including neural crest-related populations. |

Detection and expression are useful in combination. They cannot independently
establish the identity, anatomical localization, or function of a cell.
