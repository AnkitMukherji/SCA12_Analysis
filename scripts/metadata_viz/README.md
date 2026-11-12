# Clinical metadata visualization scripts

These scripts analyze clinical/genetic workbooks separately from sample29 scRNA-seq. Documentation reflects current code and files on 12 November 2026; analyses were not rerun for this documentation task.

## Inputs and cleaning

| Workbook in data/metadata | Used by | Processing |
|---|---|---|
| `metadata_sca12.xlsx` (default first sheet) | Shared loader, plots 02–06, report | Strip column names; standardize Gender. Shared loader also strips participant IDs. Code explicitly sets AO and Disease duration to missing for Record ID 60; it does not exclude that participant from every analysis. |
| `358_DATA.xlsx`, sheet `Clean` | Plot 07, Sankey 08, report | Standardize first symptoms and map MRI, autonomic, bradykinesia, rigidity and plantar categories. This cohort is not automatically linked/merged with the first workbook. |

Required first-workbook fields include Record ID, Participant ID, Gender, Age, AO, Disease duration, LOWER ALLELE, UPPER ALLELE, SARA Score, INAS Score, FTMTRS Grand Total, TETRAS Grand Total and ADL Total. The second workbook requires the category columns named in `load_clean_data()`.

Unrecognized mapped categories generally become Unknown. Parasympathetic mild-to-moderate and sympathetic moderate-to-severe loss both map to Moderate Loss. Bradykinesia/rigidity recognize exact string keys `0`/`1`; other representations may become Unknown. First symptoms are stringified/title-cased without the same unknown mapping. Check these code choices when workbook encoding changes.

## Script and output inventory

Existing figures are in [plots/metadata_plots](../../plots/metadata_plots). Current script defaults differ from this location; see the path notes below.

| Script | Computation | Named outputs |
|---|---|---|
| `data_loader.py` | First-workbook cleaning; helper subsets requiring all five clinical scores or both allele sizes. Direct execution prints a loading summary, not a cleaned workbook. | Console summary |
| `02_plot_demographics.py` | Available-value age/onset distributions and gender counts/proportions. | `01_age_distribution.png`, `02_gender_distribution.png` |
| `03_plot_genetics.py` | Allele distributions/scatter on complete allele pairs. Its scatter rule calls LOWER ALLELE >=43 biallelic expansion, otherwise heterozygous expansion. This documents code behavior, not an independently reviewed diagnostic standard. | `03_allele_distribution.png`, `04_allele_scatter.png` |
| `04_plot_anticipation.py` | Complete upper-allele/onset pairs; regression, Pearson r, R², p-value and equation. SARA colours points when available; missing SARA is grey. | `05_genetic_anticipation.png` |
| `05_plot_progression.py` | Complete duration/SARA pairs; regression and Pearson statistics. Upper allele colours points when available; missing allele is grey. | `06_clinical_progression.png` |
| `06_plot_symptom_concordance.py` | Complete cases across five clinical scales; Pearson correlation matrix, not a p-value matrix. | `07_symptom_concordance.png` |
| `07_plot_proportions.py` | Second-workbook category counts and percentages. | `08a_first_symptom.png`, `08b_mri_categories.png`, `08c_autonomic_activity.png`, `08d_extrapyramidal.png`, `08e_plantar.png` |
| `08_plot_sankey.py` | Imports plot 07's loader; counts adjacent pairs across first symptom → MRI → parasympathetic → sympathetic → bradykinesia → plantar. Combines mild/moderate autonomic categories. Rigidity is not a Sankey stage. | `09_sankey_diagram.html`, `09_sankey_diagram.png` |
| `generate_report.py` | Builds a Word report from both cohorts, calculated summaries/regressions and existing figures. Includes fixed narrative as well as computed values; does not regenerate plots. | `observations_from_metadata.docx` |

Most Matplotlib figures use 600 dpi. Sankey HTML is interactive; PNG export uses Plotly's static-image backend, Kaleido and its compatible browser runtime. These plotting scripts do not export analysis CSV tables.

## Denominators and interpretation

Demographics use available values; allele plots require both sizes; each regression requires its two variables; scale correlations require all five scores. Category percentages generally use the entire cleaned second cohort, retaining Unknown rather than counting it as clinical absence. Plantar percentages divide by summed category counts.

Upper allele versus onset age is a cross-sectional association, not parent–offspring anticipation analysis. Duration versus SARA is not an individual longitudinal progression rate. Sankey links represent co-occurring categories, not observed temporal transitions or causal flows.

Some titles/report passages hard-code cohort sizes such as 74, 80 or 356, axis limits and narrative results. They need review after input changes. The correlation colour scale is fixed to 0.6–1.0; numeric values can fall outside the displayed colour range.

## Existing path and report mismatches

Many defaults use `os.path.dirname(os.path.dirname(__file__))`. From `scripts/metadata_viz` this resolves to **D:/sca12_analysis/scripts**, not the root:

| Coded default | Current layout / implication |
|---|---|
| `scripts/data/metadata/...` | Workbooks are under `data/metadata`; default loading can fail. |
| `scripts/plots/...` | Existing figures are under `plots/metadata_plots`. |
| `scripts/results/observations_from_metadata.docx` | Current report destination, not root results. |

Changing the working directory alone cannot fix paths based on `__file__`. Before a full rerun, align defaults with the project layout (the root from these scripts is `Path(__file__).resolve().parents[2]`) or pass explicit paths where supported. An explicit output path alone does not fix default workbook loading. Python scripts were not modified by this documentation update.

The report also looks for `07_clinical_proportions_dashboard.png` and `08_sankey_diagram.png`; current scripts produce five separate `08a`–`08e` images and `09_sankey_diagram.png`. Its existence checks can omit these figures. Fixed interpretive prose needs review rather than being treated as automatically verified by report generation.

## Running and dependencies

Imports require pandas, NumPy, Matplotlib, seaborn, SciPy and an Excel reader such as openpyxl. Sankey needs Plotly/Kaleido; the report needs python-docx. No dedicated portable metadata lockfile is present here; root requirements include machine-specific build paths.

After resolving path mismatches, run plotting entrypoints numerically, then `generate_report.py` if needed. Numbering starts at 02; `data_loader.py` is the shared utility, not a missing 01 stage. These older scripts use functions; sample29 uses sequential blocks.

For a read-only first-workbook load with an explicit path, from the root in a suitably provisioned Python session:

```python
from scripts.metadata_viz.data_loader import load_and_clean_data
df = load_and_clean_data('data/metadata/metadata_sca12.xlsx')
```

This cleans in memory without modifying the workbook or regenerating plots. See the [project overview](../../README.md) for the separate single-cell workflow.
