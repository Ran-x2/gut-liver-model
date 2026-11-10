"""Assemble text responses and paired legends for an existing completed package."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

def write_nichenet_input_readme(output):
    """Write a practical guide matching the user's result/support/audit layout."""
    folder = output / "06_nichenet_inputs"
    counts_path = folder / "output" / "nichenet_input_counts.tsv"
    if not counts_path.exists():
        return
    with counts_path.open(encoding="utf-8-sig", newline="") as handle:
        cases = list(csv.DictReader(handle, delimiter="\t"))

    def link(label, relative):
        if (folder / relative).is_file():
            return f"[{label}]({relative})"
        return f"{label} (not exported)"

    guide = """# NicheNet inputs: where to find what you need

This analysis documents the genes supplied to NicheNet for reviewer comment 6.
There are eight cases: two donors, each with gut T cells, gut macrophages,
liver T cells, and liver macrophages.

| Folder | What belongs here | When to open it |
|---|---|---|
| `output/` | Final count summary and reviewer response | **Start here** |
| `gene_lists/` | Exact selected receiver genes and any retained matching details | When inspecting individual genes |
| `backgrounds/` | Expressed-gene backgrounds and any retained detection details | When reproducing the scoring inputs |
| `audit/` | Input checks, reference fingerprints, and retained historical material | Technical traceability only |

## The two main deliverables

- [output/nichenet_input_counts.tsv](output/nichenet_input_counts.tsv): the eight-row summary of input genes, receiver cells, and background sizes.
- [output/response_06.txt](output/response_06.txt): the current reviewer-response text. The same answer is included in [responses_06_to_10.txt](../responses_06_to_10.txt).

## Understanding the supporting filenames

The existing prefixes identify the case: `3` = Donor 1 (AJD3280), `4` = Donor 2
(AJG2309); `gutT` = gut T cells, `gutMac` = gut macrophages, `liverT` = liver
T cells, and `liverMac` = liver macrophages. For example,
`3_gutT_in_nichenet.tsv` contains Donor 1 gut T-cell genes present in the
NicheNet target reference. A `.tsv` is a table with tab-separated columns.

| Filename ending | Contents | Role |
|---|---|---|
| `*_in_nichenet.tsv` | Selected genes present in the NicheNet target reference used here | **Final receiver gene set** supplied for ligand prioritization |
| `*_pre_nichenet.tsv` | Genes selected before checking reference availability | Optional intermediate detail |
| `*_not_in_nichenet.tsv` | Selected genes absent from that specific reference | Optional reference-coverage detail |
| `*_expressed_background.tsv` | Gene names detected in at least 1% of that donor's receiver cells | Reusable scoring input |
| `*_expression_detection.tsv` | Per-gene detection measurements, including zero-expression genes | Optional intermediate detail |

The NicheNet reference is a precomputed ligand-target matrix. Matching checks
whether selected gene identifiers occur in its target-gene rows; it does not
fit a new model or recalculate differential expression. Before-matching genes
split into genes in the reference and genes not in it. A header-only
`*_not_in_nichenet.tsv` means every selected gene was present in the reference.

Selection reused the saved DEG results with `q_value <= 0.05` and
`normalized_effect <= -0.5` for Donor 1 or `>= 0.5` for Donor 2. This is the saved
Monocle3 normalized effect, not a log2 fold change. Gene-list columns contain
the gene identifier, original comparison term, q value, and signed effect.
`directional_effect` in the final list orients positive values toward the named
donor. The comparison term stays the same because both lists come from opposite
directions of the same donor comparison.

In detection tables, `n_cells_expressing` counts cells with a positive source
count; `n_cells` is the total receiver-cell count; `fraction_expressing` divides
the former by the latter. `mean_source_count` averages saved source counts across
all those cells, including zeros. Background membership uses the fraction,
not the mean. For 1,717 cells the 1% rule requires at least 18 expressing cells;
for 67 or 15 cells, one expressing cell passes. Background genes need not be
DEGs and are not limited to receptor genes.

Detailed before-matching, absent-gene, and detection exports are **off by default**
in the updated script. Existing files are retained as supporting records.
Set `NICHE_EXPORT_INPUT_DETAILS=1` only when those detailed exports are needed.

## Audit files

These are not additional biological results:

"""
    audits = [
        ("original_input_manifest.tsv", "Original input paths, sizes, and SHA-256 fingerprints."),
        ("reference_nichenet_manifest.tsv", "Recovered reference paths and fingerprints; exact historical version identity is not assumed."),
        ("input_schema_audit.tsv", "Checks of genes, cells, donor labels, duplicate genes, and invalid values."),
        ("cached_background_concordance.tsv", "Comparisons with available original cached backgrounds; only cases with caches are included."),
    ]
    if (folder / "audit" / "response_06_draft.txt").exists():
        audits.append(("response_06_draft.txt", "Retained historical draft; use output/response_06.txt for the current answer. The updated script does not generate this draft."))
    for name, explanation in audits:
        guide += f"- {link(name, 'audit/' + name)}: {explanation}\n"
    guide += """
A fingerprint identifies exact file contents. A manifest records provenance;
its original paths describe where files were read, not extra copies to inspect.

## Counts and exact files

The count table retains its existing headers: `DEGs_before_model` means before
NicheNet reference matching, `DEGs_after_model` means genes in the reference,
and `removed_by_model` means genes not in it. `background_model_genes` is the
reference-matched background size. `unique_DEGs_after_model` checks uniqueness.
`positive_genes_outside_background` counts selected genes outside that background;
`scored_gene_universe` is their union, and `positive_prevalence` is the selected
fraction of that universe. `low_information` flags fewer than 20 input genes;
it is not a significance test. Donor 1 gut macrophages have only six input genes,
so interpretation of that case is limited.

Numbers below are gene counts linked to the corresponding files where exported.

| Donor and receiver | Before matching | In NicheNet | Not in NicheNet | Expressed background | Detection detail |
|---|---:|---:|---:|---:|---|
"""
    for case in cases:
        stem = case["case_id"]
        before = link(f"{int(case['DEGs_before_model']):,}", f"gene_lists/{stem}_pre_nichenet.tsv")
        present = link(f"{int(case['DEGs_after_model']):,}", f"gene_lists/{stem}_in_nichenet.tsv")
        absent = link(f"{int(case['removed_by_model']):,}", f"gene_lists/{stem}_not_in_nichenet.tsv")
        background = link(f"{int(case['expressed_background_genes']):,}", f"backgrounds/{stem}_expressed_background.tsv")
        detection = link("Full table", f"backgrounds/{stem}_expression_detection.tsv")
        guide += f"| {case['donor_label']}, {case['receiver_label']} | {before} | {present} | {absent} | {background} | {detection} |\n"
    (folder / "README.md").write_text(guide, encoding="utf8")

def write_nichenet_activity_readme(output):
    """Describe the final results and the separate technical audit records."""
    folder = output / "07_nichenet_activity"
    results = folder / "output"
    if not results.exists():
        return
    rankings = results / "ligand_activity_rankings_all_donors_and_cell_types.tsv"
    rows = None
    if rankings.exists():
        with rankings.open(encoding="utf-8-sig", newline="") as handle:
            rows = sum(1 for _ in csv.DictReader(handle, delimiter="\t"))

    def link(label, relative):
        if (folder / relative).is_file():
            return f"[{label}]({relative})"
        return f"{label} (not exported)"

    guide = """# NicheNet ligand activity: final results and audit records

**Start in `output/`.** It contains the reviewer response, two figures with
legends, one complete ranking table, and one summary table. `audit/` contains
technical checks and software records; you do not need it for routine review.

## Final results

| Result | What it tells you | Files |
|---|---|---|
"""
    guide += "| Reviewer comment 7 response | Current response text; also included in the combined responses file one folder above | " + link("Read response", "output/reviewer_comment_07_response.txt") + " |\n"
    figures = [
        ("Figure R7A", "The 15 highest-AUROC ligands for each donor and cell population, with the same AUROC scale in every panel", "Figure_R7A_top15_ligands_shared_AUROC_scale"),
        ("Figure R7B", "The distribution of all candidate-ligand AUROC scores for each donor and cell population", "Figure_R7B_all_ligand_AUROC_distributions"),
    ]
    for title, description, stem in figures:
        files = [link(label, f"output/{stem}{extension}") for label, extension in [("PNG figure", ".png"), ("Legend", "_legend.txt")]]
        guide += f"| {title} | {description} | {' / '.join(files)} |\n"
    guide += "| Complete ligand rankings | Every scored ligand, its scores, and its rank within its donor/cell-type case | " + link("All donors and cell types", "output/ligand_activity_rankings_all_donors_and_cell_types.tsv") + " |\n"
    guide += "| Top-ligand and score-range summary | Eight rows showing candidate counts, the top ligand, and the range of the top 15 AUROC scores | " + link("Read summary", "output/top_ligands_and_score_ranges_by_donor_and_cell_type.tsv") + " |\n"
    if rows is not None:
        guide += f"\nThe complete ranking table currently contains **{rows:,} rows** across eight donor/cell-type cases.\n"
    guide += """
Each figure has one 600-dpi PNG and one text legend. All figure text uses
Arial, the title is centered, and there is no subtitle.

## Reading the tables

The complete rankings are sorted by decreasing AUROC within each case.
`test_ligand` is the ligand name, `auroc` its activity score, and `rank_auroc`
its rank in that case. Other columns retain AUPR, corrected AUPR, Pearson
correlation, and the recorded score origin. `case_id` identifies the donor and
receiver cell population; `receiver` and `donor_code` repeat those components.

The retained data codes are `3` = Donor 1 (AJD3280), `4` = Donor 2 (AJG2309),
`gutT` = gut T cells, `gutMac` = gut macrophages, `liverT` = liver T cells, and
`liverMac` = liver macrophages. For example, `4_liverT` is Donor 2 liver T cells.
These are case identifiers, not file versions.

In the summary, `top15_min`, `top15_max`, and `top15_range` describe the spread
of the 15 highest AUROC scores. Small differences among candidates should not
be interpreted as identifying a uniquely causal ligand. Donor 1 gut macrophages
have only six input genes, which limits interpretation of that case.

## Technical audit files

These records are separated from the deliverables:

"""
    audits = [
        ("cached_ligand_score_validation.tsv", "Checks used to decide whether original cached ligand scores could be reused."),
        ("donor2_background_correction_summary.tsv", "Summary of changes after using Donor 2's own expressed-gene background; the retained case code is 4."),
        ("nichenet_activity_execution_log.txt", "Execution messages from the scoring workflow."),
        ("nichenet_function_source_manifest.tsv", "Exact source commit, URLs, and SHA-256 fingerprints of the NicheNet functions."),
        ("nichenet_software_and_analysis_scope.txt", "Software provenance, historical-version qualifications, and the scope of the analysis."),
        ("R_session_information.txt", "R and package-version information recorded during the scoring run."),
    ]
    for cell_type in ["gut_T_cells", "liver_T_cells", "liver_macrophages"]:
        name = f"donor2_{cell_type}_original_vs_corrected_ligand_scores.tsv"
        if (folder / "audit" / name).exists():
            audits.append((name, "Detailed original-versus-corrected scores for this Donor 2 cell population."))
    for name, explanation in audits:
        guide += f"- {link(name, 'audit/' + name)}: {explanation}\n"
    guide += """

## What future runs produce

Scripts overwrite these canonical filenames. The eight separate ranking files
are no longer exported because all their rows are in the complete ranking table.
The superseded response draft is not regenerated. Detailed
original-versus-corrected comparison exports are off by default; enable
`NICHE_EXPORT_SCORE_COMPARISONS=1` only when those audit details are needed.
The existing comparison files have been retained in `audit/` for traceability.
"""
    if (folder / "audit" / "retained_redundant_exports").exists():
        guide += '\n\nThe duplicate per-case tables and old response draft are grouped in\n`audit/retained_redundant_exports/`; skip this folder for routine review.\nAutomatic review blocked their deletion, so they are preserved outside the\nresults. The scripts do not regenerate these redundant exports.\n'
    (folder / "README.md").write_text(guide, encoding="utf8")


def write_receiver_targets_readme(output):
    """Explain the compact receiver-target deliverables without creating audits."""
    folder = output / "08_receiver_targets"
    results = folder / "output"
    if not results.exists():
        return

    def count_rows(name):
        with (results / name).open(encoding="utf-8-sig", newline="") as handle:
            return sum(1 for _ in csv.DictReader(handle, delimiter="\t"))

    guide = """# Receiver-target support: final results

**Open `output/` for the five final deliverables.** This folder addresses reviewer
comment 8: which observed upregulated receiver genes are linked to the prioritized
NicheNet ligands?

| File | Purpose |
|---|---|
| [reviewer_comment_08_response.txt](output/reviewer_comment_08_response.txt) | Current reviewer response; also included in the combined responses file one folder above |
| [Figure_R8_ligand_links_to_upregulated_receiver_genes.png](output/Figure_R8_ligand_links_to_upregulated_receiver_genes.png) | Figure showing receiver-gene support for the three highest-ranked ligands in each donor/cell-type case |
| [Figure_R8_ligand_links_to_upregulated_receiver_genes_legend.txt](output/Figure_R8_ligand_links_to_upregulated_receiver_genes_legend.txt) | Selection rules, visual encodings, and interpretation limits for the figure |
| [upregulated_receiver_gene_links_for_top15_ligands.tsv](output/upregulated_receiver_gene_links_for_top15_ligands.tsv) | Complete matching ligand-to-receiver-gene links for the top 15 ligands in each case |
| [upregulated_receiver_gene_counts_for_top15_ligands.tsv](output/upregulated_receiver_gene_counts_for_top15_ligands.tsv) | Number of matching receiver genes for every prioritized ligand, including zero-support cases |

## What the figure and tables contain

For each ligand, the analysis checks which of its top 250 NicheNet prior targets
are also in the actual directional receiver DEG set. The complete link table
covers the top 15 ligands per case. The figure displays the top three ligands
and at most five matching genes per ligand, ordered by prior regulatory potential.
The figure uses the complete tables; a separate duplicate plotting table is not
needed. Its PNG is 600 dpi, with Arial text, a centered title, and no subtitle.

`ligand` and `target` name the predicted link. `weight` is the NicheNet prior
regulatory potential. `q_value` and `normalized_effect` are saved receiver DE
results. `directional_effect` orients positive effects toward the indicated
donor. In the count table, `AUROC` is the ligand's ranking score and
`receiver_DEGs_among_top250_prior_targets` is the number of matching receiver
genes. `input_DEGs` is the size of that case's input gene set.

The retained table identifiers are `3` = Donor 1 (AJD3280), `4` = Donor 2
(AJG2309), `gutT` = gut T cells, `gutMac` = gut macrophages, `liverT` = liver
T cells, and `liverMac` = liver macrophages. For example, `4_liverT` means
Donor 2 liver T cells.

These receiver DEGs were also used to prioritize the ligands. Their agreement
with predicted targets supports transcriptional consistency; it does not
independently establish ligand expression, secretion, or causality.
"""
    links = results / "upregulated_receiver_gene_links_for_top15_ligands.tsv"
    counts = results / "upregulated_receiver_gene_counts_for_top15_ligands.tsv"
    if links.exists() and counts.exists():
        guide += f"\nCurrent results contain **{count_rows(links.name):,} gene links** and **{count_rows(counts.name):,} ligand/case count rows**.\n"
    guide += "\nScripts overwrite these canonical files and do not regenerate the old response draft or duplicate plotting subset.\n"
    (folder / "README.md").write_text(guide, encoding="utf8")


def write_annotation_readme(output):
    """Keep the annotation deliverables distinct from plotting inputs and checks."""
    folder = output / "09_annotation"
    if not (folder / "output").exists():
        return
    guide = """# Cell-type annotation: where to find what you need

**Start in `output/`.** It contains the final figures, annotation evidence,
Methods text, references, and reviewer response.

| Folder | Contents |
|---|---|
| `output/` | Final figures, shared figure legend, annotation evidence table and legend, proposed Methods, references, and reviewer response |
| `plot_data/` | One reusable table containing the exact marker-expression values behind all four figures |
| `audit/` | Source/schema checks, software and figure provenance, deposited clustering calls, and detailed reproduction instructions |

## Final figures

| Figure | Populations shown | PNG |
|---|---|---|
| R9A | Gut immune cells | [Figure_R9A_gut_immune_cell_markers.png](output/Figure_R9A_gut_immune_cell_markers.png) |
| R9B | Gut epithelial, endothelial, and stromal cells | [Figure_R9B_gut_epithelial_endothelial_stromal_markers.png](output/Figure_R9B_gut_epithelial_endothelial_stromal_markers.png) |
| R9C | Liver immune and other hematopoietic cells | [Figure_R9C_liver_immune_and_hematopoietic_markers.png](output/Figure_R9C_liver_immune_and_hematopoietic_markers.png) |
| R9D | Liver epithelial, endothelial, and stromal cells | [Figure_R9D_liver_epithelial_endothelial_stromal_markers.png](output/Figure_R9D_liver_epithelial_endothelial_stromal_markers.png) |

All four figures use [one shared legend](output/Figures_R9A_to_R9D_cell_type_marker_legend.txt).
Circle area shows the fraction of cells expressing each marker. Color shows
mean stored log-expression relative to that gene's maximum population mean
within the panel. Populations retain their existing annotations. Figures are
600-dpi PNGs with Arial text, centered titles, and no subtitles.

## Other final deliverables

| File | Purpose |
|---|---|
| [reviewer_comment_09_response.txt](output/reviewer_comment_09_response.txt) | Current reviewer response; also included in the combined responses file one folder above |
| [cell_type_annotation_methods.txt](output/cell_type_annotation_methods.txt) | Proposed manuscript Methods subsection |
| [cell_type_annotation_evidence_by_tissue_and_donor.csv](output/cell_type_annotation_evidence_by_tissue_and_donor.csv) | Original and displayed labels, stored cluster identifiers, donor-specific cell counts, and marker genes examined |
| [cell_type_annotation_evidence_by_tissue_and_donor_legend.txt](output/cell_type_annotation_evidence_by_tissue_and_donor_legend.txt) | Definitions and interpretation of the evidence-table fields |
| [cell_type_annotation_method_references.txt](output/cell_type_annotation_method_references.txt) | References supporting the Methods terminology |

The evidence table includes all three source donors. The dotplots show Donors
1 and 2, whose cells were used in the model. Marker lists record genes examined
in the original annotation work; they are not mandatory-positive gating rules.
The figures document existing labels and do not independently validate cell
function or TCR specificity.

The table preserves original source-label strings. The source name
"Peritoneal Epithelial Cells" denotes portal endothelial cells, which is the
corrected display label. Cell membership is unchanged.

The later liver notebook distinguishes CD8-alpha-alpha and CD8-alpha-beta MAIT
cells. The separate `analysis/revision/annotation/refine_liver_mait.py` script
can replay those restricted clustering calls on an H5AD's saved neighbor graph.
The checked H5AD gives different subclusters from the notebook's saved plot,
so the replay's candidate identities are not applied to the final figures.
The current figures retain the saved CD8 MAIT annotation. Run instructions and
the cluster comparison are in `audit/`.

## Supporting data and audit records

- [plot_data/marker_expression_by_cell_population.csv](plot_data/marker_expression_by_cell_population.csv): one row per plotted population/marker combination, with panel, tissue, source label, cell count, mean stored log-expression, expressing-cell fraction, and the relative value used for color.
- [audit/annotation_sources_and_figure_provenance.json](audit/annotation_sources_and_figure_provenance.json): source inputs, expression fields, diagnostic scope, and figure metadata.
- [audit/annotation_input_schema_checks.json](audit/annotation_input_schema_checks.json): input structure and source-field inspection records.
- [audit/deposited_clustering_calls_and_parameters.csv](audit/deposited_clustering_calls_and_parameters.csv): literal notebook calls, parameters, and cell locations; exploratory calls may be included.
- [audit/annotation_reproduction_instructions.txt](audit/annotation_reproduction_instructions.txt): detailed runtime and input requirements.
- [audit/liver_MAIT_refinement.json](audit/liver_MAIT_refinement.json): separate MAIT replay parameters, candidate cluster counts, and marker summaries; not accepted cell identities.

The plotting table uses `panel` values A-D to identify the four figures; it
replaces the four separate plotting-value files. The data are retained to allow
style changes without rereading the H5ADs or rerunning upstream analysis.

To redraw only the figures from this saved table, run:

```powershell
python "K:/gut-liver-model/analysis/revision/annotation/annotation_markers.py" --root "RUNTIME_ROOT" --output "K:/gut-liver-revision-output/09_annotation" --plots-only
```

Use the configured runtime directory for `RUNTIME_ROOT`. Scripts overwrite the
same filenames and keep a single shared figure caption.
"""
    plot_table = folder / "plot_data" / "marker_expression_by_cell_population.csv"
    if plot_table.exists():
        with plot_table.open(encoding="utf-8-sig", newline="") as handle:
            count = sum(1 for _ in csv.DictReader(handle))
        guide += f"\nThe current plotting table contains **{count:,} population/marker rows** across the four panels.\n"
    # Source-inspection and reproduction records may be retained from earlier work.
    guide = "\n".join(
        line for line in guide.splitlines()
        if not line.startswith("- [audit/") or (folder / line.split("](", 1)[1].split(")", 1)[0]).is_file()
    ) + "\n"
    (folder / "README.md").write_text(guide, encoding="utf8")


def write_cellphonedb_readme(output):
    folder = output / "10_cellphonedb"
    if not (folder / "output").exists():
        return
    guide = """# CellPhoneDB interaction counts

**Start in `output/`.** Donor 1 is AJD3280; Donor 2 is AJG2309.

| Folder | Contents |
|---|---|
| `output/` | Final count figure, its scientific legend, four count matrices, and reviewer response |
| `plot_data/` | Cell-population display names and presence/absence needed to redraw the figure |
| `source_figures/` | Four historical molecular-interaction heatmaps and their shared scientific legend |
| `audit/` | Source records, count checks, saved notebook tables, and lower-resolution source copies |

| Final file in `output/` | Contents |
|---|---|
| [Figure_R10_cellphonedb_interaction_counts_by_tissue_and_donor.png](output/Figure_R10_cellphonedb_interaction_counts_by_tissue_and_donor.png) | Four-panel gut/liver count heatmap with one shared scale and exact numbers |
| [Figure_R10_cellphonedb_interaction_counts_by_tissue_and_donor_legend.txt](output/Figure_R10_cellphonedb_interaction_counts_by_tissue_and_donor_legend.txt) | Methods, colors, abbreviations, and interpretation limits |
| [gut_donor_1_AJD3280_cellphonedb_interaction_counts.csv](output/gut_donor_1_AJD3280_cellphonedb_interaction_counts.csv) | Gut, Donor 1: 15 x 15 count matrix |
| [gut_donor_2_AJG2309_cellphonedb_interaction_counts.csv](output/gut_donor_2_AJG2309_cellphonedb_interaction_counts.csv) | Gut, Donor 2: 15 x 15 count matrix |
| [liver_donor_1_AJD3280_cellphonedb_interaction_counts.csv](output/liver_donor_1_AJD3280_cellphonedb_interaction_counts.csv) | Liver, Donor 1: 15 x 15 count matrix |
| [liver_donor_2_AJG2309_cellphonedb_interaction_counts.csv](output/liver_donor_2_AJG2309_cellphonedb_interaction_counts.csv) | Liver, Donor 2: 15 x 15 count matrix |
| [reviewer_comment_10_response.txt](output/reviewer_comment_10_response.txt) | Current response, also included in the combined responses one folder above |

Matrix counts sum both ordered partner assignments. The original numeric
values and labels are preserved. Self-pairs and absent populations retain
their stored zeros in the CSVs but appear gray in the figure. Liver Donor 2
has no hepatocytes or glial-like cells in the saved metadata. These are
predicted interaction counts, not calibrated signaling strength or a test of
donor differences. PEC means portal endothelial cells.

The final figure is a 600-dpi PNG with Arial text, centered panel titles, and
no subtitles. `source_figures/` preserves historical images unchanged: the
best available originals for gut Donors 1/2 and liver Donor 1, and the smaller
notebook image for liver Donor 2. Their formatting and independently selected
interaction panels differ from the final count figure.

To rebuild from the configured saved inputs, run from `K:/gut-liver-model`:

```powershell
python analysis/revision/build_package.py --steps cellphonedb
```

The scripts reuse saved counts without rerunning CellPhoneDB inference.
They overwrite the same files; source HTML export is optional and disabled
by default. Technical recovery limitations are recorded in `audit/`.
"""
    (folder / "README.md").write_text(guide, encoding="utf8")


def assemble(output):
    folders = {6: "06_nichenet_inputs", 7: "07_nichenet_activity", 8: "08_receiver_targets", 9: "09_annotation", 10: "10_cellphonedb"}
    responses = {
6: '''Response to comment 6

We used a separate directional receiver gene set for each donor in four populations: gut T cells, gut macrophages, liver T cells and liver macrophages. From the original saved donor-DE tables, genes were selected using q_value <= 0.05 and normalized_effect <= -0.5 for Donor 1 or >= +0.5 for Donor 2, then restricted to genes represented in the recovered NicheNet ligand-target matrix.

Receiver                 Donor 1    Donor 2
Gut T cells                  484        162
Gut macrophages                6         42
Liver T cells              1,133        823
Liver macrophages          1,068      1,232

The accompanying eight-row table gives counts before and after model matching, receiver cell numbers, background sizes and exact gene lists. These numbers refer to genes supplied for ligand prioritization, not the number of displayed ligands or genes passing the different volcano-plot threshold. The original DEG models were not refitted. The six Donor-1 gut-macrophage genes are ribosomal/translation-related, and interpretation of that panel is correspondingly limited.
''',
7: '''Response to comment 7

We agree that closely spaced ligand activities should not be interpreted as evidence for uniquely dominant drivers. NicheNet measures how well a ligand's prior target program recovers a receiver gene set; similar scores can arise from overlapping target programs. We provide the top-ranked ligand activities on a common AUROC scale, complete score distributions and full numerical rankings (Figures R7A and R7B). This presentation makes the narrow separation among T-cell candidates explicit.

The donor-specific expressed backgrounds have been corrected, while retaining the original receiver DEG results. We describe the results as candidate upstream regulators and predicted signaling relationships, rather than uniquely established causal drivers. In particular, a high activity score for the six-gene Donor-1 gut-macrophage set is not evidence for a specific immune signaling mechanism.
''',
8: '''Response to comment 8

We provide a receiver-target support panel and full ligand-target table (Figure R8). Predicted targets are restricted to the actual directional receiver DEG set, and the display includes the saved Monocle3 normalized_effect and q_value together with the ligand-target prior regulatory potential. This replaces the previous export of model-global top targets with links supported by observed receiver upregulation.

We distinguish ligand expression in possible sender populations, receptor expression in receivers, and downstream receiver-gene upregulation. The new panel provides transcriptional consistency with the predictions; because these receiver DEGs were used in ligand prioritization, it is not independent validation of ligand expression, secretion or causal signaling. Existing ligand-expression displays and the new receiver-target panel address different parts of the proposed signaling relationship.
''',
9: '''Response to comment 9

We provide a dedicated Methods subsection, a cluster-to-annotation evidence table and four marker-expression dotplots (Figures R9A-D). The Methods describe donor-specific quality control, doublet removal, integration, broad clustering, manual marker-guided annotation and lineage refinement. The evidence table links source labels to stored cluster identifiers, marker panels and donor-specific cell counts. Dot size shows the fraction expressing each marker and color shows per-gene relative mean stored log-expression.

These materials document the existing annotations without new clustering or reannotation. They distinguish initial clustering from lineage-specific refinement, three-donor preprocessing from two-donor presentation, and marker inspection from automated label transfer. The accompanying table retains the original source labels and the labels used in the figures.
''',
10: '''Response to comment 10

We provide a four-panel quantitative interaction-count heatmap for gut and liver from each donor (Figure R10). Each retained interaction count is printed, and all panels share one color scale. Self-pairs and cell populations absent from a donor's data are distinguished from tested pairs with zero retained interactions. Exact numerical matrices were recovered from complete saved tables and checked against available original outputs; CellPhoneDB inference was not rerun.

The legend states the expression/significance criteria, additional percentile filter and symmetrization rule. These values count retained predicted interaction entries, rather than calibrated signaling strength, and donor comparisons are descriptive. Original molecular-level heatmaps are supplied separately as optional supporting material with their own legend and provenance.
'''}
    output.mkdir(parents=True, exist_ok=True)
    combined = []
    for number, text in responses.items():
        folder = output / folders[number]
        if not folder.exists():
            continue
        if number in (6, 7, 8, 9, 10):
            folder = folder / "output"
            folder.mkdir(parents=True, exist_ok=True)
        response_name = f"reviewer_comment_{number:02d}_response.txt" if number in (7, 8, 9, 10) else f"response_{number:02d}.txt"
        (folder / response_name).write_text(text, encoding="utf8")
        combined.append(text.strip())
    (output / "responses_06_to_10.txt").write_text("\n\n" + ("\n\n" + "=" * 72 + "\n\n").join(combined) + "\n", encoding="utf8")
    readme = '''Reviewer response package 6-10

Start with responses_06_to_10.txt. Individual answers are in the numbered folders.

For comments 6-10, start with the README.md in the corresponding numbered folder. Final results and individual responses are in each output/ subfolder.

Figure R7A: 07_nichenet_activity/output/Figure_R7A_top15_ligands_shared_AUROC_scale
Figure R7B: 07_nichenet_activity/output/Figure_R7B_all_ligand_AUROC_distributions
Figure R8: 08_receiver_targets/output/Figure_R8_ligand_links_to_upregulated_receiver_genes.png
Figures R9A-D: 09_annotation/output/ (see 09_annotation/README.md for panel filenames)
Figure R10: 10_cellphonedb/output/Figure_R10_cellphonedb_interaction_counts_by_tissue_and_donor.png

Primary figures are exported as 600-dpi PNG only, with Arial text, centered titles, and no subtitles. Separate text legends accompany the figures. Exact count, expression, target-support and provenance tables remain in their respective folders. Historical molecular PNGs are supplied separately and were not enlarged to invent resolution.

R7-R10 are response-package identifiers. Insert the supplied Methods, tables and figures and assign final supplementary numbers before submitting the response. The original Word manuscript has not been overwritten. Review manuscript_methods_corrections.txt before insertion.

Final annotation figures retain the saved cell labels. A separate restricted MAIT-clustering replay is recorded in 09_annotation/audit/; its candidate labels have not been applied to final results. No donor-DE refitting, new normalization or CellPhoneDB inference was performed. The original DEG tables are reused; NicheNet backgrounds/ranking/target exports were corrected. Original source data remain unchanged.
'''
    (output / "00_README.txt").write_text(readme, encoding="utf8")
    write_nichenet_input_readme(output)
    write_nichenet_activity_readme(output)
    write_receiver_targets_readme(output)
    write_annotation_readme(output)
    write_cellphonedb_readme(output)
    files = []
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "output_manifest.json":
            files.append({"path": str(path.relative_to(output)), "bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    (output / "output_manifest.json").write_text(json.dumps(files, indent=2), encoding="utf8")
    print("Assembled", len(combined), "responses and", len(files), "output files.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    assemble(parser.parse_args().output)
