"""Write annotation revision prose and source-parameter evidence; no data analysis."""
import ast
import argparse
import csv
import json
from pathlib import Path

parser = argparse.ArgumentParser(description="Write annotation response/Methods text and deposited clustering-call evidence; no data analysis.")
parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1],
                    help="Runtime root containing the source repo or deposited notebooks.")
parser.add_argument("--gut", type=Path, default=Path(r"K:\model paper\v2\gut_annotated.h5ad"),
                    help="Gut source path recorded in the source-selection notes; no H5AD is read by this script.")
parser.add_argument("--liver", type=Path, default=Path(r"K:\model paper\liver_annotated-001.h5ad"),
                    help="Liver source path recorded in the source-selection notes; no H5AD is read by this script.")
parser.add_argument("--output", type=Path, help="Output directory; defaults to ROOT/output/09_annotation.")
parser.add_argument("--notebook-root", type=Path,
                    help="Directory containing gut_process.ipynb and liver_process.ipynb; defaults to ROOT/repo if present, otherwise ROOT.")
parser.add_argument("--gut-notebook", type=Path,
                    help="Explicit gut processing notebook path; overrides --notebook-root for this file.")
parser.add_argument("--liver-notebook", type=Path,
                    help="Explicit liver processing notebook path; overrides --notebook-root for this file.")
args = parser.parse_args()
ROOT = args.root.expanduser().resolve()
OUT = (args.output or ROOT / "output" / "09_annotation").expanduser().resolve()
NOTEBOOK_ROOT = (args.notebook_root or (ROOT / "repo" if (ROOT / "repo").is_dir() else ROOT)).expanduser().resolve()
RESULTS = OUT / "output"
AUDIT = OUT / "audit"
for directory in (RESULTS, AUDIT):
    directory.mkdir(parents=True, exist_ok=True)

TEXTS = {
    'cell_type_annotation_methods.txt': """Proposed Methods subsection: Single-cell clustering and cell-type annotation

Colon and liver datasets were processed and annotated separately after integration across donors. Donor-specific expression was normalized using SCTransform v2 with S-phase score, G2/M-phase score and mitochondrial fraction included as covariates. Integration used 3,000 selected features. Neighborhood graphs were constructed using 50 principal components and 15 neighbors. Initial broad Leiden clustering used a resolution of 0.3. UMAP embeddings were computed with PAGA initialization, a minimum distance (min_dist) of 0.001, and spread of 1.5. Selected broad populations were further examined by lineage-specific and restricted subclustering; the corresponding notebook parameters are supplied with the analysis code. The refinement resolutions were population dependent and should be distinguished from the initial broad-clustering resolution.

Cell identities were assigned manually by evaluating coherent combinations of canonical marker expression, inspecting gene-expression maps and dotplots, and examining cluster-level marker profiles. Broad epithelial, stromal, endothelial and immune identities were assigned first, followed by more detailed annotations within selected populations. Examples of markers examined included EPCAM and epithelial differentiation markers such as LGR5, SMOC2, ASCL2, MUC2, CHGA and POU2F3; PECAM1, CD34 and MADCAM1 for endothelial populations; ACTA2 and stromal markers including PDGFRA, PDGFRB and NRG1; and SOX10/S100B for enteric glial populations. Liver-specific annotations were informed by published human liver single-cell profiles, including MacParland et al. (2018), and expression patterns involving ALB/GSTA2/AKR1C1, KRT19/SOX9, endothelial markers and ACTA2/COL1A1.

Immune annotations were evaluated using T-lineage genes CD3D, CD3E, CD3G and CD247; B/plasma-cell markers including CD19, MS4A1 and immunoglobulin genes; and myeloid/dendritic-cell markers including CD14, ITGAM, CD68, FCGR3A, CD163, MRC1, CD1C and CLEC4C. Natural-killer and innate-lymphoid annotations were assessed using cytotoxic or NK-receptor expression, including NKG7 and NCR-family genes, and ILC-associated markers such as KIT, RORA and RORC, interpreted jointly with T-lineage expression. T-cell refinement considered TRDC/TRGC-family expression, CD4/CD8A/CD8B, and state-associated markers including CCR7, SELL, TCF7, IL7R, S1PR1, KLF2, ITGAE, ITGA1, KLRG1 and cytotoxic genes. FOXP3 and associated regulatory markers informed regulatory T-cell labels, and TRAV1-2 contributed to the existing MAIT-cell annotations. Annotation relied on the overall marker pattern; no universal single-marker positivity rule or automated reference-label transfer was applied.

The supplementary annotation table provides the original broad and fine labels, stored cluster identifiers, markers examined and per-donor cell counts. The accompanying dotplots summarize the existing annotations for Donors 1 and 2, whose cells were used to construct the model. Stored log-transformed expression was used directly for visualization, with no new normalization, clustering, annotation or differential-expression analysis. Circle area represents the proportion of cells with expression greater than zero. Color represents the mean stored log expression divided by the largest population mean for the same gene within that panel. These displays document the transcriptional evidence supporting the existing labels.
""",
    'Figures_R9A_to_R9D_cell_type_marker_legend.txt': """Supplementary Figure. Marker expression supporting the existing colon and liver cell annotations.

(A) Gut immune populations. (B) Gut epithelial, endothelial and stromal populations. (C) Liver immune and other hematopoietic populations. (D) Liver epithelial, endothelial and stromal populations. Data from Donors 1 and 2 were pooled within each tissue: 11,654 gut cells (5,026 and 6,628 from Donors 1 and 2) and 14,036 liver cells (7,081 and 6,955). Panels A-D contain 4,957, 6,697, 13,155 and 881 cells, respectively. The number of cells assigned to each population is shown beside its label.

Circle area is proportional to the percentage of cells with stored log-transformed expression greater than zero. Color shows relative mean expression: the mean stored log-transformed expression for a population divided by the largest population mean for that gene within the displayed panel (0-1). Thus color summarizes the distribution of each marker across populations and does not compare absolute expression between different genes or tissues. Values were obtained directly from the saved log-expression layer; no additional normalization or imputation was applied.

Populations are labeled using the existing source annotations. The accompanying table provides population labels, cluster identifiers, donor-specific cell counts, and display abbreviations. LSEC, liver sinusoidal endothelial cells; PEC, portal endothelial cells; ILC, innate lymphoid cells; NK, natural killer cells; TCM, central memory T cells; TEM, effector memory T cells; TRM, tissue-resident memory T cells; DN, CD4/CD8 double negative; Treg, regulatory T cells; MAIT, mucosal-associated invariant T cells. The liver MAIT row uses the available source annotation without a CD8-alpha-alpha/CD8-alpha-beta subdivision. These descriptive marker displays do not constitute independent functional or TCR-specificity validation.
""",
    'cell_type_annotation_evidence_by_tissue_and_donor_legend.txt': """Supplementary Table. Existing annotation labels, stored cluster identifiers and marker evidence.

Each row represents one existing high-resolution cell label in a tissue. Exact source label strings and their associated broad labels are retained alongside the figure display label. The table lists observed global/refined cluster identifiers stored in the supplied objects, cell counts for all three donors, the number of cells included in the Donor 1/2 marker plots, and marker genes examined in the source analysis. Stored broad cluster identifiers do not uniquely specify every later lineage-specific split; the deposited processing notebooks provide those additional manual mapping steps. Marker lists describe genes examined when interpreting populations and are not mandatory-positive gating rules. Zero counts indicate that a source population was absent from the indicated donor; populations absent from both Donors 1 and 2 are not shown in the marker figure. Cell labels and cluster assignments are those present in the supplied objects.

The original source label "Peritoneal Epithelial Cells" denotes portal endothelial cells.
""",
    'cell_type_annotation_method_references.txt': """References supporting the proposed methods terminology

The manuscript already cites the integration, Scanpy and liver annotation references below. Keep its established reference numbering when inserting the supplied paragraph.

Stuart T et al. Comprehensive Integration of Single-Cell Data. Cell177:1888-1902.e21 (2019). https://doi.org/10.1016/j.cell.2019.05.031
Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology19:15 (2018). https://doi.org/10.1186/s13059-017-1382-0
MacParland SA et al. Single cell RNA sequencing of human liver reveals distinct intrahepatic macrophage populations. Nature Communications9:4383 (2018). https://doi.org/10.1038/s41467-018-06318-7
Hafemeister C, Satija R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biology20:296 (2019). https://doi.org/10.1186/s13059-019-1874-1

Official SCT field definitions (supporting source-code interpretation, not new normalization): https://satijalab.org/seurat/reference/sctransform
Official SCT integration workflow: https://github.com/satijalab/seurat/blob/main/vignettes/seurat5_sctransform_v2_vignette.Rmd

No external atlas-based label-transfer procedure or new mandatory marker threshold was introduced. Population marker genes and label mappings are taken from the original processing/visualization notebooks and supplied annotations.
""",
}

for name, text in TEXTS.items():
    (RESULTS / name).write_text(text, encoding="utf8")

rows = []
for filename, explicit in [("gut_process.ipynb", args.gut_notebook),
                           ("liver_process.ipynb", args.liver_notebook)]:
    notebook_path = explicit.expanduser().resolve() if explicit else NOTEBOOK_ROOT / filename
    nb = json.loads(notebook_path.read_text(encoding="utf8"))
    for cell_index, cell in enumerate(nb["cells"], 1):
        source = "".join(cell.get("source", []))
        for line in source.splitlines():
            if "sc.tl.leiden(" not in line or line.lstrip().startswith("#"):
                continue
            tree = ast.parse(line.strip())
            calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
                     and isinstance(n.func, ast.Attribute) and n.func.attr == "leiden"]
            for call in calls:
                keywords = {kw.arg: ast.unparse(kw.value) for kw in call.keywords}
                rows.append({"source_notebook": notebook_path.name, "one_based_cell_number": cell_index,
                             "object": ast.unparse(call.args[0]),
                             "resolution": keywords.get("resolution", "default"),
                             "restrict_to": keywords.get("restrict_to", "none"),
                             "key_added": keywords.get("key_added", "default"),
                             "source_call": line.strip(),
                             "evidence_scope": "Literal deposited code; not rerun; exploratory calls may be present"})
with (AUDIT / "deposited_clustering_calls_and_parameters.csv").open("w", encoding="utf-8-sig", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader(); writer.writerows(rows)
print(f"Wrote {len(TEXTS)} text files and {len(rows)} source clustering-call records.")
