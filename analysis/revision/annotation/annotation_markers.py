"""Render descriptive annotation evidence from saved labels and expression.

No normalization, clustering, differential expression, or annotation is run.
Read only the designated H5ADs; separate final results, plotting data, and audits.
"""
from __future__ import annotations

import collections
import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path

def parse_cli():
    parser = argparse.ArgumentParser(
        description="Redraw annotation marker dotplots from saved H5AD labels and log expression; no new normalization, clustering or DE."
    )
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1],
                        help="Runtime root; its python-packages directory is added to the Python import path.")
    parser.add_argument("--gut", type=Path, default=Path(r"K:\model paper\v2\gut_annotated.h5ad"),
                        help="Gut H5AD with existing celltype_lowres/highres and layers/logcounts.")
    parser.add_argument("--liver", type=Path, default=Path(r"K:\model paper\liver_annotated-001.h5ad"),
                        help="Liver H5AD with existing celltype_lowres/highres and layers/logcounts.")
    parser.add_argument("--output", type=Path,
                        help="Output directory; defaults to ROOT/output/09_annotation.")
    parser.add_argument("--plots-only", action="store_true", help="Redraw from plot_data/marker_expression_by_cell_population.csv without reading H5AD data.")
    return parser.parse_args()

ARGS = parse_cli()
ROOT = ARGS.root.expanduser().resolve()
sys.path.insert(0, str(ROOT / "python-packages"))
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, colors

OUT = (ARGS.output or ROOT / "output" / "09_annotation").expanduser().resolve()
RESULTS = OUT / "output"
PLOT_DATA = OUT / "plot_data"
AUDIT = OUT / "audit"
for directory in (RESULTS, PLOT_DATA, AUDIT):
    directory.mkdir(parents=True, exist_ok=True)
PANEL_FILENAMES = {
    "A": "Figure_R9A_gut_immune_cell_markers",
    "B": "Figure_R9B_gut_epithelial_endothelial_stromal_markers",
    "C": "Figure_R9C_liver_immune_and_hematopoietic_markers",
    "D": "Figure_R9D_liver_epithelial_endothelial_stromal_markers",
}
PLOTTING_TABLE = PLOT_DATA / "marker_expression_by_cell_population.csv"
windows_fonts = Path(os.environ.get("WINDIR", r"C:\Windows")) / "Fonts"
for font in ["arial.ttf", "arialbd.ttf", "ariali.ttf"]:
    font_path = windows_fonts / font
    if font_path.is_file():
        font_manager.fontManager.addfont(font_path)
plt.rcParams.update({"font.family": "Arial", "font.size": 9.5,
                     "figure.facecolor": "white", "axes.facecolor": "white",
                     "savefig.facecolor": "white", "axes.linewidth": 0.7})

CONFIG = {
    "gut": {
        "path": ARGS.gut.expanduser().resolve(),
        "reference": "X",
        "immune_lowres": {"T Cells", "NK", "ILC", "B Cells", "Plasma Cells",
                          "Monocytes", "Macrophages", "Dendritic Cells"},
        "panels": [
            ("A", "gut_immune", "Gut: immune populations",
             "CD3D CD3E CD3G CD247 CD4 CD8A CD8B CCR7 IL7R KLF2 ITGAE FOXP3 KLRG1 NKG7 TRAV1-2 TRDC RORA KIT NCR3 CD19 IGHA1 ITGAM CD68 CD86 FCGR3A CD163 MRC1 CD1C".split()),
            ("B", "gut_structural", "Gut: epithelial, endothelial and stromal populations",
             "EPCAM LGR5 SMOC2 ASCL2 SLC26A2 SPIB MUC2 CHGA POU2F3 CCL21 PECAM1 CD34 MADCAM1 NRG1 PDGFRA PDGFRB ACTA2 DCLK1 BMP7 WNT5A SOX10 S100B".split()),
        ],
    },
    "liver": {
        "path": ARGS.liver.expanduser().resolve(),
        "reference": "layers/counts",
        "immune_lowres": {"T Cells", "NK", "ILC", "B Cells", "Plasma Cells",
                          "Macrophages", "Dendritic Cells", "Mast Progenitor", "Erthyroid"},
        "panels": [
            ("C", "liver_immune", "Liver: immune and hematopoietic populations",
             "CD3D CD3E CD3G CD247 CD4 CD8A CD8B CCR7 IL7R KLF2 ITGAE ITGA1 KLRG1 IFNG NKG7 TRAV1-2 TRDC TBX21 RORA RORC KIT IL23R EOMES NCR3 CD19 CD27 IGHA1 ITGAM CD68 CD86 FCGR3A CD163 TNF MRC1 IL1B IL1R1 CD1C CLEC4C TPSAB1 GATA2 HBB CA1".split()),
            ("D", "liver_structural", "Liver: epithelial, endothelial and stromal populations",
             "EPCAM SOX9 KRT19 HBB CA1 S100B GSTA2 AKR1C1 ALB GATA2 CLEC14A MGP CCL14 KIT CD34 PECAM1 RAMP3 LIFR ACTA2 COL1A1".split()),
        ],
    },
}

ORDER = {
    "gut": ["TCRab DN Mobile TRM", "TCRab CD4 Mobile TRM", "TCRab CD4 TRM",
            "TCRab CD4 FOXP3+ Treg", "TCRab CD8 Mobile TRM", "TCRab CD8 TRM",
            "TCRab CD8ab TRM", "TCRab CD8aa MAIT", "TCRgd TRM", "ILC1", "ILC3",
            "NK", "Naive B Cells", "Memory B Cells", "Plasma Cells",
            "CD16- Conventional Monocytes", "M2-like CD16+ Macrophages", "cDC2",
            "Intestinal Stem Cells", "Transit Amplifying Cells", "Colonocytes",
            "SPIB+ Cells", "Goblet Cells", "Enteroendocrine Cells", "Tuft Cells",
            "MAdCAM1- Endothelial Cells", "MAdCAM1+ Endothelial Cells",
            "Lymphatic Endothelial Cells", "Fibroblastic Reticular Cells",
            "Fibroblast", "Telocytes", "Vascular Smooth Muscle Cells",
            "Muscularis Smooth Muscle Cells", "Enteric Glial Cells"],
    "liver": ["TCRab DNT", "TCRab CD4 TCM", "TCRab CD4 TRM", "TCRab CD8 TCM",
              "TCRab CD8 TEM", "TCRab CD8 TRM", "TCRab CD8 effector",
              "TCRab CD8 MAIT", "TCRgd TRM", "TCRgd effector", "ILC1-like Cells",
              "KIT+ ILC3", "KIT- ILC3-like Cells", "Regulatory ILC", "NK", "IL7R+ NK",
              "Naive B Cells", "IgM Naive B", "Memory B Cells", "Plasma Cells",
              "CD16- Macrophage", "M1-like CD16+ Macrophage", "M2-like CD16+ Macrophage",
              "M2-like CD16- Macrophage", "Immature DC", "pDC", "Mast Progenitor",
              "Erthyroid", "Hepatocytes Type 6", "Cholangiocytes",
              "Liver Sinusoidal Endothelial Cells zone 1", "Liver Sinusoidal Endothelial Cells zone 2",
              "Peritoneal Epithelial Cells", "Stellate Cells", "Glial-like Cells"],
}

DISPLAY = {
    "TCRab DN Mobile TRM": "TCRab DN mobile TRM",
    "TCRab CD4 Mobile TRM": "TCRab CD4 mobile TRM",
    "TCRab CD8 Mobile TRM": "TCRab CD8 mobile TRM",
    "TCRab CD4 FOXP3+ Treg": "TCRab CD4 FOXP3+ Treg",
    "M2-like CD16+ Macrophages": "M2-like CD16+ macrophages",
    "CD16- Conventional Monocytes": "CD16- monocytes",
    "Liver Sinusoidal Endothelial Cells zone 1": "LSEC zone 1",
    "Liver Sinusoidal Endothelial Cells zone 2": "LSEC zone 2",
    "MAdCAM1- Endothelial Cells": "MAdCAM1- endothelial cells",
    "MAdCAM1+ Endothelial Cells": "MAdCAM1+ endothelial cells",
    "Lymphatic Endothelial Cells": "Lymphatic endothelial cells",
    "Fibroblastic Reticular Cells": "Fibroblastic reticular cells",
    "Vascular Smooth Muscle Cells": "Vascular smooth muscle cells",
    "Muscularis Smooth Muscle Cells": "Muscularis smooth muscle cells",
    "Peritoneal Epithelial Cells": "Portal endothelial cells",
    "Erthyroid": "Erythroid",
}

def decode(value):
    if isinstance(value, bytes):
        return value.decode("utf8")
    if isinstance(value, np.generic):
        return value.item()
    return value

def read_column(obj):
    if isinstance(obj, h5py.Dataset):
        return np.array([decode(v) for v in obj[()]], dtype=object)
    cats = [decode(v) for v in obj["categories"][()]]
    return np.array([cats[int(i)] if int(i) >= 0 else None
                     for i in obj["codes"][()]], dtype=object)

def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def source_markers(label, tissue):
    """Genes examined in original marker panels; not an automated assignment rule."""
    if label.startswith("TCR"):
        return "CD3D;CD3E;CD3G;CD247;CD4;CD8A;CD8B;TRDC;TRAV1-2;CCR7;IL7R;KLF2;ITGAE;FOXP3"
    if "Macrophage" in label:
        return "CD68;FCGR3A;CD163;MRC1;CD86;ITGAM"
    if "Monocyte" in label:
        return "CD14;ITGAM;CD68;FCGR3A"
    if label in {"pDC", "Immature DC", "cDC2"}:
        return "CD1C;CLEC4C;CD68;CD86"
    if "ILC" in label or "NK" in label:
        return "CD3E;NKG7;NCR3;RORA;RORC;KIT;IL23R;EOMES"
    if "B Cell" in label or label == "IgM Naive B":
        return "CD19;CD27;IGHM;IGHD;MS4A1"
    if "Plasma" in label:
        return "IGHA1;CD19;CD27"
    if label == "Mast Progenitor":
        return "TPSAB1;GATA2;KIT"
    if label == "Erthyroid":
        return "HBB;CA1;GATA2"
    if "Glial" in label:
        return "S100B;SOX10" if tissue == "gut" else "S100B"
    if "Hepatocyte" in label:
        return "ALB;GSTA2;AKR1C1"
    if label == "Cholangiocytes":
        return "EPCAM;SOX9;KRT19"
    if "Endothelial" in label or "Sinusoidal" in label:
        return "PECAM1;CD34;MADCAM1;CCL21" if tissue == "gut" else "PECAM1;CD34;CLEC14A;MGP;CCL14;RAMP3;LIFR"
    if label == "Peritoneal Epithelial Cells":
        return "PECAM1;CD34;CLEC14A;MGP;CCL14;EPCAM;KRT19"
    if label in {"Stellate Cells", "Fibroblast", "Telocytes", "Fibroblastic Reticular Cells"} or "Smooth Muscle" in label:
        return "ACTA2;COL1A1;PDGFRA;PDGFRB;NRG1;BMP7;WNT5A"
    return "EPCAM;LGR5;SMOC2;ASCL2;SLC26A2;SPIB;MUC2;CHGA;POU2F3"

def aggregate(tissue, cfg):
    with h5py.File(cfg["path"], "r") as f:
        obs = f["obs"]
        high = read_column(obs["celltype_highres"])
        low = read_column(obs["celltype_lowres"])
        cluster_columns = {key: read_column(obs[key]) for key in ["leiden", "leiden2", "leiden3"] if key in obs}
        batch = read_column(obs["batch"])
        genes = read_column(f["var"][decode(f["var"].attrs["_index"])])
        selected = np.isin(batch, ["3", "4"])
        available = set(high[selected])
        groups = [g for g in ORDER[tissue] if g in available]
        groups += sorted(available - set(groups))
        group_map = {g: i for i, g in enumerate(groups)}
        group_index = np.array([group_map[h] if take else -1
                                for h, take in zip(high, selected)], dtype=np.int32)
        gene_list = list(dict.fromkeys(gene for panel in cfg["panels"] for gene in panel[3]))
        missing = set(gene_list) - set(genes)
        if missing:
            raise ValueError(f"Missing requested {tissue} markers: {sorted(missing)}")
        requested = {gene: i for i, gene in enumerate(gene_list)}
        col_map = np.array([requested.get(gene, -1) for gene in genes], dtype=np.int32)
        sums = np.zeros((len(groups), len(gene_list)), dtype=np.float64)
        positive = np.zeros_like(sums, dtype=np.int64)
        n_cells = np.bincount(group_index[selected], minlength=len(groups))
        layer = f["layers/logcounts"]
        reference = f[cfg["reference"]]
        indptr = layer["indptr"][:]
        if not np.array_equal(indptr, reference["indptr"][:]):
            raise ValueError("Cached log expression and reference have different sparse structures")
        n_compared = 0
        max_error = 0.0
        n_reference_mismatches = 0
        for start in range(0, len(high), 256):
            end = min(start + 256, len(high))
            if not selected[start:end].any():
                continue
            first, last = int(indptr[start]), int(indptr[end])
            indices = layer["indices"][first:last]
            reference_indices = reference["indices"][first:last]
            if not np.array_equal(indices, reference_indices):
                raise ValueError("Cached log expression and reference have different gene indices")
            local_rows = np.repeat(np.arange(start, end), np.diff(indptr[start:end + 1]))
            columns = col_map[indices]
            keep = (columns >= 0) & selected[local_rows]
            if not keep.any():
                continue
            values = layer["data"][first:last][keep]
            raw_reference = reference["data"][first:last][keep]
            if not np.isfinite(values).all() or (values < 0).any():
                raise ValueError(f"Non-finite or negative selected marker values in {tissue}")
            expected = np.log1p(raw_reference)
            if not np.isfinite(expected).all():
                raise ValueError(f"Non-finite marker reference values in {tissue}")
            error = float(np.max(np.abs(values - expected)))
            max_error = max(max_error, error)
            n_compared += len(values)
            n_reference_mismatches += int(np.sum(np.abs(values - expected) > 2e-5))
            gi = group_index[local_rows[keep]]
            ci = columns[keep]
            np.add.at(sums, (gi, ci), values)
            np.add.at(positive, (gi, ci), values > 0)
        means = sums / n_cells[:, None]
        fractions = positive / n_cells[:, None]
        metadata_rows = []
        for label in sorted(set(high)):
            lowres_values = sorted(set(low[high == label]))
            metadata_rows.append({
                "tissue": tissue, "source_highres_label": label,
                "figure_display_label": DISPLAY.get(label, label),
                "source_lowres_labels": ";".join(lowres_values),
                "source_global_leiden_ids": ";".join(sorted(set(cluster_columns.get("leiden", np.array([""] * len(high)))[high == label]))),
                "source_refined_leiden2_ids": ";".join(sorted(set(cluster_columns.get("leiden2", np.array([""] * len(high)))[high == label]))),
                "source_refined_leiden3_ids": ";".join(sorted(set(cluster_columns.get("leiden3", np.array([""] * len(high)))[high == label]))),
                "donor_1_batch_3_n": int(np.sum((high == label) & (batch == "3"))),
                "donor_2_batch_4_n": int(np.sum((high == label) & (batch == "4"))),
                "donor_3_batch_5_n": int(np.sum((high == label) & (batch == "5"))),
                "figure_n_donors_1_2": int(np.sum((high == label) & selected)),
                "marker_genes_examined_in_source_panels": source_markers(label, tissue),
                "annotation_evidence_type": "Existing manually assigned label; marker expression shown descriptively",
                "source_h5ad": str(cfg["path"]),
                "note": 'The original source label "Peritoneal Epithelial Cells" denotes portal endothelial cells.' if label == "Peritoneal Epithelial Cells" else "No marker threshold or new assignment applied; genes listed are examined markers, not required-positive rules",
            })
        group_lowres = {g: set(low[(high == g) & selected]) for g in groups}
        summary = {
            "tissue": tissue, "path": str(cfg["path"]),
            "total_cells_in_source": len(high), "figure_cells_donors_1_2": int(selected.sum()),
            "figure_donor_counts": {b: int(np.sum(batch == b)) for b in ["3", "4"]},
            "total_genes_in_source": len(genes), "expression_layer": "layers/logcounts",
            "log1p_reference": cfg["reference"], "validated_selected_stored_entries": n_compared,
            "max_absolute_log1p_validation_error": max_error,
            "reference_log1p_mismatch_entries": n_reference_mismatches,
            "reference_diagnostic_used_for_plot": False,
            "expression_processing": "Stored finite nonnegative log-transformed expression in layers/logcounts; no new normalization. Original source code assigns SCT@data to this layer. Current X/counts may have been mutated downstream, so their relationship is recorded separately rather than used for plotting.",
            "no_clustering_DE_or_reannotation": True,
        }
        return groups, gene_list, means, fractions, n_cells, group_lowres, metadata_rows, summary

def render_panel(tissue, cfg, panel, data, plotted_rows=None):
    groups, gene_list, means, fractions, n_cells, group_lowres, _, _ = data
    tag, stem, title, markers = panel
    immune = stem.endswith("immune")
    selected_groups = [g for g in groups if bool(group_lowres[g] & cfg["immune_lowres"]) == immune]
    rows = [groups.index(g) for g in selected_groups]
    cols = [gene_list.index(g) for g in markers]
    panel_mean = means[np.ix_(rows, cols)]
    panel_frac = fractions[np.ix_(rows, cols)]
    maximum = panel_mean.max(axis=0)
    relative = np.divide(panel_mean, maximum, out=np.zeros_like(panel_mean), where=maximum > 0)
    height = max(4.45, 2.95 + 0.235 * len(rows))
    width = 8.2
    fig = plt.figure(figsize=(width, height))
    # Dimensions in inches preserve readable type across separately sized panels.
    left, right, bottom, top = 2.65, 0.20, 2.10, 0.70
    ax = fig.add_axes([left / width, bottom / height,
                       (width-left-right) / width, (height-bottom-top) / height])
    x, y = np.meshgrid(np.arange(len(markers)), np.arange(len(rows)))
    column_pitch_pt = (width-left-right) * 72 / len(markers)
    max_area = min(96.0, (column_pitch_pt * 0.80) ** 2)
    ax.scatter(x.ravel(), y.ravel(), s=(panel_frac * max_area).ravel(),
               c=relative.ravel(), cmap="Blues", norm=colors.Normalize(0, 1),
               linewidths=0.15, edgecolors="#31516c")
    ax.set_xlim(-0.55, len(markers)-0.45)
    ax.set_ylim(len(rows)-0.45, -0.55)
    ax.set_xticks(np.arange(len(markers)), markers, rotation=90, fontsize=9.2, fontstyle="italic")
    labels = [f"{DISPLAY.get(g,g)} (n={int(n_cells[groups.index(g)]):,})" for g in selected_groups]
    ax.set_yticks(np.arange(len(rows)), labels, fontsize=9.2)
    ax.tick_params(axis="both", length=0, pad=4)
    ax.set_xticks(np.arange(-0.5, len(markers), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="#e8edf0", linewidth=0.5)
    ax.tick_params(which="minor", length=0)
    ax.set_axisbelow(True)
    for side in ax.spines.values():
        side.set_color("#687580")
        side.set_linewidth(0.6)
    fig.text(0.025, 1-0.24/height, tag, fontsize=17, fontweight="bold", va="top")
    fig.text(0.5, 1-0.24/height, title, fontsize=10.4, fontweight="bold", va="top", ha="center", fontfamily="Arial")
    # Legends sit below the gene labels, so they do not squeeze the matrix.
    legend_ax = fig.add_axes([left/width, 0.08, 2.55/width, 0.32/height])
    legend_ax.axis("off")
    legend_ax.set_xlim(0, 5.8); legend_ax.set_ylim(-0.85, 1.0)
    for i, frac in enumerate([0.10, 0.25, 0.50, 0.75, 1.0]):
        legend_ax.scatter(i+0.5, 0, s=max_area*frac, c="#6b879e", linewidths=0.2, edgecolors="#31516c")
        legend_ax.text(i+0.5, -0.65, f"{int(frac*100)}", fontsize=8.5, ha="center", va="top")
    fig.text(left/width, 0.08+0.43/height, "Cells with expression > 0 (%)", fontsize=9.2)
    color_ax = fig.add_axes([(left+3.0)/width, 0.08, 2.15/width, 0.105/height])
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(0,1), cmap="Blues"),
                       cax=color_ax, orientation="horizontal", ticks=[0,0.5,1])
    bar.ax.tick_params(labelsize=8.5, length=2, pad=2)
    bar.outline.set_linewidth(0.5)
    if bar.solids is not None:
        bar.solids.set_rasterized(False)
    fig.text((left+3.0)/width, 0.08+0.43/height,
             "Relative mean expression", fontsize=9.2)
    fig.text((left+3.0)/width, 0.08+0.26/height,
             "Mean / gene maximum", fontsize=8.3, color="#45515c")
    base = RESULTS / PANEL_FILENAMES[tag]
    fig.savefig(base.with_suffix(".png"), dpi=600)
    plt.close(fig)
    quantitative = []
    for ri,g in enumerate(selected_groups):
        for ci,gene in enumerate(markers):
            quantitative.append({"panel":tag,"tissue":tissue,"source_highres_label":g,
                                 "gene":gene,"n_cells":int(n_cells[groups.index(g)]),
                                 "mean_stored_log_expression":float(panel_mean[ri,ci]),
                                 "fraction_expression_gt_0":float(panel_frac[ri,ci]),
                                 "relative_mean_for_color":float(relative[ri,ci])})
    if plotted_rows is not None:
        plotted_rows.extend(quantitative)
    return {"panel":tag,"stem":stem,"figure_file":"output/" + PANEL_FILENAMES[tag] + ".png","n_cells":int(n_cells[rows].sum()),
            "n_populations":len(rows),"n_markers":len(markers),
            "markers":markers,"source_labels":selected_groups,"width_inches":width,
            "height_inches":height,"resolution_dpi":600,"font":"Arial",
            "color_encoding":"mean stored log expression / maximum population mean for that gene within this panel",
            "size_encoding":"circle area proportional to fraction with stored log expression > 0"}

def render_saved_panels():
    """Render the saved rectangular plotting tables without changing their values."""
    panels = []
    with PLOTTING_TABLE.open(encoding="utf-8-sig", newline="") as handle:
        all_records = list(csv.DictReader(handle))
    for tissue, cfg in CONFIG.items():
        for panel in cfg["panels"]:
            tag, stem, _, markers = panel
            records = [row for row in all_records if row["panel"] == tag and row["tissue"] == tissue]
            groups = list(dict.fromkeys(row["source_highres_label"] for row in records))
            indexed = {(row["source_highres_label"], row["gene"]): row for row in records}
            assert len(indexed) == len(records) == len(groups) * len(markers)
            means = np.array([[float(indexed[g, marker]["mean_stored_log_expression"]) for marker in markers] for g in groups])
            fractions = np.array([[float(indexed[g, marker]["fraction_expression_gt_0"]) for marker in markers] for g in groups])
            n_cells = np.array([int(indexed[g, markers[0]]["n_cells"]) for g in groups])
            group_lowres = {g: cfg["immune_lowres"] if stem.endswith("immune") else set() for g in groups}
            data = (groups, markers, means, fractions, n_cells, group_lowres, [], {})
            panels.append(render_panel(tissue, cfg, panel, data))
    report_path = AUDIT / "annotation_sources_and_figure_provenance.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    report.update({"formats": ["png"], "font": "Arial", "title_alignment": "center", "subtitles": False, "display_label_corrections": {"Peritoneal Epithelial Cells": "Portal endothelial cells"}})
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print("Redrew four annotation PNGs from saved plotting values.")


def main():
    if ARGS.plots_only:
        render_saved_panels()
        return
    all_rows, sources, panels, plotted_rows = [], [], [], []
    for tissue,cfg in CONFIG.items():
        data = aggregate(tissue,cfg)
        all_rows.extend(data[-2]); sources.append(data[-1])
        for panel in cfg["panels"]:
            panels.append(render_panel(tissue,cfg,panel,data,plotted_rows=plotted_rows))
    write_csv(RESULTS / "cell_type_annotation_evidence_by_tissue_and_donor.csv", all_rows)
    write_csv(PLOTTING_TABLE, plotted_rows)
    report = {"sources":sources,"panels":panels,
              "display_label_corrections":{"Peritoneal Epithelial Cells":"Portal endothelial cells"},
              "formats":["png"],"font":"Arial","title_alignment":"center","subtitles":False,
              "scope":"Descriptive redraw of existing annotations and saved expression only. No normalization, reannotation, clustering, DE, or inference.",
              "notes":["Only requested marker entries for donors1/2 were validated; unused invalid expression entries were not imputed or plotted.",
                       "Fine MAIT labels are retained from these source objects, without creating the later CD8aa/CD8ab liver split.",
                       "Marker evidence table lists genes examined in original panels; it is not a required-positive classifier."]}
    (AUDIT / "annotation_sources_and_figure_provenance.json").write_text(json.dumps(report,indent=2),encoding="utf8")
    print(json.dumps(report,indent=2))

if __name__ == "__main__":
    main()
