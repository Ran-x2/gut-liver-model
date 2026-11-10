"""Shared names and display labels for the four saved CellPhoneDB count cases."""
from pathlib import Path

CASES = [
    ("gut_3", "A", "gut", 1, "AJD3280"),
    ("gut_4", "B", "gut", 2, "AJG2309"),
    ("liver_3", "C", "liver", 1, "AJD3280"),
    ("liver_4", "D", "liver", 2, "AJG2309"),
]
FIGURE_STEM = "Figure_R10_cellphonedb_interaction_counts_by_tissue_and_donor"
PRESENCE_FILE = "cell_population_display_names_and_presence.csv"
COUNT_FIGURE_LEGEND = """Figure R10. Quantitative CellPhoneDB interaction counts across donor-matched gut and liver cell populations.

(A) AJD3280 gut, (B) AJG2309 gut, (C) AJD3280 liver, and (D) AJG2309 liver. Each cell shows the number of retained significant molecular-interaction entries between the indicated cell populations, summed over both ordered partner assignments (A|B and B|A). The matrices are therefore symmetric; an interaction significant for both assignments can contribute to both constituent counts. The same linear count scale is used for all four panels, with exact counts printed in the cells. Gray diagonal cells denote excluded self-pairs. Other gray rows and columns denote populations absent from the saved donor-specific metadata; in panel D these are hepatocytes and glial-like cells. Tested population pairs with no retained interaction are shown as numeric zero rather than gray.

The original CellPhoneDB statistical analysis used database release v4.1.0, 1,000 permutations, an expression-fraction cutoff of 0.30, nominal permutation P < 0.05, no subsampling, seed 42, and six threads. The original plotting workflow first transformed significant mean-expression values as log(1 + value), then retained entries strictly above the 5th percentile of positive transformed values within each donor/tissue analysis. Retained entries were counted per ordered cell-type pair and the two assignments were summed. The heatmaps show these stored matrices before the additional count < 5 suppression used when drawing the original R chord diagrams. These counts summarize candidate transcript-based interactions and do not measure physical signaling, calibrated communication strength, or statistically significant differences between donors.

DC, dendritic cells; FRC, fibroblastic reticular cells; ILC, innate lymphoid cells; NK, natural killer cells; intestinal epi., intestinal epithelial cells; enteric glia, enteric glial cells; stellate, stellate cells; PEC, portal endothelial cells; LSEC, liver sinusoidal endothelial cells. Plasma denotes plasma cells.
"""
DISPLAY = {
    "T Cells": "T cells", "NK": "NK", "ILC": "ILC", "B Cells": "B cells",
    "Plasma Cells": "Plasma", "Monocytes": "Monocytes", "Macrophages": "Macrophages",
    "Dendritic Cells": "DC", "Endothelial Cells": "Endothelial", "Telocytes": "Telocytes",
    "Fibroblastic Reticular Cells": "FRC", "Fibroblast": "Fibroblasts",
    "Smooth Muscle Cells": "Smooth muscle", "Intestinal Epithelial Cells": "Intestinal epi.",
    "Enteric Glial Cells": "Enteric glia", "Erthyroid": "Erythroid",
    "Mast Progenitor": "Mast progenitor", "PEC": "PEC", "LSEC": "LSEC",
    "Stellate Cells": "Stellate", "Cholangiocytes": "Cholangiocytes",
    "Hepatocytes": "Hepatocytes", "Glial-like Cells": "Glial-like",
}


def count_table_path(folder: Path, key: str) -> Path:
    case = next(case for case in CASES if case[0] == key)
    _, _, tissue, donor_number, donor_id = case
    return folder / "output" / f"{tissue}_donor_{donor_number}_{donor_id}_cellphonedb_interaction_counts.csv"
