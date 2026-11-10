# Reproducibility and scope

Original scientific inputs are preserved. The file-move map records old/new paths; historical notebooks retain their original code/output cells with provenance notes. New revision scripts use explicit configured inputs and export source hashes, count tables, metric/source versions and figure legends.

The original receiver objects establish the exact cell cohorts. The gut T input excludes NK-labelled cells; a literal replay of an incomplete broad-label recode would not reproduce it. Use recovered receiver objects for reproduction. Gene sets retain saved `q_value <= 0.05` and directional `normalized_effect` thresholds, followed by matching to the recorded model.

Scientific corrections are limited to donor-4 backgrounds, explicit ligand ranking, receiver-DEG-restricted target links, and comparable presentation. The original DEG results are not refitted. Cached scores are accepted only after candidate/background and numerical concordance checks. Undefined Pearson correlations are compared using matching NA patterns.

The manuscript's `log2(coefficient/intercept)` formula does not describe Monocle3 normalized_effect. The original adjustment method should not be asserted as BH solely from a column named q_value; retain saved values and consult the pinned source note. Marker terminology and initial/refined clustering discrepancies are recorded in the output package.

Historical preprocessing and bulk notebooks require their original data/environment and may retain original local paths. They were not scientifically rerun during cleanup. Bulk design formulas and deposited numeric matrices remain unchanged. The corrected RSEM collector is verified with a two-sample gene-reordering fixture; this code fix does not establish that deposited matrices were generated incorrectly.

The CellPhoneDB count figure is a descriptive visualization of retained predicted interaction entries, not calibrated signaling strength or a statistical donor comparison. Source table recovery, original-file checks and filtering/symmetry rules are provided with the figure.
