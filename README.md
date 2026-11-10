# Gut liver model

Analysis code for **Autologous Human Gut–Liver Microphysiological Systems Link Tissue-Resident Immune Architecture to Donor-Specific Challenge Responses**.
The associated [preprint](https://doi.org/10.1101/2025.05.05.652244) appeared under an earlier title.

## Repository layout

| Directory | Contents |
|---|---|
| `analysis/scrna` | Historical tissue preprocessing, annotation, marker figures and donor DE steps |
| `analysis/nichenet` | Historical preparation, corrected ligand ranking and expression figures |
| `analysis/cellphonedb` | Original inference notebooks and chord plots |
| `analysis/bulk_rnaseq` | QC, quantification, named-column RSEM collection, DE/enrichment |
| `analysis/deconvolution` | Mixture and study/public reference-profile preparation |
| `analysis/revision` | Reproducible response package for reviewer comments 6–10 |
| `data/inputs`, `data/derived` | Existing bulk matrices and generated deconvolution inputs |
| `references` | Adapter and gene-annotation resources |
| `results` | Existing exported results |
| `config`, `docs` | Local paths, workflow/provenance and original-to-new file map |

`docs/file_moves.tsv` accounts for all 114 original tracked files. Original versions remain in Git history at baseline `c21587a`. Large external H5ADs and receiver objects remain outside Git.

## Build the revision package

Copy `config/paths.example.json` to `config/paths.local.json` and fill the source paths. The local file is ignored by Git. A configured local file is already supplied in this working checkout.

```powershell
python analysis/revision/build_package.py --config config/paths.local.json --stage-only
python analysis/revision/build_package.py --config config/paths.local.json
```

Use `--steps nichenet`, `--steps annotation`, or `--steps cellphonedb` to run one component. `--install-deps` prepares local dependencies when required; it does not install into the global Python/R environments.

The package reuses the four original receiver DEG tables and original serialized receiver cohorts. It reconstructs detection backgrounds and corrects donor-specific NicheNet scoring/target exports. The main build does not refit DE, rerun normalization/clustering, or rerun CellPhoneDB inference. Annotation panels summarize saved labels and stored marker expression; CellPhoneDB count matrices are recovered exactly from baseline saved tables and checked against available originals.

Outputs include text responses, separate scientific legends, PNG figures, and supporting tables. Figures use Arial text, centered titles and no subtitles. Final deliverables are in each analysis folder's `output/`; technical checks are in `audit/`. Supplementary numbers in the response package are working identifiers until inserted into the manuscript.

The separate [liver MAIT refinement script](analysis/revision/annotation/refine_liver_mait.py) replays the notebook's restricted Leiden calls using an H5AD's saved graph. Use Python 3.9 with [its pinned requirements](analysis/revision/annotation/requirements_mait_refinement.txt). It records cluster counts and marker expression in `09_annotation/audit/`; `--export-proposed-labels` optionally exports candidate cell labels there. These labels require cluster-identity review before use: a rerun on the current H5AD does not reproduce the original notebook's second-pass cluster structure. The script leaves source H5ADs and final annotation figures unchanged.

## Historical analysis steps

Numbering is local to each workflow. In `analysis/scrna`, step 03 prepares donor objects, step 04 fits historical donor DE, and step 05 plots saved results. The revision uses recovered original objects/tables and therefore skips steps 03–04. The preprocessing notebooks retain original scientific cells and saved outputs, including historical local paths and non-linear execution state; they are not validated as fresh end-to-end runs. See `docs/REPRODUCIBILITY.md` before attempting a historical rerun.

The RSEM collector now reads named `expected_count` and `TPM` fields consistently and aligns by gene ID. Existing deposited matrices are preserved unchanged.

## Interpretation and licensing

Ligand activities are predictions, not unique causal drivers. Receiver DEGs used for prioritization provide consistency evidence rather than independent signaling validation. Two donors support descriptions of these individuals rather than population-level donor inference.

The original repository license is retained. Pinned upstream NicheNet source in `analysis/revision/nichenet/vendor` retains its GPL-3 license and notices; see that directory before redistributing third-party code.
