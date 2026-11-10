# Bulk RNA sequencing workflow

1. `analysis/bulk_rnaseq/01_fastqc.sh`: inspect raw read quality.
2. `02_trim_adapters.sh`: adapter trimming with the reference under `references/adapters`.
3. `03_quantify_star_rsem.slurm`: historical cluster quantification. Adapt scheduler and reference locations to the available cluster. Alternative STAR/Salmon examples are separate under `alternatives`.
4. `04_collect_rsem_matrices.R INPUT_DIR OUTPUT_DIR`: collect every sample's named expected_count and TPM columns into separate matrices, aligning gene IDs. Existing outputs are not overwritten unless explicitly requested.
5. `05_differential_expression_and_enrichment.R`: historical comparisons and enrichment; requires the original packages and comparison inputs. Run from the repository root or set GUT_LIVER_REPO_ROOT. Outputs route to results/bulk_rnaseq or GUT_LIVER_BULK_OUTPUT.

The existing raw-count, TPM and metadata inputs are under data/inputs/bulk_rnaseq. They are unchanged by this revision. No bulk reanalysis was required for comments 6–10.
