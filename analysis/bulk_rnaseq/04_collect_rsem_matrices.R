# Collect RSEM .genes.results into separate expected-count and TPM matrices.
# Replaces two identical historical scripts that mixed columns 5 and 6.
# Existing deposited matrices are never modified by this cleanup.
# Usage: Rscript analysis/bulk_rnaseq/04_collect_rsem_matrices.R INPUT_DIR OUTPUT_DIR [--overwrite]

collect_rsem_matrices <- function(input_dir, output_dir, overwrite = FALSE) {
  if (!dir.exists(input_dir)) stop("RSEM input directory does not exist: ", input_dir)
  files <- sort(list.files(input_dir, pattern = "\\.genes\\.results$", full.names = TRUE))
  if (!length(files)) stop("No .genes.results files found in ", input_dir)
  samples <- sub("\\.genes\\.results$", "", basename(files))
  if (anyDuplicated(samples)) stop("Sample names are not unique.")
  read_sample <- function(path) {
    tab <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    required <- c("gene_id", "expected_count", "TPM")
    if (!all(required %in% names(tab))) stop("Missing required RSEM columns in ", path)
    if (!nrow(tab) || anyNA(tab$gene_id) || any(!nzchar(tab$gene_id)) || anyDuplicated(tab$gene_id))
      stop("gene_id values must be non-empty and unique in ", path)
    for (field in c("expected_count", "TPM")) {
      if (!is.numeric(tab[[field]]) || any(!is.finite(tab[[field]])) || any(tab[[field]] < 0))
        stop("Invalid non-negative numeric ", field, " in ", path)
    }
    tab
  }
  tables <- lapply(files, read_sample)
  genes <- as.character(tables[[1]]$gene_id)
  for (i in seq_along(tables)) {
    if (!setequal(genes, as.character(tables[[i]]$gene_id)))
      stop("Gene sets differ between samples; refusing positional alignment: ", files[[i]])
  }
  collect_field <- function(field) {
    columns <- lapply(tables, function(tab) tab[[field]][match(genes, tab$gene_id)])
    result <- as.data.frame(columns, check.names = FALSE)
    names(result) <- samples
    rownames(result) <- genes
    result
  }
  expected_counts <- collect_field("expected_count")
  tpm <- collect_field("TPM")
  destinations <- file.path(output_dir, c("raw_counts_STAR_RSEM.csv", "TPM_counts_STAR_RSEM.csv", "rsem_input_files.tsv"))
  if (!overwrite && any(file.exists(destinations)))
    stop("Output already exists; choose a new output directory or explicitly pass --overwrite.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(expected_counts, destinations[[1]], row.names = TRUE)
  write.csv(tpm, destinations[[2]], row.names = TRUE)
  sources <- data.frame(sample = samples, source_file = normalizePath(files, winslash = "/"))
  write.table(sources, destinations[[3]], sep = "\t", row.names = FALSE, quote = TRUE)
  invisible(list(expected_counts = expected_counts, tpm = tpm, sources = sources))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2L || length(args) > 3L || (length(args) == 3L && args[[3]] != "--overwrite"))
    stop("Usage: Rscript 04_collect_rsem_matrices.R INPUT_DIR OUTPUT_DIR [--overwrite]")
  collect_rsem_matrices(args[[1]], args[[2]], overwrite = "--overwrite" %in% args)
}
