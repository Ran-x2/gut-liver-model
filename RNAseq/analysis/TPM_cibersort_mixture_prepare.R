install.packages("e1071")
BiocManager::install("preprocessCore")
BiocManager::install("limma")
devtools::install_github("BaderLab/HumanLiver")
library(HumanLiver)
viewHumanLiver()
library(dplyr)
library(tibble)
setwd('C:/Users/16220/Documents/GitHub/gut-liver-model/RNAseq/')

metadata = read.csv('analysis/data/metadata.csv')
metadata$Study.Design <- sub("-[0-9]+$", "", metadata$Study.Design)

metadata$Donor <- factor(metadata$Donor)
metadata$Cell.Type <- factor(metadata$Cell.Type)
metadata$Study.Design.Batch <- factor(metadata$Study.Design.Batch)
metadata$Study.Design <- factor(metadata$Study.Design)
metadata$Treatment <- factor(metadata$Treatment)
metadata$Cell.Type.Study.Design = factor(paste0(metadata$Cell.Type,'__',metadata$Study.Design))
metadata$Treatment <- relevel(metadata$Treatment, ref = "NO TREATMENT")
metadata$Study.Design <- relevel(metadata$Study.Design, ref = "ISOLATION (ISO) PURE")

raw_counts = read.csv('analysis/data/TPM_counts_STAR_RSEM.csv',row.names = 1)
raw_counts = raw_counts[,metadata$Sample.Name]
unique_celltypes <- unique(metadata$Cell.Type)
unique_donor <- unique(metadata$Donor)

ref = read.csv('analysis/features_10x_2020A.tsv',sep = '\t',header = 1)
lookup = ref[,2]
names(lookup) = ref[,1]

# Find duplicated values in lookup
dup_vals <- lookup[duplicated(lookup)]
# For each duplicated value, append the name (key) to the value with "_"
for (key in names(dup_vals)) {
  val <- lookup[key]
  lookup[key] <- paste(val, key, sep = "_")
}
lookup[duplicated(lookup)]

split_names <- strsplit(row.names(raw_counts), "_")
ensg_id = sapply(split_names, function(x) x[1])
gene_names = sapply(split_names, function(x) x[2])

choose_id <- function(x, y) {
  if (length(x) == 2 && x[1] %in% names(y)) {
    return(y[[x[[1]]]])
  } else if (length(x) == 2) {
    return(paste0(c(x[[2]], x[[1]]), collapse = "_"))
  } else {
    return(x[[1]])
  }
}
mapped_gene_names <- mapply(choose_id, split_names, MoreArgs = list(y = lookup))

mapped_gene_names[duplicated(mapped_gene_names)]
row.names(raw_counts) <- mapped_gene_names
write.csv(raw_counts,'analysis/data/symbol_mappedTPM_counts.csv')


for (celltype in unique_celltypes) {
subset_metadata <- metadata %>%
  filter(Cell.Type == celltype)
subset_counts <- raw_counts[, subset_metadata$Sample.Name]
subset_filename = paste('deconvolution/input_for_cibersortx/',celltype,'_TPM.tsv',sep ='')
subset_filename = gsub(" ", "_", subset_filename)
subset_filename = gsub("\\+", 'pos', subset_filename)
print(subset_filename)
subset_counts = add_column(subset_counts, GeneSymbol = row.names(subset_counts), .before = 1)
write.table(subset_counts,subset_filename, sep = "\t",row.names = FALSE,quote=FALSE)
}

for (donor in unique_donor) {
for (celltype in unique_celltypes) {
subset_metadata <- metadata %>%
  filter(Donor == donor, Cell.Type == celltype)
subset_counts <- raw_counts[, subset_metadata$Sample.Name]
subset_filename = paste('deconvolution/input_for_cibersortx/by_donor/',donor,'_',celltype,'_TPM.tsv',sep ='')
subset_filename = gsub(" ", "_", subset_filename)
subset_filename = gsub("\\+", 'pos', subset_filename)
print(subset_filename)
subset_counts = add_column(subset_counts, GeneSymbol = row.names(subset_counts), .before = 1)
write.table(subset_counts,subset_filename, sep = "\t",row.names = FALSE,quote=FALSE)
}
}
