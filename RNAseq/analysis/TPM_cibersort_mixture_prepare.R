install.packages("e1071")
BiocManager::install("preprocessCore")
BiocManager::install("limma")
devtools::install_github("BaderLab/HumanLiver")
library(HumanLiver)
viewHumanLiver()
library(dplyr)

metadata = read.csv('D:/Trapecar/metadata.csv')
metadata$Study.Design <- sub("-[0-9]+$", "", metadata$Study.Design)

metadata$Donor <- factor(metadata$Donor)
metadata$Cell.Type <- factor(metadata$Cell.Type)
metadata$Study.Design.Batch <- factor(metadata$Study.Design.Batch)
metadata$Study.Design <- factor(metadata$Study.Design)
metadata$Treatment <- factor(metadata$Treatment)
metadata$Cell.Type.Study.Design = factor(paste0(metadata$Cell.Type,'__',metadata$Study.Design))
metadata$Treatment <- relevel(metadata$Treatment, ref = "NO TREATMENT")
metadata$Study.Design <- relevel(metadata$Study.Design, ref = "ISOLATION (ISO) PURE")

raw_counts = read.csv('D:/Trapecar/TPM_counts_STAR_RSEM.csv',row.names = 1)
raw_counts = raw_counts[,metadata$Sample.Name]
unique_celltypes <- unique(metadata$Cell.Type)
unique_donor <- unique(metadata$Donor)

choose_id <- function(x) {
  if (length(x) == 2) {
    return(x[2])
  } else {
    return(x[1])
  }
}
split_names <- strsplit(row.names(raw_counts), "_")
string_array <- sapply(split_names, choose_id)
raw_counts['Gene.ID'] = string_array
summed_counts <- raw_counts %>%
  group_by(Gene.ID) %>%
  summarise(across(everything(), sum, na.rm = TRUE))
write.csv(summed_counts,'D:/Trapecar/summed_TPM_counts.csv')


for (donor in unique_donor) {
for (celltype in unique_celltypes) {
subset_metadata <- metadata %>%
  filter(Donor == donor, Cell.Type == celltype)

subset_counts <- summed_counts[, c('Gene.ID',subset_metadata$Sample.Name)]
subset_filename = paste('D:/Trapecar/',donor,celltype,'.tsv',sep ='')
subset_filename = gsub(" ", "_", subset_filename)
subset_filename = gsub("\\+", 'pos', subset_filename)
print(subset_filename)
write.table(subset_counts,subset_filename, sep = "\t",row.names = FALSE,quote=FALSE)
}
}
##############################################################
library(data.table)
liver_data = read.csv('D:/Trapecar/GSE115469_Data.csv',row.names = 1)
annotation_liver = read.csv('D:/Trapecar/GSE115469_CellClusterType.txt', sep = '\t')
colnames(liver_data) ==annotation_liver$CellName
liver_data = data.frame(liver_data, check.names = FALSE)
colnames(liver_data) = annotation_liver$CellType
names(liver_data) <- sub("\\.\\d+$", "", names(liver_data))
colnames(liver_data)
liver_data_hep = liver_data[,startsWith(colnames(liver_data),'Hep')]
write.table(liver_data_hep,'D:/Trapecar/deconvolute/hep_sig.txt', sep = "\t",row.names = TRUE,quote=FALSE)
###############################################################
