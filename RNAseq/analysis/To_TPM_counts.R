base_dir <- "D:/RSEM_result"
samples <- list.files(base_dir)
files <- file.path(base_dir, samples)
names(files) <- samples
files

sample_name = strsplit(samples,'.genes.results')[1]
sample_df = read.csv(files[1],sep = '\t')
expected_counts = as.data.frame(sample_df[,5])
# TPM_counts = as.data.frame(sample_df[,6])
row.names(expected_counts) = sample_df[[1]]
colnames(expected_counts)=sample_name


for (quant in samples[2:length(samples)]){
  sample_name = strsplit(quant,'.genes.results')[1]
  sample_add_df = read.csv(file.path(base_dir, quant),sep = '\t')
  colnames(sample_add_df)[6] = sample_name
  expected_counts[sample_name[[1]]] = sample_add_df[6]
}

write.csv(expected_counts,'D:/Trapecar/expected_counts_STAR_RSEM.csv')