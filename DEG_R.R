library(monocle3)
library(Seurat)

setwd('E:/AAA_Labwork/Tcell tissues/v2/')
srat = readRDS('gutT')

cds <- as.cell_data_set(srat)
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(srat[["RNA"]])
  
# Run DE analysis
gene_fits <- fit_models(cds, model_formula_str = "~DonorID",cores = 5)
fit_coefs <- coefficient_table(gene_fits)

# Adjust condition in the filter
terms <- fit_coefs %>% filter(term == 'DonorIDDonor AJG2309')
terms = terms %>% select(gene_short_name, term, q_value, normalized_effect)

write.csv(terms,'gutT_DEG.csv')


# gutMac ------------------------------------------------------------------


srat = readRDS('gutMac')

cds <- as.cell_data_set(srat)
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(srat[["RNA"]])

# Run DE analysis
gene_fits <- fit_models(cds, model_formula_str = "~DonorID",cores = 5)
fit_coefs <- coefficient_table(gene_fits)

# Adjust condition in the filter
terms <- fit_coefs %>% filter(term == 'DonorIDDonor AJG2309')
terms = terms %>% select(gene_short_name, term, q_value, normalized_effect)

write.csv(terms,'gutMac_DEG.csv')

# liverT ------------------------------------------------------------------
srat = readRDS('liverT')

cds <- as.cell_data_set(srat)
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(srat[["RNA"]])

# Run DE analysis
gene_fits <- fit_models(cds, model_formula_str = "~DonorID",cores = 5)
fit_coefs <- coefficient_table(gene_fits)

# Adjust condition in the filter
terms <- fit_coefs %>% filter(term == 'DonorIDDonor AJG2309')
terms = terms %>% select(gene_short_name, term, q_value, normalized_effect)

write.csv(terms,'liverT_DEG.csv')


# liverMac ------------------------------------------------------------------

srat = readRDS('liverMac')

cds <- as.cell_data_set(srat)
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(srat[["RNA"]])

# Run DE analysis
gene_fits <- fit_models(cds, model_formula_str = "~DonorID",cores = 5)
fit_coefs <- coefficient_table(gene_fits)

# Adjust condition in the filter
terms <- fit_coefs %>% filter(term == 'DonorIDDonor AJG2309')
terms = terms %>% select(gene_short_name, term, q_value, normalized_effect)

write.csv(terms,'liverMac_DEG.csv')