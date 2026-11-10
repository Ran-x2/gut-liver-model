# Corrected historical source; configure the personal/reference paths below explicitly.
# Donor-4 background, sorted top-15 ligands, and DEG-restricted reported targets are fixed.
# Historical results are retained; this script was not rerun during repository cleanup.

top_receiver_deg_targets <- function(ligands, receiver_degs, target_matrix, n = 5L) {
  genes <- intersect(unique(receiver_degs), rownames(target_matrix))
  ligands <- intersect(ligands, colnames(target_matrix))
  if (!length(genes) || !length(ligands)) return(character())
  unique(unlist(lapply(ligands, function(ligand) {
    weights <- target_matrix[genes, ligand, drop = TRUE]
    names(head(sort(weights, decreasing = TRUE), n))
  }), use.names = FALSE))
}

library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggplot2)
library(car)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
# BiocManager::install("GenomicRanges")
# BiocManager::install("rtracklayer")
library(GenomicRanges)
library(rtracklayer)

setwd('E:/AAA_Labwork/Tcell tissues/v2')
lr_network <- readRDS('F:/CRC/Nichenet/lr_network_human_21122021.rds')
ligand_target_matrix <- readRDS('F:/CRC/Nichenet/ligand_target_matrix_nsga2r_final.rds')
weighted_networks <- readRDS('F:/CRC/Nichenet/weighted_networks_nsga2r_final.rds')
gr_network = readRDS(url("https://zenodo.org/record/7074291/files/gr_network_human_21122021.rds"))
sig_network = readRDS(url("https://zenodo.org/record/7074291/files/signaling_network_human_21122021.rds"))

filenames = c('gutT','gutMac','liverT','liverMac')
for (file in filenames) {
  
  d3_bg_filename = paste0('3_',file,'_receptors.csv')
  d4_bg_filename = paste0('4_',file,'_receptors.csv')
  deg_filename = paste0(file,'_DEG.csv') 
    
  # ligands influence d3 compare to d4
  background_expressed_genes<- read.table(d3_bg_filename, sep = ',', header = 1, row.names = 1)
  potential_ligands <- lr_network %>% filter(to %in% background_expressed_genes$X0) %>% pull(from) %>% unique()
  
  DE_table_receiver = read.table(deg_filename, sep = ',', header = 1, row.names = 1)
  geneset_oi <- DE_table_receiver %>% filter(q_value <= 0.05 & normalized_effect <= -0.5) %>% pull(gene_short_name)
  geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                                 background_expressed_genes = background_expressed_genes$X0,
                                                 ligand_target_matrix = ligand_target_matrix,
                                                 potential_ligands = potential_ligands)
  
  sorted_table <- ligand_activities[order(ligand_activities$auroc,decreasing = TRUE), ]
  write.table(sorted_table, paste0('3_', file,'_nichenet.csv'), row.names = FALSE)
  
  top_15_ligands = head(sorted_table$test_ligand, 15L)
  # Report supported targets only within the direction-specific receiver DEG set.
  all_targets <- top_receiver_deg_targets(top_15_ligands, geneset_oi, ligand_target_matrix)
  write.table(all_targets, paste0('3_', file,'_nichenet_targets.csv'), row.names = FALSE)  
  
  #ligands influence d4 compare to d3
  background_expressed_genes<- read.table(d4_bg_filename, sep = ',', header = 1, row.names = 1)
  potential_ligands <- lr_network %>% filter(to %in% background_expressed_genes$X0) %>% pull(from) %>% unique()
  
  DE_table_receiver = read.table(deg_filename, sep = ',', header = 1, row.names = 1)
  geneset_oi <- DE_table_receiver %>% filter(q_value <= 0.05 & normalized_effect >= 0.5) %>% pull(gene_short_name)
  geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                                 background_expressed_genes = background_expressed_genes$X0,
                                                 ligand_target_matrix = ligand_target_matrix,
                                                 potential_ligands = potential_ligands)
  
  sorted_table <- ligand_activities[order(ligand_activities$auroc,decreasing = TRUE), ]
  write.table(sorted_table, paste0('4_', file,'_nichenet.csv'), row.names = FALSE)
  
  top_15_ligands = head(sorted_table$test_ligand, 15L)
  # Report supported targets only within the direction-specific receiver DEG set.
  all_targets <- top_receiver_deg_targets(top_15_ligands, geneset_oi, ligand_target_matrix)
  write.table(all_targets, paste0('4_', file,'_nichenet_targets.csv'), row.names = FALSE) 
}
