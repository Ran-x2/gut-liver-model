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
  d4_bg_filename = paste0('3_',file,'_receptors.csv')
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
  
  top_15_ligands = ligand_activities$test_ligand[1:15]
  target_gene_union <- as.data.frame( matrix(NA, nrow = 5, ncol = 15))
  colnames(target_gene_union) <- top_15_ligands
  for (ligand in top_15_ligands) {
      top_5_targets =ligand_target_matrix[,ligand][order(ligand_target_matrix[,ligand],decreasing = TRUE) ][1:5]
      target_gene_union[ligand] = names(top_5_targets)
    }
  all_targets <- unique(as.vector(as.matrix(target_gene_union)))
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
  
  top_15_ligands = ligand_activities$test_ligand[1:15]
  target_gene_union <- as.data.frame( matrix(NA, nrow = 5, ncol = 15))
  colnames(target_gene_union) <- top_15_ligands
  for (ligand in top_15_ligands) {
    top_5_targets =ligand_target_matrix[,ligand][order(ligand_target_matrix[,ligand],decreasing = TRUE) ][1:5]
    target_gene_union[ligand] = names(top_5_targets)
  }
  all_targets <- unique(as.vector(as.matrix(target_gene_union)))
  write.table(all_targets, paste0('4_', file,'_nichenet_targets.csv'), row.names = FALSE) 
}
