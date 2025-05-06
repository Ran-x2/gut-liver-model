# Load Library ------------------------------------------------------
organism = "org.Hs.eg.db"
library(qvalue)
library(enrichR)
library(organism, character.only = TRUE)
library(enrichplot)
library(clusterProfiler)
library(msigdbr)
library(ggplot2)
library(car)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
library(ggrepel)
library(pheatmap)
library(devtools)
install_github("stephens999/ashr")
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
dbs <- c("Reactome_2022","WikiPathway_2023_Human","KEGG_2021_Human","GO_Biological_Process_2023","GO_Molecular_Function_2023","GO_Cellular_Component_2023")
setwd('E:/AAA_Labwork/Tcell model RNAseq')
# Read Data ------------------------------------------------------
metadata = read.csv('metadata.csv')
metadata$Study.Design <- sub("-[0-9]+$", "", metadata$Study.Design)

metadata$Cell.Type <- factor(metadata$Cell.Type)
metadata$Treatment <- factor(metadata$Treatment)
metadata$Treatment <- relevel(metadata$Treatment, ref = "NO TREATMENT")
metadata$Study.Design <- factor(metadata$Study.Design)
metadata$Study.Design <- relevel(metadata$Study.Design, ref = "ISOLATION (ISO) PURE")

raw_counts = read.csv('raw_counts_STAR_RSEM.csv',row.names = 1)
raw_counts = raw_counts[,metadata$Sample.Name]
raw_counts = round(raw_counts)

choose_id <- function(x) if (length(x) == 2) x[2] else x[1]
easy_search <- function(x) if (length(x) == 2) x[2] else ""
choose_ENSG_id <- function(x) x[1]
check_strings <- function(name, x) grepl(name, x, fixed = TRUE)
choose_ENSG_id <- function(x) {
  return(x[1])
}
evaldiv = function(x) eval(parse(text = x ))

# Compare Treatments ------------------------------------------------------
dir.create('donor_celltype_compareTREAT')
setwd('donor_celltype_compareTREAT')
studydesign <- 'INTERACTION (INT) GUT-LIVER-IMMUNE CELLS'
subset_metadata <- metadata %>%
  filter(Study.Design == studydesign)
unique_donor <- unique(metadata$Donor)
unique_celltypes <- unique(subset_metadata$Cell.Type)
unique_treatments <- unique(subset_metadata$Treatment)
unique_design = unique(subset_metadata$Study.Design)
for (donor in unique_donor) {
  for (celltype in unique_celltypes) {
    organism = "org.Hs.eg.db"
    subset_metadata <- metadata %>%
      filter(Donor == donor, Cell.Type == celltype, Study.Design == studydesign)
    # we will go check here: this is the step where the samples used for in-group comparison is selected.
    subset_counts <- raw_counts[, subset_metadata$Sample.Name]
    print(subset_metadata)
    dds <- DESeqDataSetFromMatrix(countData = subset_counts, colData = subset_metadata, design = ~ Treatment)
    dds <- dds[rowSums(counts(dds)) >= 10, ]
    dds <- DESeq(dds)
    #dds <- DESeq(dds, minReplicatesForReplace=Inf)
    saveRDS(dds, paste0(donor, celltype, "_grouped.rds"))

    normalized_counts <- counts(dds, normalized = TRUE)
    genes_to_plot <- c('CYP2A7','CYP3A7','CYP2A6','ACSS2','SEC16B','GSTA2','BCHE','EPCAM','ASCL2','CLDN10','KRT20','MUC2','SLC26A2',
                       'CA1','HES1','SPIB','CHGA','CD27','NCAM1','CD86','CD1C','HLA.DPB1','IGHM',
                       'IGHA1','CD19','CD20','CD3G','CD3D','CD3E','TRDC','CD8A','KLRB1','TRAV1.2','KIT','PTPRC')
    split_names <- strsplit(row.names(normalized_counts), "_")
    string_array <- sapply(split_names, choose_id)
    rownames(normalized_counts) = string_array
    # Subset normalized counts for these genes
    normalized_subset <- normalized_counts[rownames(normalized_counts) %in% genes_to_plot, ]
    #breaks <- seq(0, 2000, length.out = 21)
    # Plot the heatmap
    dds_heat = pheatmap(normalized_subset,
                        #breaks = breaks,
                        #color = colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1),
                        scale = "none", # Normalize by row (gene-wise) to better visualize differences
                        cluster_rows = TRUE,
                        cluster_cols = TRUE,
                        show_rownames = TRUE,
                        show_colnames = TRUE)
    ggsave(paste0(donor,'_',celltype,'_heatmap.png'), plot = dds_heat, width = 10, height = 10)

    res_names <- resultsNames(dds)

    for (res_name in res_names) {
      if (res_name != "Intercept") {
        res = results(dds,name = res_name)
        res_name_sanitized <- gsub("c.", "_", res_name)
    #     res.shrink <- lfcShrink(dds = dds, coef = res_name, type = "ashr")
    #     res_for_vp <- res.shrink
        res_for_vp = res
    #     split_names <- strsplit(row.names(res.shrink), "_")
        split_names <- strsplit(row.names(res_for_vp), "_")

        ENSGID <- sapply(split_names, choose_ENSG_id)
        string_array <- sapply(split_names, choose_id)
        easy_search_array <- sapply(split_names, easy_search)
    #     res.shrink$GeneSymbol <- easy_search_array
    #     res.shrink$ENSGID <- ENSGID
        res$GeneSymbol <- easy_search_array
        res$ENSGID <- ENSGID

        row.names(res_for_vp) <- string_array
        volcano_labels <- ifelse(res_for_vp$padj <= 0.05 & abs(res_for_vp$log2FoldChange) > 1, rownames(res_for_vp), '')
        volcano_labels[is.na(volcano_labels)] <- ""
        volcano_labels <- ifelse(!startsWith(string_array, 'ENSG'), volcano_labels, '')

        volcanoplot <- EnhancedVolcano(
          res_for_vp,
          lab = '',  # No labels here to avoid duplication
          x = 'log2FoldChange',
          y = 'padj',
          ylab = bquote(~-Log[10] ~ italic(Padj)),
          pCutoff = 0.05,
          FCcutoff = 1,
          pointSize = 1.0,
          labSize = 0.5,
          drawConnectors = TRUE,
          widthConnectors = 0.5,
          max.overlaps = 30
        )

        volcanoplot <- volcanoplot + geom_text_repel(
          aes(label = volcano_labels),
          size = 2,  # Adjust size as needed
          segment.color = 'grey50',
          color = 'black',
          bg.color = 'cyan',  # Background color for the outline
          bg.r = 0.1,  # Radius for the outline
          max.overlaps = 70,  # Allow more overlaps
          min.segment.length = 0.5  # Minimum segment length for label connectors
        )

    volcano_filename <- paste0( donor, celltype, "_", res_name_sanitized, "volcanoplot.png")
    ggsave(volcano_filename, plot = volcanoplot, width = 5, height = 7)

    # resOrdered <- res.shrink[order(res.shrink$padj), ]
    resOrdered <- res[order(res$padj), ]
    
    filename <- paste( donor, celltype, res_name_sanitized, "TreatmentCompare.csv", sep = "_")
    write.csv(as.data.frame(resOrdered), file = filename)
    
    significant_results <- subset(resOrdered, padj <= 0.05)
    significant_results<- significant_results[order(significant_results$log2FoldChange), ]
    
    filename <- paste( donor, celltype, res_name_sanitized, "TreatmentCompare_significantONLY_sortedFC.csv", sep = "_")
    write.csv(as.data.frame(significant_results), file = filename)
    
    gene_list <- significant_results$log2FoldChange
    names(gene_list) <- rownames(significant_results)
    
    gene_list <- sort(gene_list, decreasing = TRUE)
    ensg_list = gene_list
    
    names(gene_list) <- sapply(strsplit(names(gene_list), "_"), choose_id)
    gene_list <- gene_list[!startsWith(names(gene_list), 'ENSG')]
    
    names(ensg_list) <- sapply(strsplit(names(ensg_list), "_"),choose_ENSG_id)
    print(ensg_list)
    
    if (length(gene_list) == 0) {
      message("No significant genes found for ", donor,  celltype, " ", res_name_sanitized)
      next
    }
    
    try({
    
      uplist = names(gene_list[gene_list>0])
      downlist = names(gene_list[gene_list<0])
      UPenriched <- enrichr(uplist, dbs)
      Sys.sleep(5)
      DOWNenriched = enrichr(downlist, dbs)
      Sys.sleep(5)
      ALLenriched = enrichr(names(gene_list), dbs)
      Sys.sleep(5)
      for (DB_name in dbs) {
        sig_terms_UP = UPenriched[[DB_name]][UPenriched[[DB_name]]$Adjusted.P.value<0.05,]
        filename <- paste( donor,  celltype, res_name_sanitized, DB_name, "UP_enrichr_significantONLY.csv", sep = "_")
        write.csv(sig_terms_UP, filename)
    
        Sys.sleep(1)
        sig_terms_down = DOWNenriched[[DB_name]][DOWNenriched[[DB_name]]$Adjusted.P.value<0.05,]
        filename <- paste( donor,  celltype, res_name_sanitized, DB_name, "DOWN_enrichr_significantONLY.csv", sep = "_")
        write.csv(sig_terms_down, filename)
    
        Sys.sleep(1)
        sig_terms_ALL = ALLenriched[[DB_name]][ALLenriched[[DB_name]]$Adjusted.P.value<0.05,]
        filename <- paste( donor,  celltype, res_name_sanitized, DB_name, "ALL_enrichr_significantONLY.csv", sep = "_")
        write.csv(sig_terms_ALL, filename)
        Sys.sleep(1)
      }
    
      m_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>%
        dplyr::select(gs_description, ensembl_gene )
      C7_GSEA <- GSEA(ensg_list, minGSSize = 5, pAdjustMethod = "BH", pvalueCutoff = 0.05, TERM2GENE = m_t2g)
      gse_df <- as.data.frame(C7_GSEA)
      gsea_csv_filename <- paste( donor, celltype, res_name_sanitized, "C7_GSEAResults.csv", sep = "_")
      write.csv(gse_df, file = gsea_csv_filename)
    
      m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>%
        dplyr::select(gs_description, ensembl_gene )
      C5_GSEA <- GSEA(ensg_list, minGSSize = 5, pAdjustMethod = "BH", pvalueCutoff = 0.05, TERM2GENE = m_t2g)
      gse_df <- as.data.frame(C5_GSEA)
      gsea_csv_filename <- paste( donor, celltype, res_name_sanitized, "C5_GSEAResults.csv", sep = "_")
      write.csv(gse_df, file = gsea_csv_filename)
    
      gse <- gseGO(geneList = ensg_list,
                   ont = "ALL",
                   keyType = "ENSEMBL",
                   nPerm = 10000,
                   minGSSize = 3,
                   maxGSSize = 800,
                   pvalueCutoff = 0.05,
                   verbose = TRUE,
                   OrgDb = organism,
                   pAdjustMethod = "BH")
      gse_df <- as.data.frame(gse)
      gsea_csv_filename <- paste(donor, celltype, res_name_sanitized, "GSEGO_GSEAResults.csv", sep = "_")
      write.csv(gse_df, file = gsea_csv_filename)
    
    }, silent = TRUE)
    
      }
    }
  }
}



# Compare Devices ------------------------------------------------------
dir.create('E:/AAA_Labwork/Tcell model RNAseq/donor_celltype_compareDEVICE')
setwd('E:/AAA_Labwork/Tcell model RNAseq/donor_celltype_compareDEVICE')
treatment = 'NO TREATMENT'
subset_metadata <- metadata %>%
  filter(Treatment == treatment)
unique_donor <- unique(metadata$Donor)
unique_celltypes <- unique(subset_metadata$Cell.Type)
unique_treatments <- unique(subset_metadata$Treatment)
unique_design = unique(subset_metadata$Study.Design)
for (donor in unique_donor) {
  for (celltype in unique_celltypes) {
    organism = "org.Hs.eg.db"
    subset_metadata <- metadata %>%
      filter(Donor == donor, Cell.Type == celltype, Treatment == treatment)
    
    if (length(unique(subset_metadata$Study.Design)) < 2) {
      print(paste('No second design for ', donor, celltype))
      next
    } else {
      print(paste('Analyze ', donor, celltype))
    }
    print(subset_metadata)
    subset_counts <- raw_counts[, subset_metadata$Sample.Name]
    
    dds <- DESeqDataSetFromMatrix(countData = subset_counts, colData = subset_metadata, design = ~ Study.Design)
    dds <- dds[rowSums(counts(dds)) >= 10, ]
    dds <- DESeq(dds)
    #dds <- DESeq(dds, minReplicatesForReplace=Inf)
    saveRDS(dds, paste0(donor, celltype, "_grouped.rds"))
    normalized_counts <- counts(dds, normalized = TRUE)
    genes_to_plot <- c('CYP2A7','CYP3A7','CYP2A6','ACSS2','SEC16B','GSTA2','BCHE','EPCAM','ASCL2','CLDN10','KRT20','MUC2','SLC26A2',
                       'CA1','HES1','SPIB','CHGA','CD27','NCAM1','CD86','CD1C','HLA.DPB1','IGHM',
                       'IGHA1','CD19','CD20','CD3G','CD3D','CD3E','TRDC','CD8A','KLRB1','TRAV1.2','KIT','PTPRC')
    split_names <- strsplit(row.names(normalized_counts), "_")
    string_array <- sapply(split_names, choose_id)
    rownames(normalized_counts) = string_array
    # Subset normalized counts for these genes
    normalized_subset <- normalized_counts[rownames(normalized_counts) %in% genes_to_plot, ]
    #breaks <- seq(0, 2000, length.out = 21)
    # Plot the heatmap
    dds_heat = pheatmap(normalized_subset, 
             #breaks = breaks,
             #color = colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1),
             scale = "none", # Normalize by row (gene-wise) to better visualize differences
             cluster_rows = TRUE, 
             cluster_cols = TRUE, 
             show_rownames = TRUE, 
             show_colnames = TRUE)
    ggsave(paste0(donor,'_',celltype,'_heatmap.png'), plot = dds_heat, width = 10, height = 10)
    
    res_names <- resultsNames(dds)
    
    for (res_name in res_names) {
      if (res_name != "Intercept") {
        res = results(dds,name = res_name)
        res_name_sanitized <- gsub("c.", "_", res_name)
        #     res.shrink <- lfcShrink(dds = dds, coef = res_name, type = "ashr")
        #     res_for_vp <- res.shrink
        res_for_vp = res
        #     split_names <- strsplit(row.names(res.shrink), "_")
        split_names <- strsplit(row.names(res_for_vp), "_")
        
        ENSGID <- sapply(split_names, choose_ENSG_id)
        string_array <- sapply(split_names, choose_id)
        easy_search_array <- sapply(split_names, easy_search)
        #     res.shrink$GeneSymbol <- easy_search_array
        #     res.shrink$ENSGID <- ENSGID
        res$GeneSymbol <- easy_search_array
        res$ENSGID <- ENSGID
        
        row.names(res_for_vp) <- string_array
        volcano_labels <- ifelse(res_for_vp$padj <= 0.05 & abs(res_for_vp$log2FoldChange) > 1, rownames(res_for_vp), '')
        volcano_labels[is.na(volcano_labels)] <- ""
        volcano_labels <- ifelse(!startsWith(string_array, 'ENSG'), volcano_labels, '')
        
        volcanoplot <- EnhancedVolcano(
          res_for_vp,
          lab = '',  # No labels here to avoid duplication
          x = 'log2FoldChange',
          y = 'padj',
          ylab = bquote(~-Log[10] ~ italic(Padj)),
          pCutoff = 0.05,
          FCcutoff = 1,
          pointSize = 1.0,
          labSize = 0.5,
          drawConnectors = TRUE,
          widthConnectors = 0.5,
          max.overlaps = 30
        )
        
        volcanoplot <- volcanoplot + geom_text_repel(
          aes(label = volcano_labels),
          size = 2,  # Adjust size as needed
          segment.color = 'grey50',
          color = 'black',
          bg.color = 'cyan',  # Background color for the outline
          bg.r = 0.1,  # Radius for the outline
          max.overlaps = 70,  # Allow more overlaps
          min.segment.length = 0.5  # Minimum segment length for label connectors
        )
        
        volcano_filename <- paste0( donor, celltype, "_", res_name_sanitized, "DesignCompare_volcanoplot.png")
        ggsave(volcano_filename, plot = volcanoplot, width = 5, height = 7)
        
        # resOrdered <- res.shrink[order(res.shrink$padj), ]
        resOrdered <- res[order(res$padj), ]
        
        filename <- paste( donor, celltype, res_name_sanitized, "DesignCompare.csv", sep = "_")
        write.csv(as.data.frame(resOrdered), file = filename)
        
        significant_results <- subset(resOrdered, padj <= 0.05)
        significant_results<- significant_results[order(significant_results$log2FoldChange), ]
        
        filename <- paste( donor, celltype, res_name_sanitized, "DesignCompare_significantONLY_sortedFC.csv", sep = "_")
        write.csv(as.data.frame(significant_results), file = filename)
        
        gene_list <- significant_results$log2FoldChange
        names(gene_list) <- rownames(significant_results)
        
        gene_list <- sort(gene_list, decreasing = TRUE)
        ensg_list = gene_list
        
        names(gene_list) <- sapply(strsplit(names(gene_list), "_"), choose_id)
        gene_list <- gene_list[!startsWith(names(gene_list), 'ENSG')]
        
        names(ensg_list) <- sapply(strsplit(names(ensg_list), "_"),choose_ENSG_id)
        print(ensg_list)
        
        if (length(gene_list) == 0) {
          message("No significant genes found for ", donor,  celltype, " ", res_name_sanitized)
          next
        }
        
        try({
          
          uplist = names(gene_list[gene_list>0])
          downlist = names(gene_list[gene_list<0])
          UPenriched <- enrichr(uplist, dbs)
          Sys.sleep(5)
          DOWNenriched = enrichr(downlist, dbs)
          Sys.sleep(5)
          ALLenriched = enrichr(names(gene_list), dbs)
          Sys.sleep(5)
          for (DB_name in dbs) {
            sig_terms_UP = UPenriched[[DB_name]][UPenriched[[DB_name]]$Adjusted.P.value<0.05,]
            filename <- paste( donor,  celltype, res_name_sanitized, DB_name, "UP_enrichr_significantONLY.csv", sep = "_")
            write.csv(sig_terms_UP, filename)
            
            Sys.sleep(1)
            sig_terms_down = DOWNenriched[[DB_name]][DOWNenriched[[DB_name]]$Adjusted.P.value<0.05,]
            filename <- paste( donor,  celltype, res_name_sanitized, DB_name, "DOWN_enrichr_significantONLY.csv", sep = "_")
            write.csv(sig_terms_down, filename)
        
            Sys.sleep(1)
            sig_terms_ALL = ALLenriched[[DB_name]][ALLenriched[[DB_name]]$Adjusted.P.value<0.05,]
            filename <- paste( donor,  celltype, res_name_sanitized, DB_name, "ALL_enrichr_significantONLY.csv", sep = "_")
            write.csv(sig_terms_ALL, filename)
            Sys.sleep(1)
          }
          
          m_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>%
            dplyr::select(gs_description, ensembl_gene )
          C7_GSEA <- GSEA(ensg_list, minGSSize = 5, pAdjustMethod = "BH", pvalueCutoff = 0.05, TERM2GENE = m_t2g)
          gse_df <- as.data.frame(C7_GSEA)
          gsea_csv_filename <- paste( donor, celltype, res_name_sanitized, "C7_GSEAResults.csv", sep = "_")
          write.csv(gse_df, file = gsea_csv_filename)

          m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>%
            dplyr::select(gs_description, ensembl_gene )
          C5_GSEA <- GSEA(ensg_list, minGSSize = 5, pAdjustMethod = "BH", pvalueCutoff = 0.05, TERM2GENE = m_t2g)
          gse_df <- as.data.frame(C5_GSEA)
          gsea_csv_filename <- paste( donor, celltype, res_name_sanitized, "C5_GSEAResults.csv", sep = "_")
          write.csv(gse_df, file = gsea_csv_filename)

          gse <- gseGO(geneList = ensg_list,
                       ont = "ALL",
                       keyType = "ENSEMBL",
                       nPerm = 10000,
                       minGSSize = 3,
                       maxGSSize = 800,
                       pvalueCutoff = 0.05,
                       verbose = TRUE,
                       OrgDb = organism,
                       pAdjustMethod = "BH")
          gse_df <- as.data.frame(gse)
          gsea_csv_filename <- paste(donor, celltype, res_name_sanitized, "GSEGO_GSEAResults.csv", sep = "_")
          write.csv(gse_df, file = gsea_csv_filename)

        }, silent = TRUE)
        
      }
    }
  }
}

# Compare Devices For GUT EPI, GUT EPI and GUT CD45, GUT EPI (with liver and no CD45+), and GUT EPI and GUT CD45 with full interacting chambers ------------------------------------------------------
dir.create('E:/AAA_Labwork/Tcell model RNAseq/donor_celltype_compareDEVICE_EPIi')
setwd('E:/AAA_Labwork/Tcell model RNAseq/donor_celltype_compareDEVICE_EPIi')
treatment = 'NO TREATMENT'
subset_metadata <- metadata %>%
  filter(Treatment == treatment)
unique_donor <- unique(metadata$Donor)
unique_celltypes = c("GUT EPI and GUT EPI CD45+", "GUT EPI","GUT EPI CD45+")
unique_design = unique(subset_metadata$Study.Design)
for (donor in unique_donor) {
  subset_metadata <- metadata %>%
    filter(Donor == donor, Cell.Type %in% unique_celltypes, Treatment == treatment)
  print(subset_metadata)
  subset_metadata$Study.Design <- factor(subset_metadata$Study.Design)
  subset_metadata$Study.Design <- relevel(subset_metadata$Study.Design, ref = "ISOLATION (ISO) PURE")#"INTERACTION (INT) GUT-LIVER-IMMUNE CELLS")
  #subset_metadata$Study.Design <- relevel(subset_metadata$Study.Design, ref = "INTERACTION (INT) GUT-LIVER-IMMUNE CELLS")
  
  subset_counts <- raw_counts[, subset_metadata$Sample.Name]
  
  dds <- DESeqDataSetFromMatrix(countData = subset_counts, colData = subset_metadata, design = ~ Study.Design)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds)
  #dds <- DESeq(dds, minReplicatesForReplace=Inf)
  saveRDS(dds, paste0(donor, "_EPI_grouped.rds"))
  normalized_counts <- counts(dds, normalized = TRUE)
  genes_to_plot <- c('CYP2A7','CYP3A7','CYP2A6','ACSS2','SEC16B','GSTA2','BCHE','EPCAM','ASCL2','CLDN10','KRT20','MUC2','SLC26A2',
                     'CA1','HES1','SPIB','CHGA','CD27','NCAM1','CD86','CD1C','HLA.DPB1','IGHM',
                     'IGHA1','CD19','CD20','CD3G','CD3D','CD3E','TRDC','CD8A','KLRB1','TRAV1.2','KIT','PTPRC')
  split_names <- strsplit(row.names(normalized_counts), "_")
  string_array <- sapply(split_names, choose_id)
  rownames(normalized_counts) = string_array
  # Subset normalized counts for these genes
  normalized_subset <- normalized_counts[rownames(normalized_counts) %in% genes_to_plot, ]
  #breaks <- seq(0, 2000, length.out = 21)
  # Plot the heatmap
  dds_heat = pheatmap(normalized_subset, 
                      #breaks = breaks,
                      #color = colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1),
                      scale = "none", # Normalize by row (gene-wise) to better visualize differences
                      cluster_rows = TRUE, 
                      cluster_cols = TRUE, 
                      show_rownames = TRUE, 
                      show_colnames = TRUE)
  ggsave(paste0(donor,'_EPIi_heatmap.png'), plot = dds_heat, width = 10, height = 10)
  
  res_names <- resultsNames(dds)
  
  for (res_name in res_names) {
    if (res_name != "Intercept") {
      
      res = results(dds,name = res_name)
      res_name_sanitized <- gsub("c.", "_", res_name)
      #     res.shrink <- lfcShrink(dds = dds, coef = res_name, type = "ashr")
      #     res_for_vp <- res.shrink
      res_for_vp = res
      #     split_names <- strsplit(row.names(res.shrink), "_")
      split_names <- strsplit(row.names(res_for_vp), "_")
      
      ENSGID <- sapply(split_names, choose_ENSG_id)
      string_array <- sapply(split_names, choose_id)
      easy_search_array <- sapply(split_names, easy_search)
      #     res.shrink$GeneSymbol <- easy_search_array
      #     res.shrink$ENSGID <- ENSGID
      res$GeneSymbol <- easy_search_array
      res$ENSGID <- ENSGID
      
      row.names(res_for_vp) <- string_array
      volcano_labels <- ifelse(res_for_vp$padj <= 0.05 & abs(res_for_vp$log2FoldChange) > 1, rownames(res_for_vp), '')
      volcano_labels[is.na(volcano_labels)] <- ""
      volcano_labels <- ifelse(!startsWith(string_array, 'ENSG'), volcano_labels, '')
      
      volcanoplot <- EnhancedVolcano(
        res_for_vp,
        lab = '',  # No labels here to avoid duplication
        x = 'log2FoldChange',
        y = 'padj',
        ylab = bquote(~-Log[10] ~ italic(Padj)),
        pCutoff = 0.05,
        FCcutoff = 1,
        pointSize = 1.0,
        labSize = 0.5,
        drawConnectors = TRUE,
        widthConnectors = 0.5,
        max.overlaps = 30
      )
      
      volcanoplot <- volcanoplot + geom_text_repel(
        aes(label = volcano_labels),
        size = 2,  # Adjust size as needed
        segment.color = 'grey50',
        color = 'black',
        bg.color = 'cyan',  # Background color for the outline
        bg.r = 0.1,  # Radius for the outline
        max.overlaps = 70,  # Allow more overlaps
        min.segment.length = 0.5  # Minimum segment length for label connectors
      )
      
      volcano_filename <- paste0( donor,  "_", res_name_sanitized, "EPIi_DesignCompare_volcanoplot.png")
      ggsave(volcano_filename, plot = volcanoplot, width = 5, height = 7)
      
      # resOrdered <- res.shrink[order(res.shrink$padj), ]
      resOrdered <- res[order(res$padj), ]
      
      filename <- paste( donor, res_name_sanitized, "EPIi_DesignCompare.csv", sep = "_")
      write.csv(as.data.frame(resOrdered), file = filename)
      
      significant_results <- subset(resOrdered, padj <= 0.05)
      significant_results<- significant_results[order(significant_results$log2FoldChange), ]
      
      filename <- paste( donor, res_name_sanitized, "EPIi_DesignCompare_significantONLY_sortedFC.csv", sep = "_")
      write.csv(as.data.frame(significant_results), file = filename)
      
      gene_list <- significant_results$log2FoldChange
      names(gene_list) <- rownames(significant_results)
      
      gene_list <- sort(gene_list, decreasing = TRUE)
      ensg_list = gene_list
      
      names(gene_list) <- sapply(strsplit(names(gene_list), "_"), choose_id)
      gene_list <- gene_list[!startsWith(names(gene_list), 'ENSG')]
      
      names(ensg_list) <- sapply(strsplit(names(ensg_list), "_"),choose_ENSG_id)
      print(ensg_list)
      
      if (length(gene_list) == 0) {
        message("No significant genes found for ", donor, " EPIi_", res_name_sanitized)
        next
      }
      
      try({
        
        uplist = names(gene_list[gene_list>0])
        downlist = names(gene_list[gene_list<0])
        UPenriched <- enrichr(uplist, dbs)
        Sys.sleep(5)
        DOWNenriched = enrichr(downlist, dbs)
        Sys.sleep(5)
        ALLenriched = enrichr(names(gene_list), dbs)
        Sys.sleep(5)
        for (DB_name in dbs) {
          sig_terms_UP = UPenriched[[DB_name]][UPenriched[[DB_name]]$Adjusted.P.value<0.05,]
          filename <- paste( donor, res_name_sanitized, DB_name, "EPIi_UP_enrichr_significantONLY.csv", sep = "_")
          write.csv(sig_terms_UP, filename)
          
          Sys.sleep(1)
          sig_terms_down = DOWNenriched[[DB_name]][DOWNenriched[[DB_name]]$Adjusted.P.value<0.05,]
          filename <- paste( donor, res_name_sanitized, DB_name, "EPIi_DOWN_enrichr_significantONLY.csv", sep = "_")
          write.csv(sig_terms_down, filename)
          
          Sys.sleep(1)
          sig_terms_ALL = ALLenriched[[DB_name]][ALLenriched[[DB_name]]$Adjusted.P.value<0.05,]
          filename <- paste( donor, res_name_sanitized, DB_name, "EPIi_ALL_enrichr_significantONLY.csv", sep = "_")
          write.csv(sig_terms_ALL, filename)
          Sys.sleep(1)
        }
        
        m_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>%
          dplyr::select(gs_description, ensembl_gene )
        C7_GSEA <- GSEA(ensg_list, minGSSize = 5, pAdjustMethod = "BH", pvalueCutoff = 0.05, TERM2GENE = m_t2g)
        gse_df <- as.data.frame(C7_GSEA)
        gsea_csv_filename <- paste( donor, celltype, res_name_sanitized, "C7_GSEAResults.csv", sep = "_")
        write.csv(gse_df, file = gsea_csv_filename)

        m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>%
          dplyr::select(gs_description, ensembl_gene )
        C5_GSEA <- GSEA(ensg_list, minGSSize = 5, pAdjustMethod = "BH", pvalueCutoff = 0.05, TERM2GENE = m_t2g)
        gse_df <- as.data.frame(C5_GSEA)
        gsea_csv_filename <- paste( donor, celltype, res_name_sanitized, "C5_GSEAResults.csv", sep = "_")
        write.csv(gse_df, file = gsea_csv_filename)
        
        gse <- gseGO(geneList = ensg_list,
                     ont = "ALL",
                     keyType = "ENSEMBL",
                     nPerm = 10000,
                     minGSSize = 3,
                     maxGSSize = 800,
                     pvalueCutoff = 0.05,
                     verbose = TRUE,
                     OrgDb = organism,
                     pAdjustMethod = "BH")
        gse_df <- as.data.frame(gse)
        gsea_csv_filename <- paste(donor, celltype, res_name_sanitized, "GSEGO_GSEAResults.csv", sep = "_")
        write.csv(gse_df, file = gsea_csv_filename)

      }, silent = TRUE)
      
    }
  }
}



# Plot heatmap ------------------------------------------------------------
HEP_dds = readRDS('D3LIVER CD45+_grouped.rds')
normalized_counts <- counts(HEP_dds, normalized = TRUE)
# List of genes of interest
genes_to_plot <- c('CYP2A7','CYP3A7','CYP2A6','ACSS2','SEC16B','GSTA2','BCHE','EPCAM','ASCL2','CLDN10','KRT20','MUC2','SLC26A2',
                   'CA1','HES1','SPIB','CHGA','CD27','NCAM1','CD86','CD1C','HLA.DPB1','IGHM',
                   'IGHA1','CD19','CD20','CD3G','CD3D','CD3E','TRDC','CD8A','KLRB1','TRAV1.2','KIT','PTPRC')
split_names <- strsplit(row.names(normalized_counts), "_")
string_array <- sapply(split_names, choose_id)
rownames(normalized_counts) = string_array
# Subset normalized counts for these genes
normalized_subset <- normalized_counts[rownames(normalized_counts) %in% genes_to_plot, ]
library(pheatmap)
#breaks <- seq(0, 2000, length.out = 21)
# Plot the heatmap
pheatmap(normalized_subset, 
         #breaks = breaks,
         #color = colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1),
         scale = "none", # Normalize by row (gene-wise) to better visualize differences
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE)

dds = readRDS('D3GUT LP CD45+_grouped.rds')

# Plot dotplot ------------------------------------------------------------
#setwd('DEG/donor_celltype_compareDESIGN')
setwd('DEG/donor_celltype_compareTREAT')
enrichfiles = list.files(path = 'DEG/donor_celltype_compareDESIGN', 
                         pattern = "*enrichr.*\\.csv$")

enrichfiles = list.files(path = 'DEG/donor_celltype_compareTREAT', 
                         pattern = "*enrichr.*\\.csv$")

# Function to safely evaluate the "Overlap" column
evaldiv <- function(x) {
  eval(parse(text = x))
}

for (enrichfile in enrichfiles) {
  
  try({
  enrichdf = read.csv(enrichfile, row.names = 1)
  
  # Skip the file if the dataframe is empty
  if (nrow(enrichdf) == 0) {
    next
  }
  
  # Limit to first 8 terms if there are more than 8
  if (nrow(enrichdf) > 8) {
    enrichdf = enrichdf[1:8, ]
  }
  
  # Convert Overlap to numeric
  enrichdf$Overlap_num = sapply(enrichdf$Overlap, evaldiv)
  
  # Create the dot plot
  dotplot_enrichr = ggplot(enrichdf, aes(x = -log10(Adjusted.P.value), y = Term)) +
    geom_point(aes(size = Overlap_num, color = Odds.Ratio)) +
    scale_color_gradient(low = "blue", high = "red") +  # Customize color scale
    theme_minimal() +
    labs(
      title = "Dot Plot for GSEA Results",
      x = "-log10(Adjusted P-Value)",
      y = "Term",
      size = "Set Overlap",
      color = "Odds Ratio"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Generate output filename
  output_filename <- gsub(".csv", "_dotplot.png", enrichfile)
  
  # Save the plot
  ggsave(output_filename, plot = dotplot_enrichr, width = 14, height = 6)
}, silent = TRUE)
}

# Plot dotplot C5 C7 and GSEGO ------------------------------------------------------------
setwd('DEG/donor_celltype_compareDESIGN')
enrichfiles = list.files(path = 'DEG/donor_celltype_compareDESIGN', 
                         pattern = "*GSEAResults*\\.csv$",
                         recursive =TRUE)

# Function to safely evaluate the "Overlap" column
evaldiv <- function(x) {
  eval(parse(text = x))
}

for (enrichfile in enrichfiles) {
  
  try({
    enrichdf = read.csv(enrichfile)
    
    # Skip the file if the dataframe is empty
    if (nrow(enrichdf) == 0) {
      next
    }
    
    # Limit to first 8 terms if there are more than 8
    if (nrow(enrichdf) > 8) {
      enrichdf = enrichdf[1:8, ]
    }
    
    # Create the dot plot
    dotplot_enrichr = ggplot(enrichdf, aes(x = -log10(p.adjust), y = ID)) +
      geom_point(aes(color = enrichmentScore)) +
      scale_color_gradient(low = "blue", high = "red") +  # Customize color scale
      theme_minimal() +
      labs(
        title = "Dot Plot for GSEA Results",
        x = "-log10(Adjusted P-Value)",
        y = "Term",
        size = "Set Overlap",
        color = "Odds Ratio"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Generate output filename
    output_filename <- gsub(".csv", "_dotplot.png", enrichfile)
    
    # Save the plot
    ggsave(output_filename, plot = dotplot_enrichr, width = 14, height = 6)
  }, silent = TRUE)
}

# Select Terms ------------------------------------------------------------
# Define helper function to evaluate the "Overlap" column
# Define helper function to evaluate the "Overlap" column

enrichfiles = list.files(pattern = "Reactome.*\\.csv$", recursive = TRUE)
interested_terms = c("T cell","lymphocyte","leukocyte","antigen","chemo","CD","receptor","cytokine","effector",
                     "nterleukin","MHC","IL ","Th1","Th2","stress","Stress","signaling","Signaling")
unwanted_terms = c("SARS","B cell","immunoglobulin","Cancer","cancer","IL 24","By Interleukins","Immune System","Treg","Family")

evaldiv <- function(x) {
  eval(parse(text = x))
}

for (enrichfile in enrichfiles) {
  
  try({
    enrichdf = read.csv(enrichfile)
    # Generate output filename
    output_filename = gsub("_enrichr_significantONLY", "", enrichfile)
    output_filename = gsub("_2023_Human", "", output_filename)
    output_filename <- gsub(".csv", "_barplot.jpg", output_filename)
    
    pattern <- paste(interested_terms, collapse = "|")
    unwanted_pattern <- paste(unwanted_terms, collapse = "|")
    
    # Filter the dataframe to keep rows where the 'Description' column contains any of the interested terms
    filtered_gsea_df <- enrichdf[grepl(pattern, enrichdf$Term, ignore.case = TRUE), ]
    filtered_gsea_df <- filtered_gsea_df[!grepl(unwanted_pattern, filtered_gsea_df$Term, ignore.case = TRUE), ]
    
    if (nrow(filtered_gsea_df) < 2) {
      filtered_gsea_df = enrichdf
    }
    
    # Clean up the 'Term' column
    filtered_gsea_df$Term <- gsub("R-HSA-\\d+", "", filtered_gsea_df$Term)
    filtered_gsea_df$Term <- trimws(filtered_gsea_df$Term)
    filtered_gsea_df$Overlap_num = sapply(filtered_gsea_df$Overlap, evaldiv)
    filtered_gsea_df$Term <- gsub("(.{30})\\s", "\\1\n", filtered_gsea_df$Term)
    
    if (nrow(filtered_gsea_df) > 10) {
      filtered_gsea_df = filtered_gsea_df[1:10, ]
    }
    
    # Set color and axis position based on whether the output file contains "UP"
    if (grepl("UP", output_filename)) {
      filtered_gsea_df$color <- 'red'
      x_position <- "bottom"
      y_position <- "left"
    } else {
      filtered_gsea_df$color <- 'blue'
      filtered_gsea_df$Combined.Score <- filtered_gsea_df$Combined.Score * -1
      
      x_position <- "top"
      y_position <- "right"
    }
    
    filtered_gsea_df$log.padj = round(-log10(filtered_gsea_df$Adjusted.P.value),3)
    
    barplot_enrichr = ggplot(filtered_gsea_df, aes(x = Combined.Score, y = reorder(Term, Combined.Score))) +
      geom_bar(aes(fill = color), stat = "identity") +  # Use 'fill' based on the color column
      # geom_text(aes(label = log.padj),               # Add overlap number labels # Position label in the middle of the bar
      #           hjust = -0.3,                           # Adjust horizontal alignment of labels
      #           size = 3, fontface  = "bold") +                             # Set label font size
      scale_fill_identity() +  # Use the exact colors from the 'color' column without mapping them
      theme_minimal() +
      labs(
        title = paste0("Bar Plot for GSEA Results - ", gsub("_", " ", gsub("_dotplot.jpg", "", output_filename))),
        x = "Enrichment",
        y = "Terms"
      ) +
      scale_x_continuous(position = x_position) +  # Move x-axis to the top if necessary
      scale_y_discrete(position = y_position) +  # Move y-axis to the right if necessary
      theme(
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),  # Adjust x-axis text font size and rotation
        axis.text.y = element_text(size = 12, face = "bold"),         # Adjust y-axis text font size
        plot.title = element_text(size = 16, face = "bold"),          # Adjust plot title font size and bold
        legend.text = element_text(size = 10),                        # Adjust legend text font size
        legend.title = element_text(size = 12)                        # Adjust legend title font size
      )
    
    # Save the plot
    ggsave(output_filename, plot = barplot_enrichr, width = 9, height = 5*nrow(filtered_gsea_df)/12+1)
    
  }, silent = TRUE)
}


# Combined donor test -----------------------------------------------------
setwd('DEG/celltype_compareTREAT')
studydesign <- 'INTERACTION (INT) GUT-LIVER-IMMUNE CELLS'
subset_metadata <- metadata %>%
  filter(Study.Design == studydesign)
unique_celltypes <- unique(subset_metadata$Cell.Type)
unique_treatments <- unique(subset_metadata$Treatment)
unique_design = unique(subset_metadata$Study.Design)
for (celltype in unique_celltypes) {
  organism = "org.Hs.eg.db"
  subset_metadata <- metadata %>%
    filter(Cell.Type == celltype, Study.Design == studydesign)
  # we will go check here: this is the step where the samples used for in-group comparison is selected.
  subset_counts <- raw_counts[, subset_metadata$Sample.Name]
  print(subset_metadata)
  dds <- DESeqDataSetFromMatrix(countData = subset_counts, colData = subset_metadata, design = ~ Treatment)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds)
  #dds <- DESeq(dds, minReplicatesForReplace=Inf)
  saveRDS(dds, paste0( celltype, "_grouped.rds"))
  
  normalized_counts <- counts(dds, normalized = TRUE)
  genes_to_plot <- c('CYP2A7','CYP3A7','CYP2A6','ACSS2','SEC16B','GSTA2','BCHE','EPCAM','ASCL2','CLDN10','KRT20','MUC2','SLC26A2',
                     'CA1','HES1','SPIB','CHGA','CD27','NCAM1','CD86','CD1C','HLA.DPB1','IGHM',
                     'IGHA1','CD19','CD20','CD3G','CD3D','CD3E','TRDC','CD8A','KLRB1','TRAV1.2','KIT','PTPRC')
  split_names <- strsplit(row.names(normalized_counts), "_")
  string_array <- sapply(split_names, choose_id)
  rownames(normalized_counts) = string_array
  # Subset normalized counts for these genes
  normalized_subset <- normalized_counts[rownames(normalized_counts) %in% genes_to_plot, ]
  #breaks <- seq(0, 2000, length.out = 21)
  # Plot the heatmap
  dds_heat = pheatmap(normalized_subset,
                      #breaks = breaks,
                      #color = colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1),
                      scale = "none", # Normalize by row (gene-wise) to better visualize differences
                      cluster_rows = TRUE,
                      cluster_cols = TRUE,
                      show_rownames = TRUE,
                      show_colnames = TRUE)
  ggsave(paste0(celltype,'_heatmap.png'), plot = dds_heat, width = 10, height = 10)
  
  res_names <- resultsNames(dds)
  
  for (res_name in res_names) {
    if (res_name != "Intercept") {
      
      res_name_sanitized <- gsub("c.", "_", res_name)
      res.shrink <- lfcShrink(dds = dds, coef = res_name, type = "ashr")
      res_for_vp <- res.shrink
      split_names <- strsplit(row.names(res.shrink), "_")
      
      ENSGID <- sapply(split_names, choose_ENSG_id)
      string_array <- sapply(split_names, choose_id)
      easy_search_array <- sapply(split_names, easy_search)
      res.shrink$GeneSymbol <- easy_search_array
      res.shrink$ENSGID <- ENSGID
      
      row.names(res_for_vp) <- string_array
      volcano_labels <- ifelse(res_for_vp$padj <= 0.05 & abs(res_for_vp$log2FoldChange) > 1, rownames(res_for_vp), '')
      volcano_labels[is.na(volcano_labels)] <- ""
      volcano_labels <- ifelse(!startsWith(string_array, 'ENSG'), volcano_labels, '')
      
      volcanoplot <- EnhancedVolcano(
        res_for_vp,
        lab = '',  # No labels here to avoid duplication
        x = 'log2FoldChange',
        y = 'padj',
        ylab = bquote(~-Log[10] ~ italic(Padj)),
        pCutoff = 0.05,
        FCcutoff = 1,
        pointSize = 1.0,
        labSize = 0.5,
        drawConnectors = TRUE,
        widthConnectors = 0.5,
        max.overlaps = 30
      )
      
      volcanoplot <- volcanoplot + geom_text_repel(
        aes(label = volcano_labels),
        size = 2,  # Adjust size as needed
        segment.color = 'grey50',
        color = 'black',
        bg.color = 'cyan',  # Background color for the outline
        bg.r = 0.1,  # Radius for the outline
        max.overlaps = 70,  # Allow more overlaps
        min.segment.length = 0.5  # Minimum segment length for label connectors
      )
      
      volcano_filename <- paste0(celltype, "_", res_name_sanitized, "volcanoplot.png")
      ggsave(volcano_filename, plot = volcanoplot, width = 5, height = 7)
      
      resOrdered <- res.shrink[order(res.shrink$padj), ]
      
      filename <- paste(celltype, res_name_sanitized, "TreatmentCompare.csv", sep = "_")
      write.csv(as.data.frame(resOrdered), file = filename)
      
      significant_results <- subset(resOrdered, padj <= 0.05)
      significant_results<- significant_results[order(significant_results$log2FoldChange), ]
      
      filename <- paste(celltype, "TreatmentCompare_significantONLY_sortedFC.csv", sep = "_")
      write.csv(as.data.frame(significant_results), file = filename)
      
      gene_list <- significant_results$log2FoldChange
      names(gene_list) <- rownames(significant_results)
      
      gene_list <- sort(gene_list, decreasing = TRUE)
      ensg_list = gene_list
      
      names(gene_list) <- sapply(strsplit(names(gene_list), "_"), choose_id)
      gene_list <- gene_list[!startsWith(names(gene_list), 'ENSG')]
      
      names(ensg_list) <- sapply(strsplit(names(ensg_list), "_"),choose_ENSG_id)
      print(ensg_list)
      
      if (length(gene_list) == 0) {
        message("No significant genes found for ", celltype, " ", res_name_sanitized)
        next
      }
      
      try({
        
        # uplist = names(gene_list[gene_list>0])
        # downlist = names(gene_list[gene_list<0])
        # UPenriched <- enrichr(uplist, dbs)
        # Sys.sleep(5)
        # DOWNenriched = enrichr(downlist, dbs)
        # Sys.sleep(5)
        # ALLenriched = enrichr(names(gene_list), dbs)
        # Sys.sleep(5)
        # for (DB_name in dbs) {
        #   sig_terms_UP = UPenriched[[DB_name]][UPenriched[[DB_name]]$Adjusted.P.value<0.05,]
        #   filename <- paste(celltype, res_name_sanitized, DB_name, "UP_enrichr_significantONLY.csv", sep = "_")
        #   write.csv(sig_terms_UP, filename)
        # 
        #   Sys.sleep(1)
        #   sig_terms_down = DOWNenriched[[DB_name]][DOWNenriched[[DB_name]]$Adjusted.P.value<0.05,]
        #   filename <- paste(celltype, res_name_sanitized, DB_name, "DOWN_enrichr_significantONLY.csv", sep = "_")
        #   write.csv(sig_terms_down, filename)
        # 
        #   Sys.sleep(1)
        #   sig_terms_ALL = ALLenriched[[DB_name]][ALLenriched[[DB_name]]$Adjusted.P.value<0.05,]
        #   filename <- paste( celltype, res_name_sanitized, DB_name, "ALL_enrichr_significantONLY.csv", sep = "_")
        #   write.csv(sig_terms_ALL, filename)
        #   Sys.sleep(1)
        # }
        
        # m_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>%
        #   dplyr::select(gs_description, ensembl_gene )
        # C7_GSEA <- GSEA(ensg_list, minGSSize = 5, pAdjustMethod = "BH", pvalueCutoff = 0.05, TERM2GENE = m_t2g)
        # gse_df <- as.data.frame(C7_GSEA)
        # gsea_csv_filename <- paste( donor, celltype, res_name_sanitized, "C7_GSEAResults.csv", sep = "_")
        # write.csv(gse_df, file = gsea_csv_filename)
        #
        # m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>%
        #   dplyr::select(gs_description, ensembl_gene )
        # C5_GSEA <- GSEA(ensg_list, minGSSize = 5, pAdjustMethod = "BH", pvalueCutoff = 0.05, TERM2GENE = m_t2g)
        # gse_df <- as.data.frame(C5_GSEA)
        # gsea_csv_filename <- paste( donor, celltype, res_name_sanitized, "C5_GSEAResults.csv", sep = "_")
        # write.csv(gse_df, file = gsea_csv_filename)
        #
        gse <- gseGO(geneList = ensg_list,
                     ont = "ALL",
                     keyType = "ENSEMBL",
                     nPerm = 10000,
                     minGSSize = 3,
                     maxGSSize = 800,
                     pvalueCutoff = 0.05,
                     verbose = TRUE,
                     OrgDb = organism,
                     pAdjustMethod = "BH")
        gse_df <- as.data.frame(gse)
        gsea_csv_filename <- paste(donor, celltype, res_name_sanitized, "GSEGO_GSEAResults.csv", sep = "_")
        write.csv(gse_df, file = gsea_csv_filename)
        
      }, silent = TRUE)
      
    }
  }
}

# Combined donor -- Devices For GUT EPI, GUT EPI and GUT CD45, GUT EPI (with liver and without CD45+), and GUT EPI and GUT CD45 with full interacting chambers ------------------------------------------------------

treatment = 'NO TREATMENT'
subset_metadata <- metadata %>%
  filter(Treatment == treatment)
unique_donor <- unique(metadata$Donor)
unique_celltypes = c("GUT EPI and GUT EPI CD45+", "GUT EPI","GUT EPI CD45+")
unique_design = unique(subset_metadata$Study.Design)
setwd('DEG/compareDESIGN_EPIi')

subset_metadata <- metadata %>%
  filter(Cell.Type %in% unique_celltypes, Treatment == treatment)
print(subset_metadata)
subset_metadata$Study.Design <- factor(subset_metadata$Study.Design)
subset_metadata$Study.Design <- relevel(subset_metadata$Study.Design, ref = "ISOLATION (ISO) PURE")#"INTERACTION (INT) GUT-LIVER-IMMUNE CELLS")

subset_counts <- raw_counts[, subset_metadata$Sample.Name]

dds <- DESeqDataSetFromMatrix(countData = subset_counts, colData = subset_metadata, design = ~ Study.Design)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
#dds <- DESeq(dds, minReplicatesForReplace=Inf)
saveRDS(dds, "EPI_grouped.rds")
normalized_counts <- counts(dds, normalized = TRUE)
genes_to_plot <- c('CYP2A7','CYP3A7','CYP2A6','ACSS2','SEC16B','GSTA2','BCHE','EPCAM','ASCL2','CLDN10','KRT20','MUC2','SLC26A2',
                   'CA1','HES1','SPIB','CHGA','CD27','NCAM1','CD86','CD1C','HLA.DPB1','IGHM',
                   'IGHA1','CD19','CD20','CD3G','CD3D','CD3E','TRDC','CD8A','KLRB1','TRAV1.2','KIT','PTPRC')
split_names <- strsplit(row.names(normalized_counts), "_")
string_array <- sapply(split_names, choose_id)
rownames(normalized_counts) = string_array
# Subset normalized counts for these genes
normalized_subset <- normalized_counts[rownames(normalized_counts) %in% genes_to_plot, ]
#breaks <- seq(0, 2000, length.out = 21)
# Plot the heatmap
dds_heat = pheatmap(normalized_subset, 
                    #breaks = breaks,
                    #color = colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1),
                    scale = "none", # Normalize by row (gene-wise) to better visualize differences
                    cluster_rows = TRUE, 
                    cluster_cols = TRUE, 
                    show_rownames = TRUE, 
                    show_colnames = TRUE)
ggsave(paste0('EPIi_heatmap.png'), plot = dds_heat, width = 10, height = 10)

res_names <- resultsNames(dds)

for (res_name in res_names) {
  if (res_name != "Intercept") {
    
    res_name_sanitized <- gsub("\\.", "_", res_name)
    res.shrink <- lfcShrink(dds = dds, coef = res_name, type = "ashr")
    res_for_vp <- res.shrink
    split_names <- strsplit(row.names(res.shrink), "_")
    
    ENSGID <- sapply(split_names, choose_ENSG_id)
    string_array <- sapply(split_names, choose_id)
    easy_search_array <- sapply(split_names, easy_search)
    res.shrink$GeneSymbol <- easy_search_array
    res.shrink$ENSGID <- ENSGID
    
    row.names(res_for_vp) <- string_array
    volcano_labels <- ifelse(res_for_vp$padj <= 0.05 & abs(res_for_vp$log2FoldChange) > 1, rownames(res_for_vp), '')
    volcano_labels[is.na(volcano_labels)] <- ""
    volcano_labels <- ifelse(!startsWith(string_array, 'ENSG'), volcano_labels, '')
    
    volcanoplot <- EnhancedVolcano(
      res_for_vp,
      lab = '',  # No labels here to avoid duplication
      x = 'log2FoldChange',
      y = 'padj',
      ylab = bquote(~-Log[10] ~ italic(Padj)),
      pCutoff = 0.05,
      FCcutoff = 1,
      pointSize = 1.0,
      labSize = 0.5,
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      max.overlaps = 30
    )
    
    volcanoplot <- volcanoplot + geom_text_repel(
      aes(label = volcano_labels),
      size = 2,  # Adjust size as needed
      segment.color = 'grey50',
      color = 'black',
      bg.color = 'cyan',  # Background color for the outline
      bg.r = 0.1,  # Radius for the outline
      max.overlaps = 70,  # Allow more overlaps
      min.segment.length = 0.5  # Minimum segment length for label connectors
    )
    
    volcano_filename <- paste0( celltype, "_", res_name_sanitized, "EPIi_DesignCompare_volcanoplot.png")
    ggsave(volcano_filename, plot = volcanoplot, width = 5, height = 7)
    
    resOrdered <- res.shrink[order(res.shrink$padj), ]
    
    filename <- paste( res_name_sanitized, "EPIi_DesignCompare.csv", sep = "_")
    write.csv(as.data.frame(resOrdered), file = filename)
    
    significant_results <- subset(resOrdered, padj <= 0.05)
    significant_results<- significant_results[order(significant_results$log2FoldChange), ]
    
    filename <- paste( res_name_sanitized, "EPIi_DesignCompare_significantONLY_sortedFC.csv", sep = "_")
    write.csv(as.data.frame(significant_results), file = filename)
    
    gene_list <- significant_results$log2FoldChange
    names(gene_list) <- rownames(significant_results)
    
    gene_list <- sort(gene_list, decreasing = TRUE)
    ensg_list = gene_list
    
    names(gene_list) <- sapply(strsplit(names(gene_list), "_"), choose_id)
    gene_list <- gene_list[!startsWith(names(gene_list), 'ENSG')]
    
    names(ensg_list) <- sapply(strsplit(names(ensg_list), "_"),choose_ENSG_id)
    print(ensg_list)
    
    if (length(gene_list) == 0) {
      message("No significant genes found for EPIi_", res_name_sanitized)
      next
    }
    
    try({
      
      uplist = names(gene_list[gene_list>0])
      downlist = names(gene_list[gene_list<0])
      UPenriched <- enrichr(uplist, dbs)
      Sys.sleep(5)
      DOWNenriched = enrichr(downlist, dbs)
      Sys.sleep(5)
      ALLenriched = enrichr(names(gene_list), dbs)
      Sys.sleep(5)
      for (DB_name in dbs) {
        sig_terms_UP = UPenriched[[DB_name]][UPenriched[[DB_name]]$Adjusted.P.value<0.05,]
        filename <- paste( res_name_sanitized, DB_name, "EPIi_UP_enrichr_significantONLY.csv", sep = "_")
        write.csv(sig_terms_UP, filename)
        
        Sys.sleep(1)
        sig_terms_down = DOWNenriched[[DB_name]][DOWNenriched[[DB_name]]$Adjusted.P.value<0.05,]
        filename <- paste(res_name_sanitized, DB_name, "EPIi_DOWN_enrichr_significantONLY.csv", sep = "_")
        write.csv(sig_terms_down, filename)
        
        Sys.sleep(1)
        sig_terms_ALL = ALLenriched[[DB_name]][ALLenriched[[DB_name]]$Adjusted.P.value<0.05,]
        filename <- paste( res_name_sanitized, DB_name, "EPIi_ALL_enrichr_significantONLY.csv", sep = "_")
        write.csv(sig_terms_ALL, filename)
        Sys.sleep(1)
      }
      
      m_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>%
        dplyr::select(gs_description, ensembl_gene )
      C7_GSEA <- GSEA(ensg_list, minGSSize = 5, pAdjustMethod = "BH", pvalueCutoff = 0.05, TERM2GENE = m_t2g)
      gse_df <- as.data.frame(C7_GSEA)
      gsea_csv_filename <- paste( donor, celltype, res_name_sanitized, "C7_GSEAResults.csv", sep = "_")
      write.csv(gse_df, file = gsea_csv_filename)
      
      m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>%
        dplyr::select(gs_description, ensembl_gene )
      C5_GSEA <- GSEA(ensg_list, minGSSize = 5, pAdjustMethod = "BH", pvalueCutoff = 0.05, TERM2GENE = m_t2g)
      gse_df <- as.data.frame(C5_GSEA)
      gsea_csv_filename <- paste( donor, celltype, res_name_sanitized, "C5_GSEAResults.csv", sep = "_")
      write.csv(gse_df, file = gsea_csv_filename)
      
      gse <- gseGO(geneList = ensg_list,
                   ont = "ALL",
                   keyType = "ENSEMBL",
                   nPerm = 10000,
                   minGSSize = 3,
                   maxGSSize = 800,
                   pvalueCutoff = 0.05,
                   verbose = TRUE,
                   OrgDb = organism,
                   pAdjustMethod = "BH")
      gse_df <- as.data.frame(gse)
      gsea_csv_filename <- paste(donor, celltype, res_name_sanitized, "GSEGO_GSEAResults.csv", sep = "_")
      write.csv(gse_df, file = gsea_csv_filename)
      
    }, silent = TRUE)
    
  }
}
