library(tidyverse)
library(EnhancedVolcano)
library(data.table)
library(gplots)
library(RColorBrewer)

# Input the frame that needs to be worked on
UpGeneFrame_Core50 <- read.csv("D:/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/UpGeneFrame_Core50.csv")
#binded_frame <- read.csv("D:/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/Data_Alina_MasterTable_13strains.csv")

by_gene_forlog2FC <- UpGeneFrame_Core50 %>% group_by(symbol) %>% 
  summarise(log2FoldChange, Strain_name) %>% ungroup()
head(by_gene_forlog2FC)

log2FC_476 <- by_gene_forlog2FC %>% 
  group_by(symbol) %>% 
  filter(Strain_name == "S476") %>%
  summarise(log2FC_476 = paste0(log2FoldChange)) 
log2FC_754 <- by_gene_forlog2FC %>% 
  group_by(symbol) %>% 
  filter(Strain_name == "S754") %>%
  summarise(log2FC_754 = paste0(log2FoldChange))
log2FC_755 <- by_gene_forlog2FC %>% 
  group_by(symbol) %>% 
  filter(Strain_name == "S755") %>%
  summarise(log2FC_755 = paste0(log2FoldChange))
log2FC_757 <- by_gene_forlog2FC %>% 
  group_by(symbol) %>% 
  filter(Strain_name == "S757") %>%
  summarise(log2FC_757 = paste0(log2FoldChange))
log2FC_758 <- by_gene_forlog2FC %>% 
  group_by(symbol) %>% 
  filter(Strain_name == "S758") %>%
  summarise(log2FC_758 = paste0(log2FoldChange))
log2FC_760 <- by_gene_forlog2FC %>% 
  group_by(symbol) %>% 
  filter(Strain_name == "S760") %>%
  summarise(log2FC_760 = paste0(log2FoldChange))
log2FC_761 <- by_gene_forlog2FC %>% 
  group_by(symbol) %>% 
  filter(Strain_name == "S761") %>%
  summarise(log2FC_761 = paste0(log2FoldChange))
log2FC_762 <- by_gene_forlog2FC %>% 
  group_by(symbol) %>% 
  filter(Strain_name == "S762") %>%
  summarise(log2FC_762 = paste0(log2FoldChange))
log2FC_763 <- by_gene_forlog2FC %>% 
  group_by(symbol) %>% 
  filter(Strain_name == "S763") %>%
  summarise(log2FC_763 = paste0(log2FoldChange))
log2FC_764 <- by_gene_forlog2FC %>% 
  group_by(symbol) %>% 
  filter(Strain_name == "S764") %>%
  summarise(log2FC_764 = paste0(log2FoldChange))
log2FC_765 <- by_gene_forlog2FC %>% 
  group_by(symbol) %>% 
  filter(Strain_name == "S765") %>%
  summarise(log2FC_765 = paste0(log2FoldChange))
log2FC_766 <- by_gene_forlog2FC %>% 
  group_by(symbol) %>% 
  filter(Strain_name == "S766") %>%
  summarise(log2FC_766 = paste0(log2FoldChange))
log2FC_768 <- by_gene_forlog2FC %>% 
  group_by(symbol) %>% 
  filter(Strain_name == "S768") %>%
  summarise(log2FC_768 = paste0(log2FoldChange))
log2FC_769 <- by_gene_forlog2FC %>% 
  group_by(symbol) %>% 
  filter(Strain_name == "S769") %>%
  summarise(log2FC_769 = paste0(log2FoldChange))


df_list_UPGeneDF <- list(log2FC_476, log2FC_754, log2FC_755, log2FC_757,log2FC_758,
                log2FC_760, log2FC_761, log2FC_762, log2FC_763,log2FC_764,
                log2FC_765,log2FC_766,log2FC_768,log2FC_769)

# Merge all the df into a single dataframe as per individual gene - with all values per gene being in one row
merged_L2FC <- Reduce(function(d1, d2) merge(d1, d2, by = "symbol", all.x = TRUE, 
                                             all.y = FALSE), df_list_UPGeneDF)

# Determining the Zscore tables for UP Genes

zs_476 = scale(merged_L2FC$log2FC_476, center = TRUE, scale = TRUE)
zs_754 = scale(merged_L2FC$log2FC_754, center = TRUE, scale = TRUE)
zs_755 = scale(merged_L2FC$log2FC_755, center = TRUE, scale = TRUE)
zs_757 = scale(merged_L2FC$log2FC_757, center = TRUE, scale = TRUE)
zs_758 = scale(merged_L2FC$log2FC_758, center = TRUE, scale = TRUE)
zs_760 = scale(merged_L2FC$log2FC_760, center = TRUE, scale = TRUE)
zs_761 = scale(merged_L2FC$log2FC_761, center = TRUE, scale = TRUE)
zs_762 = scale(merged_L2FC$log2FC_762, center = TRUE, scale = TRUE)
zs_763 = scale(merged_L2FC$log2FC_763, center = TRUE, scale = TRUE)
zs_764 = scale(merged_L2FC$log2FC_764, center = TRUE, scale = TRUE)
zs_765 = scale(merged_L2FC$log2FC_765, center = TRUE, scale = TRUE)
zs_766 = scale(merged_L2FC$log2FC_766, center = TRUE, scale = TRUE)
zs_768 = scale(merged_L2FC$log2FC_768, center = TRUE, scale = TRUE)
zs_769 = scale(merged_L2FC$log2FC_769, center = TRUE, scale = TRUE)

GeneName <- merged_L2FC$symbol #obtaining genenames from earlier df
ZS_UPGenes_table <- data.frame(GeneName, zs_476,zs_754,zs_755,zs_757,zs_758,
                               zs_760, zs_761,zs_762,zs_763,zs_764,zs_765,
                               zs_766,zs_768,zs_769)

ZS_UPGenes_table <- data.frame(ZS_UPGenes_table, row.names = 1)

#write.csv(ZS_UPGenes_table,"D:/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/ZS_UPGenes_table.csv","w")

# MEAN OF MINIMUM WAS 1.1796611

#ZS_UPGenes_table <- ZS_UPGenes_table %>% replace(is.na(.), -1.1796611)

geneExp_matrix_DE_UpGene <- as.matrix(ZS_UPGenes_table)

heatmap(geneExp_matrix_DE_UpGene)

tiff("D:/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/heatmaps/h1_UPGenes.tif", compression = "lzw")
heatmap(geneExp_matrix_DE_UpGene)
dev.off()