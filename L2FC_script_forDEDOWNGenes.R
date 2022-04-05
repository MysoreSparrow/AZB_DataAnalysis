library(tidyverse)
library(EnhancedVolcano)
library(data.table)
library(gplots)
library(RColorBrewer)

# Input the frame that needs to be worked on
DownGeneFrame_Core50 <- read.csv("D:/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/DownGeneFrame_Core50.csv")
#binded_frame <- read.csv("D:/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/Data_Alina_MasterTable_13strains.csv")

by_gene_forlog2FC <- DownGeneFrame_Core50 %>% group_by(symbol) %>% 
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


merged_L2fc1 <- merge(data.frame(log2FC_476, row.names=NULL), data.frame(log2FC_754, row.names=NULL), by = 1, all = TRUE)
merged_L2fc2 <- merge(data.frame(log2FC_755, row.names=NULL), data.frame(log2FC_757, row.names=NULL), by = 1, all = TRUE)
merged_L2fc3 <- merge(data.frame(log2FC_758, row.names=NULL), data.frame(log2FC_760, row.names=NULL), by = 1, all = TRUE)
merged_L2fc4 <- merge(data.frame(log2FC_761, row.names=NULL), data.frame(log2FC_762, row.names=NULL), by = 1, all = TRUE)
merged_L2fc5 <- merge(data.frame(log2FC_763, row.names=NULL), data.frame(log2FC_764, row.names=NULL), by = 1, all = TRUE)
merged_L2fc6 <- merge(data.frame(log2FC_765, row.names=NULL), data.frame(log2FC_766, row.names=NULL), by = 1, all = TRUE)
merged_L2fc7 <- merge(data.frame(log2FC_768, row.names=NULL), data.frame(log2FC_769, row.names=NULL), by = 1, all = TRUE)


merged_L2fc_A1 <- merge(data.frame(merged_L2fc1, row.names=NULL),
                        data.frame(merged_L2fc2, row.names=NULL),
                        by = 1, all = TRUE)
merged_L2fc_A2 <- merge(data.frame(merged_L2fc3, row.names=NULL), 
                        data.frame(merged_L2fc4, row.names=NULL), 
                        by = 1, all = TRUE)
merged_L2fc_A3 <- merge(data.frame(merged_L2fc5, row.names=NULL), 
                        data.frame(merged_L2fc6, row.names=NULL), 
                        by = 1, all = TRUE)
merged_L2fc_B1 <- merge(data.frame(merged_L2fc_A1, row.names=NULL), 
                        data.frame(merged_L2fc_A2, row.names=NULL), 
                        by = 1, all = TRUE)
merged_L2fc_B2 <- merge(data.frame(merged_L2fc_A3, row.names=NULL), 
                        data.frame(merged_L2fc7, row.names=NULL), 
                        by = 1, all = TRUE)
merged_L2FC <- merge(data.frame(merged_L2fc_B1, row.names=NULL), 
                     data.frame(merged_L2fc_B2, row.names=NULL), 
                     by = 1, all = TRUE)


# change the type of rest of columns except for symbol from chr to double
merged_L2FC <- merged_L2FC %>% mutate_at(c(2:15), as.double)
# ensure that there are no duplicate rows
merged_L2FC <- merged_L2FC %>% distinct(symbol,.keep_all=TRUE)

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
ZS_DownGenes_table <- data.frame(GeneName, zs_476,zs_754,zs_755,zs_757,zs_758,
                               zs_760, zs_761,zs_762,zs_763,zs_764,zs_765,
                               zs_766,zs_768,zs_769)

ZS_DownGenes_table <- data.frame(ZS_DownGenes_table, row.names = 1)

#write.csv(ZS_DownGenes_table,"D:/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/ZS_DownGenes_table.csv","w")

# MEAN OF MINIMUM WAS 1.1796611

#ZS_DownGenes_table <- ZS_DownGenes_table %>% replace(is.na(.), -1.1796611)

geneExp_matrix_DE_DownGene <- as.matrix(ZS_DownGenes_table)

heatmap(geneExp_matrix_DE_DownGene)

tiff("D:/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/heatmaps/h1_DownGenes.tif", compression = "lzw")
heatmap(geneExp_matrix_DE_DownGene)
dev.off()