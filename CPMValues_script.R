## Collecting the CPMr1 and CPMr23 values , gene by gene

binded_frame <- read.csv("D:/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/Data_Alina_MasterTable_13strains.csv")

by_gene <- binded_frame %>% 
           group_by(symbol) %>% 
           summarise(Strain_name, CPM_infected_R1, 
           CPM_infected_R2) %>% ungroup()
head(by_gene)

### bygene_number
bygene_476 <- by_gene %>% 
           group_by(symbol, Strain_name) %>% 
           filter(Strain_name == "S476") %>%
           summarise(cpm_476_r1 = paste0(CPM_infected_R1, collapse = ","),
                     cpm_476_r2 = paste0(CPM_infected_R2, collapse = ","))%>% 
           select(-Strain_name) 
bygene_755 <- by_gene %>% 
           group_by(symbol, Strain_name) %>% 
           filter(Strain_name == "S755") %>%
           summarise(cpm_755_r1 = paste0(CPM_infected_R1, collapse = ","),
                     cpm_755_r2 = paste0(CPM_infected_R2, collapse = ",")) %>% 
           select(-Strain_name)
bygene_754 <- by_gene %>% 
           group_by(symbol, Strain_name) %>% 
           filter(Strain_name == "S754") %>%
           summarise(cpm_754_r1 = paste0(CPM_infected_R1, collapse = ","),
                     cpm_754_r2 = paste0(CPM_infected_R2, collapse = ",")) %>% 
           select(-Strain_name)
bygene_757 <- by_gene %>% 
           group_by(symbol, Strain_name) %>% 
           filter(Strain_name == "S757") %>%
           summarise(cpm_757_r1 = paste0(CPM_infected_R1, collapse = ","),
                     cpm_757_r2 = paste0(CPM_infected_R2, collapse = ",")) %>% 
           select(-Strain_name)
bygene_758 <- by_gene %>% 
           group_by(symbol, Strain_name) %>% 
           filter(Strain_name == "S758") %>%
           summarise(cpm_758_r1 = paste0(CPM_infected_R1, collapse = ","),
                     cpm_758_r2 = paste0(CPM_infected_R2, collapse = ",")) %>% 
           select(-Strain_name)

bygene_760 <- by_gene %>%
           group_by(symbol, Strain_name) %>% 
           filter(Strain_name == "S760") %>%
           summarise(cpm_760_r1 = paste0(CPM_infected_R1, collapse = ","),
                     cpm_760_r2 = paste0(CPM_infected_R2, collapse = ",")) %>% 
           select(-Strain_name)
bygene_761 <- by_gene %>% 
           group_by(symbol, Strain_name) %>% 
           filter(Strain_name == "S761") %>%
           summarise(cpm_761_r1 = paste0(CPM_infected_R1, collapse = ","),
                     cpm_761_r2 = paste0(CPM_infected_R2, collapse = ",")) %>% 
           select(-Strain_name)
bygene_762 <- by_gene %>% 
           group_by(symbol, Strain_name) %>% 
           filter(Strain_name == "S762") %>%
           summarise(cpm_762_r1 = paste0(CPM_infected_R1, collapse = ","),
                     cpm_762_r2 = paste0(CPM_infected_R2, collapse = ",")) %>% 
           select(-Strain_name)
bygene_763 <- by_gene %>% 
           group_by(symbol, Strain_name) %>% 
           filter(Strain_name == "S763") %>%
           summarise(cpm_763_r1 = paste0(CPM_infected_R1, collapse = ","),
                     cpm_763_r2 = paste0(CPM_infected_R2, collapse = ",")) %>% 
           select(-Strain_name)
bygene_764 <- by_gene %>% 
           group_by(symbol, Strain_name) %>% 
           filter(Strain_name == "S764") %>%
           summarise(cpm_764_r1 = paste0(CPM_infected_R1, collapse = ","),
                     cpm_764_r2 = paste0(CPM_infected_R2, collapse = ",")) %>% 
           select(-Strain_name)
bygene_765 <- by_gene %>% 
           group_by(symbol, Strain_name) %>% 
           filter(Strain_name == "S765") %>%
           summarise(cpm_765_r1 = paste0(CPM_infected_R1, collapse = ","),
                     cpm_765_r2 = paste0(CPM_infected_R2, collapse = ",")) %>% 
           select(-Strain_name)
bygene_766 <- by_gene %>% 
           group_by(symbol, Strain_name) %>% 
           filter(Strain_name == "S766") %>%
           summarise(cpm_766_r1 = paste0(CPM_infected_R1, collapse = ","),
                     cpm_766_r2 = paste0(CPM_infected_R2, collapse = ",")) %>% 
           select(-Strain_name)
bygene_768 <- by_gene %>% 
           group_by(symbol, Strain_name) %>% 
           filter(Strain_name == "S768") %>%
           summarise(cpm_768_r1 = paste0(CPM_infected_R1, collapse = ","),
                     cpm_768_r2 = paste0(CPM_infected_R2, collapse = ",")) %>% 
           select(-Strain_name)
bygene_769 <- by_gene %>% 
           group_by(symbol, Strain_name) %>% 
           filter(Strain_name == "S769") %>%
           summarise(cpm_769_r1 = paste0(CPM_infected_R1, collapse = ","),
                     cpm_769_r2 = paste0(CPM_infected_R2, collapse = ",")) %>% 
           select(-Strain_name)


# Now lets try to merge the dataframe for respective genes, row-wise

# Merging two df at a time.
merged_df1 <- merge(data.frame(bygene_476, row.names=NULL), data.frame(bygene_754, row.names=NULL), by = 1, all = TRUE)
merged_df2 <- merge(data.frame(bygene_755, row.names=NULL), data.frame(bygene_757, row.names=NULL), by = 1, all = TRUE)

merged_df3 <- merge(data.frame(bygene_758, row.names=NULL), data.frame(bygene_760, row.names=NULL), by = 1, all = TRUE)
merged_df4 <- merge(data.frame(bygene_761, row.names=NULL), data.frame(bygene_762, row.names=NULL), by = 1, all = TRUE)

merged_df5 <- merge(data.frame(bygene_763, row.names=NULL), data.frame(bygene_764, row.names=NULL), by = 1, all = TRUE)
merged_df6 <- merge(data.frame(bygene_765, row.names=NULL), data.frame(bygene_766, row.names=NULL), by = 1, all = TRUE)

merged_df7 <- merge(data.frame(bygene_768, row.names=NULL), data.frame(bygene_769, row.names=NULL), by = 1, all = TRUE)


CPM_Merged12 <- merge(
                   data.frame(merged_df1, row.names=NULL), 
                   data.frame(merged_df2, row.names=NULL), all = TRUE)
CPM_Merged12

CPM_Merged34 <- merge(
                   data.frame(merged_df3, row.names=NULL), 
                   data.frame(merged_df4, row.names=NULL), all = TRUE)
CPM_Merged34

CPM_Merged56 <- merge(
                   data.frame(merged_df5, row.names=NULL), 
                   data.frame(merged_df6, row.names=NULL), all = TRUE)
CPM_Merged56

CPM_Merged567 <- merge(
                   data.frame(CPM_Merged56, row.names=NULL), 
                   data.frame(merged_df7, row.names=NULL), all = TRUE)
CPM_Merged567

CPM_MergedTable_1234 <- merge(
                   data.frame(CPM_Merged12, row.names=NULL), 
                   data.frame(CPM_Merged34, row.names=NULL), all = TRUE)
CPM_MergedTable_1234


# Finally open the complete merged table

CPM_MergedTable_ByGene <- merge(
                   data.frame(CPM_MergedTable_1234, row.names=NULL), 
                   data.frame(CPM_Merged567, row.names=NULL), all = TRUE)
CPM_MergedTable_ByGene

write.csv(CPM_MergedTable_ByGene,"D:/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/CPM_MergedTable_ByGene.csv")
