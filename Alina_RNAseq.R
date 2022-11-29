# Author: KeshavPrasadGubbi

Package_List <- c(
  "dplyr", "DESeq2", "pheatmap", "PoiClaClu", "RColorBrewer", "vsn", "EnhancedVolcano", "gplots",
  "org.Mm.eg.db", "stringr", "genefilter", "tidyverse", "AnnotationDbi", "ComplexHeatmap", "DOSE",
  "clusterProfiler", "ggrepel", "GO.db", "GOstats", "gage", "gageData", "GOSemSim", "enrichplot",
  "ggnewscale", "glue", "ggupset", "FactoMineR", "factoextra", "here", "tibble"
)
not_installed <- Package_List[!(Package_List %in% installed.packages()[, "Package"])] # Extract not installed packages
if (length(not_installed)) install.packages(not_installed) # Install the uninstalled packages
invisible(lapply(Package_List, suppressPackageStartupMessages(library), character.only = TRUE))

# File Path Declarations
here::i_am(path = "Alina_RNAseq.R")
paste0(here())

# Also Create a comparison Variable: That Could be used later for all other comparison titles using a
# glue Variable. Define the Comparison and also create the folder for saving all plots and results to be
# saved as per the comparison

Comparison <- "BL6_InfectedVsControl"
# Determine the Comparison Condition: Comment one of them out based on the comparison you are trying to run.
Comparison_Condition <- "condition"

# Folder Paths for Different Comparisons
Comparison_path <- file.path(here(), glue("{Comparison}"))

if (!dir.exists(here(Comparison_path))) { dir.create(here(Comparison_path))} else { print("Dir already exists")}
paste0(Comparison_path)

# for Example: UninfectedVSInfected is the volcano plot notation: for Deseq2 results,
# the infected is Numerator.and uninfected is Denominator. which will result in DESEQ contrast like (condition, Infected, Uninfected) and hence the infected part will come on the right side of the volcano plot.


## Creating metadata for the DGE Analysis

# Read the csv file and change the column name. the samples.csv is a list of sample names, ie, the names of bam files.

sample_ID <- read.csv(file.path(here(), "/samples.csv"))
condition <- c(
  "Infected", "Infected", "Infected", "Infected", "Infected", "Infected",
  "Infected", "Infected", "Infected", "Infected", "Infected", "Infected",
  "Infected", "Infected", "control", "control"
)
coldata <- data.frame(sample_ID, condition)
colnames(coldata) <- c("Sample_Name", "condition") # change name of one of the columns
# The metadata can be found in a df called coldata!

# Tidying up the names for plots later! First from coldata.

# tidying up the names of samples in both columns that list of samples
coldata$Sample_Name <- str_remove_all(coldata$Sample_Name, pattern = "run6_trimmed_|_.bam|_S\\d\\d|_S\\d")
coldata$condition <- as.factor(coldata$condition)

# Changing the names of samples (as per Alina)
coldata[coldata == "476_R1"] <- "T"
coldata[coldata == "754_R1"] <- "S54"
coldata[coldata == "755_R1"] <- "S55"
coldata[coldata == "757_R1"] <- "L57"
coldata[coldata == "758_R1"] <- "A58"
coldata[coldata == "760_R1"] <- "L60"
coldata[coldata == "761_R1"] <- "S61"
coldata[coldata == "762_R1"] <- "A62"
coldata[coldata == "763_R1"] <- "L63"
coldata[coldata == "764_R1"] <- "A64"
coldata[coldata == "765_R1"] <- "S65"
coldata[coldata == "766_R1"] <- "L66"
coldata[coldata == "768_R1"] <- "A68"
coldata[coldata == "769_R1"] <- "L69"
coldata[coldata == "Ctrl1_R1"] <- "C1"
coldata[coldata == "Ctrl2_R2"] <- "C2"
# convert column1 with sample names to row.names of coldata
rownames(coldata) <- coldata$Sample_Name

# Adding the groupings by Alina for further Metadata Information

coldata$Epithelial_response <- c(
  "LowInducer", "LowInducer", "HighInducer",
  "HighInducer", "LowInducer", "LowInducer",
  "HighInducer", "HighInducer", "LowInducer",
  "HighInducer", "LowInducer", "HighInducer",
  "LowInducer", "LowInducer", "NR", "NR"
)
coldata$clinical_outcome <- c(
  "symptomatic", "symptomatic", "symptomatic",
  "Lethal", "asymptomatic", "Lethal", "symptomatic",
  "asymptomatic", "Lethal", "symptomatic", "symptomatic",
  "Lethal", "asymptomatic", "Lethal", "NR", "NR"
)
coldata$microcolonies <- c(
  "Low", "Low", "Low", "High", "Low", "Low",
  "High", "High", "Low", "High", "Low", "High", "Low",
  "Low", "NR", "NR"
)
coldata$ER_microcolonies <- c(
  "LI_LM", "LI_LM", "HI_LM", "HI_HM", "LI_LM", "LI_LM",
  "HI_HM", "HI_HM", "LI_LM", "HI_HM", "LI_LM", "HI_HM",
  "LI_LM", "LI_LM", "NR", "NR"
)
coldata$phylogenomic_lineage <- c(
  "EPEC1", "EPEC10", "EPEC9", "EPEC9", "NC", "EPEC5",
  "EPEC8", "NC", "EPEC7", "NC", "EPEC2", "EPEC9",
  "EPEC2", "EPEC2", "NR", "NR"
)
coldata$phylogroup <- c(
  "B2", "A", "B2", "B2", "B1", "A", "B2", "B2", "B1", "B2", "B1",
  "B2", "B1", "B2", "NR", "NR"
)
coldata$Intimin_Type <- c(
  "alpha", "ND", "lambda", "lambda", "epsilon", "epsilon",
  "mu", "lambda", "beta", "kappa", "beta", "alpha", "beta",
  "beta", "NR", "NR"
)

# Lets Deal with the Countmatrix
# Readin  countsmatrix

countsmatrix <- as.data.frame(read.csv(file.path(here(), "/newcounts.csv")))
names(countsmatrix)[1] <- "EnsemblID" # change name of 1st column

## Removal of Gender Genes from ENSEMBL ID itself
# Filter out the other 5 known gender genes as well from other RNAseq Projects.
countsmatrix <- countsmatrix %>% filter(countsmatrix$EnsemblID != "ENSMUSG00000086503",
                                        countsmatrix$EnsemblID != "ENSMUSG00000097571",
                                        countsmatrix$EnsemblID != "ENSMUSG00000086370",
                                        countsmatrix$EnsemblID != "ENSMUSG00000031329",
                                        countsmatrix$EnsemblID != "ENSMUSG00000030057")

# Annotating and Exporting ENSEMBL ID into Gene Symbols
## Adding genes annotated from ENSEMBL ID to Gene symbols and ENTREZ Id to countsmatrix table.
symbols <- as.data.frame(mapIds(org.Mm.eg.db,
                  keys = countsmatrix$EnsemblID, # mapping ENSEMBL to Gene Symbol
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")) %>%
  rownames_to_column(var = "EID") %>% #move the ensemblID in rownames to separate column called EID
  drop_na() %>% # drop Na rows so that it reduces size of matrix. Na values arise due to 1:many mapping of ensembl.
  rename(genename = 2) # chnage the name of 2nd column to genename

countsmatrix <- countsmatrix %>%
  filter(countsmatrix$EnsemblID %in% symbols$EID) %>% # keep the genes in countmatrix that have gene names in symbols, based on matching ensemblID columns.
  mutate(genename = symbols$genename)

# Removing the duplicated genes & then these genes put into rownames for countsmatrix and drop EnsemblId column
countsmatrix <- countsmatrix[!duplicated(countsmatrix$genename), ] %>% remove_rownames() %>%
  column_to_rownames(var = "genename") %>% select(-EnsemblID) %>% as.matrix()
# the elements from Sample_Name from coldata must the the colnames of countsmatrix
colnames(countsmatrix) <- coldata$Sample_Name
# Changing countsmatrix into Matrix of numeric values so that only numeric values are present in it as an input of DESEq Object.
class(countsmatrix) <- "numeric"



