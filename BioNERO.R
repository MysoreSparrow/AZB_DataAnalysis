

Package_List <- c(
  "dplyr",
  "DESeq2",
  "pheatmap",
  "PoiClaClu",
  "RColorBrewer",
  "vsn",
  "EnhancedVolcano",
  "gplots",
  "org.Mm.eg.db",
  "stringr",
  "genefilter",
  "tidyverse",
  "AnnotationDbi",
  "ComplexHeatmap",
  "DOSE",
  "clusterProfiler",
  "ggrepel",
  "GO.db",
  "GOstats",
  "gage",
  "gageData",
  "GOSemSim",
  "enrichplot",
  "ggnewscale",
  "glue",
  "ggupset",
  "FactoMineR",
  "factoextra",
  "here",
  "tibble",
  "edgeR",
  "BioNERO",
  "readxl",
  "WGCNA"
)
not_installed <-
  Package_List[!(Package_List %in% installed.packages()[, "Package"])] # Extract not installed packages
if (length(not_installed))
  install.packages(not_installed) # Install the uninstalled packages
invisible(lapply(
  Package_List,
  suppressPackageStartupMessages(library),
  character.only = TRUE
))
set.seed(123)

# File Path Declarations
here::i_am(path = "BioNERO.R")
paste0(here())

# Metadata for the Analysis

## Creating metadata for the DGE Analysis
# Read the csv file and change the column name. the samples.csv is a list of sample names, ie, the names of bam files.

sample_ID <- read.csv(file.path(here(), "/samples.csv"))
condition <- c(
  "Infected",
  "Infected",
  "Infected",
  "Infected",
  "Infected",
  "Infected",
  "Infected",
  "Infected",
  "Infected",
  "Infected",
  "Infected",
  "Infected",
  "Infected",
  "Infected",
  "control",
  "control"
)
coldata <- data.frame(sample_ID, condition)
colnames(coldata) <-
  c("Sample_Name", "condition") # change name of one of the columns
# The metadata can be found in a df called coldata!

# Tidying up the names for plots later! First from coldata.

# tidying up the names of samples in both columns that list of samples
coldata$Sample_Name <-
  str_remove_all(coldata$Sample_Name, pattern = "run6_trimmed_|_.bam|_S\\d\\d|_S\\d")
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
  "LowInducer",
  "LowInducer",
  "HighInducer",
  "HighInducer",
  "LowInducer",
  "LowInducer",
  "HighInducer",
  "HighInducer",
  "LowInducer",
  "HighInducer",
  "LowInducer",
  "HighInducer",
  "LowInducer",
  "LowInducer",
  "NR",
  "NR"
)
# coldata$clinical_outcome <- c(
#   "symptomatic", "symptomatic", "symptomatic",
#   "Lethal", "asymptomatic", "Lethal", "symptomatic",
#   "asymptomatic", "Lethal", "symptomatic", "symptomatic",
#   "Lethal", "asymptomatic", "Lethal", "NR", "NR"
# )
#

coldata <- coldata[1:14,] # Remove C1 and C2
coldata_bionero <- coldata
coldata_bionero <- coldata_bionero[, -c(1, 2)]
head(coldata_bionero)

# ## Read-in the Effectors
effector <- as.data.frame(read_excel((
  file.path(here(), "allStrains_PresentAbsent_effectorlist_forKeshav.xlsx")),
  col_types = c(
    "text", "numeric", "numeric", "numeric", "numeric",
    "numeric", "numeric", "numeric", "numeric", "numeric",
    "numeric", "numeric", "numeric", "numeric", "numeric"
  )
))

rownames(effector) <- effector[,1] # move effectors to rownames
effector <- subset(effector, select = - Effector) #remove effector symbol column
effector <- t(effector)
head(effector)

# Replace 0, 1 , 2
effector[effector == "0"] <- "Absent"
effector[effector == "1"] <- "Present"
effector[effector == "2"] <- "Double"
#
#
# # Merge coldata and effectors into single DF
coldata_merged <- merge(coldata, effector, by = 'row.names', all = TRUE)
rownames(coldata_merged) <- coldata_merged[,1] # move effectors to rownames
coldata_merged <- subset(coldata_merged, select = - Row.names) #remove effector symbol column
coldata_merged <- coldata_merged[, -c(1, 2)]
head(coldata_merged)

# Lets Deal with the Countmatrix
# Readin  countsmatrix
countsmatrix <-
  as.data.frame(read.csv(file.path(here(), "/newcounts.csv")))
names(countsmatrix)[1] <- "EnsemblID" # change name of 1st column

## Removal of Gender Genes from ENSEMBL ID itself
# Filter out the other 5 known gender genes as well from other RNAseq Projects.
countsmatrix <- countsmatrix %>% filter(
  countsmatrix$EnsemblID != "ENSMUSG00000086503",
  countsmatrix$EnsemblID != "ENSMUSG00000097571",
  countsmatrix$EnsemblID != "ENSMUSG00000086370",
  countsmatrix$EnsemblID != "ENSMUSG00000031329",
  countsmatrix$EnsemblID != "ENSMUSG00000030057"
)

# Annotating and Exporting ENSEMBL ID into Gene Symbols
## Adding genes annotated from ENSEMBL ID to Gene symbols and ENTREZ Id to countsmatrix table
symbols <- as.data.frame(
  mapIds(
    org.Mm.eg.db,
    keys = countsmatrix$EnsemblID,
    # mapping ENSEMBL to Gene Symbol
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
) %>%
  rownames_to_column(var = "EID") %>% # move the ensemblID in rownames to separate column called EID
  drop_na() %>% # drop Na rows so that it reduces size of matrix. Na values arise due to 1:many mapping of ensembl.
  rename(genename = 2) # change the name of 2nd column to genename

countsmatrix <- countsmatrix %>%
  filter(countsmatrix$EnsemblID %in% symbols$EID) %>% # keep the genes in countmatrix that have gene names in symbols, based on matching ensemblID columns.
  mutate(genename = symbols$genename)

# Removing the duplicated genes & then these genes put into rownames for countsmatrix and drop EnsemblId column
countsmatrix <-
  countsmatrix[!duplicated(countsmatrix$genename),] %>%
  remove_rownames() %>%
  column_to_rownames(var = "genename") %>%
  select(-EnsemblID) %>%
  as.matrix()

# drop the Control Columns
countsmatrix <- countsmatrix[, 1:14]

# the elements from Sample_Name from coldata must the the colnames of countsmatrix
colnames(countsmatrix) <- rownames(coldata)
# Changing countsmatrix into Matrix of numeric values so that only numeric values are present in it as an input of DESEq Object.
class(countsmatrix) <- "numeric"
head(countsmatrix)

# countsmatrix <- t(countsmatrix)
# genomic_idx <- match(rownames(coldata_merged), colnames(countsmatrix))
# genomic_idx
# countsmatrix <- countsmatrix[ , genomic_idx]

# create SummarizedExperiment object
se = SummarizedExperiment(list(counts = countsmatrix), colData = (coldata_bionero))
# se@metadata <- as.data.frame(effector) # Adding the effectors as an metadata to se object

# Pre-Processing
# Automatic One step Processing
final_exp <- exp_preprocess(
  se,
  n = 2500,
  Zk_filtering = TRUE,
  variance_filter = TRUE,
  cor_method = "pearson",
  vstransform = TRUE
)
print(final_exp)

# Exploratory Data Analysis
(Sample_Corr_Heatmap <-
    plot_heatmap(exp = final_exp, type = "samplecor"))
# col_metadata = c(coldata$Sample_Name, coldata$Epithelial_response)))

(Expression_Heatmap = plot_heatmap(final_exp, type = "expr"))

# final_exp1 <- final_exp
# final_exp1@metadata <- as.data.frame(effector)
# (plot_PCA(final_exp1, metadata = colData(final_exp1)
#           )
#   )

# Gene CoExpression Network Inference

sft <-
  SFT_fit(final_exp, net_type = "signed hybrid", cor_method = "pearson")
sft$power
power <- sft$power
sft$plot

net <- exp2gcn(
  final_exp,
  net_type = "signed hybrid",
  SFTpower = power,
  cor_method = "pearson"
)
names(net)

# Dendro and colors
plot_dendro_and_colors(net)

# Eigengene networks
plot_eigengene_network(net)

### Number of Genes in every Module
plot_ngenes_per_module(net)

# Gene coexpression network analysis
#
# Make rows in coldata_merged with same order of row names as coldata
# row_order <- match(rownames(coldata), rownames(coldata_merged))
# row_order
# coldata_merged <- coldata_merged[row_order,]
# coldata_merged <- coldata_merged[, 1:2]
# coldata_merged <- coldata_merged[,c(3,2)]
# head(coldata_merged)


# Module Trait Correlations
# Using the SE Object
# MEtrait <- module_trait_cor(
#   exp = final_exp,
#   MEs = net$MEs,
#   cor_method = "pearson",
#   continuous_trait = FALSE,
#   transpose = FALSE
# )
# head(MEtrait)

# MEtrait <- module_trait_cor(
#   exp = final_exp,
#   MEs = net$MEs,
#   cor_method = "pearson",
#   continuous_trait = FALSE,
#   transpose = FALSE
# )
# head(MEtrait)

###########To make it work like a individual column###########
exp_df <- as.data.frame(final_exp@assays@data@listData)
dim(exp_df)
metdata_df <- as.data.frame(final_exp@colData)
# metdata_df <- as.data.frame(coldata_merged$NleA_1)
dim(metdata_df)

# MEtrait_new <- module_trait_cor(
#   exp = exp_df,
#   metadata = metdata_df,
#   MEs = net$MEs,
#   cor_method = "spearman",
#   continuous_trait = FALSE,
#   transpose = FALSE
# )
# head(MEtrait_new)

## Trying with all columns of coldata merged
##
# Create a vector from column 2 to last column
idx <- 3:4 #ncol(coldata_merged)

# Create an empty Dataframe
df = data.frame(matrix(nrow = 66, ncol = 4))

# cor_list <- lapply(idx, function(x) {
#   MEtrait <- module_trait_cor(exp = exp_df,
#                               MEs = net$MEs,
#                               metadata = coldata_merged[, x, drop = FALSE],
#                               cor_method = "pearson")
#   return(MEtrait)
# })

## Gene Significance
(gs_Reg3g <- gene_significance(exp = final_exp,
                        genes = "Reg3g",
                        show_rownames = TRUE))
# head(gs, 3)

# Visualising Module Expression Profile

plot_expression_profile(
  exp = final_exp,
  net = net,
  plot_module = TRUE,
  modulename = "darkturquoise"
)
plot_expression_profile(
  exp = final_exp,
  net = net,
  plot_module = TRUE,
  modulename = "green"
)
plot_expression_profile(
  exp = final_exp,
  net = net,
  plot_module = TRUE,
  modulename = "black"
)

# Enrichment Analysis

# functionalAnalysis <- module_enrichment(net = net,
#                                         background_genes = rownames(final_exp))

# Hub gene identification
hubs <- get_hubs_gcn(final_exp, net)
hubs

# Extracting subgraphs

edges <- get_edge_list(net, module = "green")
head(edges)

# Remove edges based on optimal scale-free topology fit
edges_filtered <-
  get_edge_list(
    net,
    module = "green",
    filter = TRUE,
    method = "pvalue",
    nSamples = ncol(final_exp)
  )

# Network Visualisation
# plot_gcn(
#   edgelist_gcn = edges_filtered,
#   net = net,
#   color_by = "module",
#   hubs = hubs
# )

# plot_gcn(
#   edgelist_gcn = edges_filtered,
#   net = net,
#   color_by = "module",
#   hubs = hubs,
#   show_labels = "tophubs",
#   top_n_hubs = 5,
#   interactive = TRUE,
#   dim_interactive = c(500, 500)
# )


# My own trial for a better module trait correlation
handleSE <- function(exp) {
  if(is(exp, "SummarizedExperiment")) {
    fexp <- SummarizedExperiment::assay(exp)
  } else {
    fexp <- exp
  }
  return(fexp)
}

pval2symbol <- function(matrix) {
  modtraitsymbol <- matrix
  modtraitsymbol[modtraitsymbol < 0.001] <- "***"
  modtraitsymbol[modtraitsymbol >= 0.001 & modtraitsymbol < 0.01] <- "**"
  modtraitsymbol[modtraitsymbol >= 0.01 & modtraitsymbol < 0.05] <- "*"
  modtraitsymbol[modtraitsymbol >= 0.05] <- ""
  return(modtraitsymbol)
}


handle_trait_type <- function(metadata, continuous_trait = FALSE) {
  if(!continuous_trait) {
    sampleinfo <- cbind(Samples = rownames(metadata), metadata)
    tmpdir <- tempdir()
    tmpfile <- tempfile(tmpdir = tmpdir, fileext = "traitmatrix.txt")
    tablesamples <- table(sampleinfo)
    write.table(
      tablesamples, file = tmpfile, quote = FALSE, sep = "\t",
      row.names = TRUE
    )
    trait <- read.csv(tmpfile, header = TRUE, sep = "\t", row.names = 1)
    unlink(tmpfile)
  } else {
    trait <- metadata
  }
  return(trait)
}

module_trait_cor_new <- function(exp, metadata, MEs, cor_method = "spearman",
                             transpose = FALSE, palette = "RdYlBu",
                             continuous_trait = FALSE,
                             cex.lab.x = 0.6, cex.lab.y = 0.6,
                             cex.text = 0.6) {

  if(is(exp, "SummarizedExperiment")) {
    metadata <- as.data.frame(SummarizedExperiment::colData(exp))
  }
  exp <- handleSE(exp)
  metadata <- metadata[colnames(exp), , drop=FALSE]
  trait <- handle_trait_type(metadata, continuous_trait)

  modtraitcor <- cor(as.matrix(MEs), trait, use = "p", method = cor_method)
  nSamples <- ncol(exp)
  modtraitpvalue <- WGCNA::corPvalueStudent(modtraitcor, nSamples)
  modtraitsymbol <- pval2symbol(modtraitpvalue)

  modtrait_long <- reshape2::melt(modtraitcor)
  pvals_long <- reshape2::melt(modtraitpvalue)
  combined_long <- merge(modtrait_long, pvals_long, by = c("Var1", "Var2"))
  colnames(combined_long) <- c("ME", "trait", "cor", "pvalue")
  combined_long$ME <- as.character(combined_long$ME)
  combined_long$trait <- as.character(combined_long$trait)

  textMatrix <- paste(signif(modtraitcor, 2), modtraitsymbol, sep = "")
  dim(textMatrix) <- dim(modtraitcor)

  if(transpose) {
    modtraitcor <- t(modtraitcor)
    textMatrix <- t(textMatrix)
    yLabels <- colnames(trait)
    xLabels <- names(MEs)
    xSymbols <- names(MEs)
    ySymbols <- NULL
    xColorLabels <- TRUE
    par(mar = c(5, 5, 1, 1))
  } else {
    par(mar = c(6, 8.5, 3, 3))
    yLabels <- names(MEs)
    xLabels <- colnames(trait)
    xColorLabels <- FALSE
    xSymbols <- NULL
    ySymbols <- names(MEs)
  }
  cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(10, palette)))(100)
  hm_plot <- WGCNA::labeledHeatmap(
    Matrix = modtraitcor, yLabels = yLabels, xLabels = xLabels,
    ySymbols = ySymbols, xSymbols = xSymbols, colorLabels = FALSE,
    colors = cols, textMatrix = textMatrix, setStdMargins = FALSE,
    cex.text = cex.text, cex.lab.x = cex.lab.x, cex.lab.y = cex.lab.y,
    zlim = c(-1,1), cex.main = 1,
    main = paste("Module-trait relationships")
  )
  return(combined_long)
}

# MEtrait_new <- module_trait_cor_new(
#   exp = exp_df,
#   metadata = metdata_df,
#   MEs = net$MEs,
#   cor_method = "spearman",
#   continuous_trait = FALSE,
#   transpose = FALSE
# )
idx <- 2:ncol(coldata_merged)

# Create an empty Dataframe
df = data.frame(matrix(nrow = 66, ncol = 4))

cor_list <- lapply(idx, function(x) {
  MEtrait_result <- module_trait_cor_new(exp = exp_df,
                              MEs = net$MEs,
                              metadata = coldata_merged[, x, drop = FALSE],
                              cor_method = "pearson")
  MEtrait_result
  return(MEtrait_result)
})


