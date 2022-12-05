# Author: KeshavPrasadGubbi

Package_List <- c(
  "dplyr", "DESeq2", "pheatmap", "PoiClaClu", "RColorBrewer", "vsn", "EnhancedVolcano", "gplots",
  "org.Mm.eg.db", "stringr", "genefilter", "tidyverse", "AnnotationDbi", "ComplexHeatmap", "DOSE",
  "clusterProfiler", "ggrepel", "GO.db", "GOstats", "gage", "gageData", "GOSemSim", "enrichplot",
  "ggnewscale", "glue", "ggupset", "FactoMineR", "factoextra", "here", "tibble", "edgeR"
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

# # *****************Now to the comparisons*************
### Reduce larger Matrix to smaller one - based on comparison
paste0(Comparison )


# **********************FUNCTIONS******************************************************************
# Function to save generic plots
saveplot <- function(plot, plotname) {
  # Function to save the plots
  extension <- ".jpeg"
  ggsave(filename = file.path(Comparison_path, paste(plotname, glue("_{Comparison}_"), extension),sep = ""),
         plot = plot, dpi = 300, width = 10, height = 10, units = "in")
  #dev.off()
  while (!is.null(dev.list()))  dev.off()
}
# **********************DESeq Analysis********************************
# Sanity Check for DDS
all(rownames(coldata) %in% colnames(countsmatrix))
ncol(countsmatrix) == nrow(coldata)
dim(countsmatrix)

# Create DEseq Object based on design that was chosen through the Comparison_Condition that was chosen at the begininning of the script run.
if (Comparison_Condition == "condition") {
  ## Creating the DESeq Data set Object
  dds <- DESeqDataSetFromMatrix(countData = countsmatrix, colData = coldata, design = ~condition)
}else {
  # dds <- DESeqDataSetFromMatrix(countData = countsmatrix, colData = coldata, design = ~Epithelial_response)
  print("Design Condition Not Found in coldata column!")
}

# Further filtering of low count genes
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep, ]
nrow(dds)
## Applying VST transformation
vsd <- vst(dds, blind = FALSE) # Apply Variance Stabilizing Transformation
vsd_coldata <- colData(vsd) # Creating a SummarizedExperiment objects
dds <- estimateSizeFactors(dds)

##############################for 2D Analysis#############################################
vsd <- varianceStabilizingTransformation(dds)
dds <- estimateDispersions(dds)
wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
summary(rv_wpn)
# q95_wpn <- quantile(rowVars(wpn_vsd), 0.95)
# normalized_input <- wpn_vsd[rv_wpn>q95_wpn,]
# head(normalized_input)
write.csv(rv_wpn, file.path(here(), "mappedcounts.csv"))
###########################################################################

### Euclidean Distance between samples
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Sample_Name
colnames(sampleDistMatrix) <- vsd$Sample_Name
colors <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(255)
(EuclideanDistanceHeatmap <- pheatmap(sampleDistMatrix,
                                      clustering_distance_rows = sampleDists,
                                      clustering_distance_cols = sampleDists,
                                      main = glue("Euclidean Distance of Samples: {Comparison}"),
                                      col = colors))
### Poisson Distance between Samples
poisd <- PoissonDistance(t(counts(dds))) # raw counts or unnormalised data
samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- dds$Sample_Name
colnames(samplePoisDistMatrix) <- dds$Sample_Name
colors <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(255)
(poisson_dist_plot <- pheatmap(samplePoisDistMatrix,
                               clustering_distance_rows = poisd$dd,
                               clustering_distance_cols = poisd$dd,
                               main = glue("Poisson Distance of Samples: {Comparison}"),
                               col = colors
))

#  **************************PCA Plot**********************************
# Functions for Plot aethetics and saving PCA Plots
color_values <- c("red", "red", "red", "red",  "black", "black", "red", "red", "red", "red", "red",  "red",
                  "red", "red", "red", "blue")
# The basic set of common aesthetic settings for PCA plots,
theme.my.own <- list(theme_bw(),
                     geom_point(size = 3),
                     coord_fixed(),
                     scale_y_continuous(breaks = seq(-100, 100, 10),
                                        sec.axis = sec_axis(~ . * 1, labels = NULL, breaks = NULL)),
                     scale_x_continuous(breaks = seq(-100, 100, 10),
                                        sec.axis = sec_axis(~ . * 1, labels = NULL, breaks = NULL)),
                     theme_classic(),
                     geom_hline(yintercept = 0, color = "gray", linewidth = 1),
                     geom_vline(xintercept = 0, color = "gray", linewidth = 1),
                     theme(text = element_text(size = 24),
                           axis.text = element_text(size = 24),
                           #legend.position = "bottom",
                           legend.position = "none", # if one wants to remove the legend
                           aspect.ratio = 1),
                     #geom_text(size = 8),
                     geom_text_repel(size = 8, min.segment.length = 0.1)
)
# PCA Plot Calculation
# Calculating all PCA Values
plotPCA_local <- function(object, intgroup = "condition", ntop = 500, returnData = TRUE, nPC = 4) {
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  ntop <- 500
  # select the ntop genes by variance
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select, ]))
  # summary(pca)
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  } else {
    colData(object)[[intgroup]]
  }
  # assembly the data for the plot
  d <- cbind(
    pca$x[, seq_len(min(nPC, ncol(pca$x))), drop = FALSE],
    data.frame(group = group, intgroup.df, name = colnames(object))
  )
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:nPC]
    # l <- list(pca,d)
    # return(l)
    return(d)
  }
}
## PCA Plot with VST Data
### Function for calculating percentvar for different variables
percentvar_calculation <- function(pcaData_variable) {
  percentvar_variable <- round(100 * attr(pcaData_variable, "percentVar"), digits = 1)
  return(percentvar_variable)
}

pcaData <- plotPCA_local(vsd, intgroup = c("condition", "Sample_Name"), returnData = T)
#pcaData <- plotPCA_local(vsd, intgroup = c("MouseType", "Sample_Name"), returnData = T)
pcaData
percentVar <- percentvar_calculation(pcaData)
percentVar

# PC Plot: PC1 vs PC2
(PCAplot_vst <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Sample_Name, label = Sample_Name)) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle(glue("PCA: {Comparison}")) +
    scale_colour_manual(values = color_values) +
    theme.my.own)
saveplot(PCAplot_vst, plotname = "PCA_PC1vsPC2")

# PCA Plot : PC3 vs PC4
(PCAplot_pc34 <- ggplot(
  pcaData,
  aes(x = PC3,y = PC4, color = Sample_Name, label = Sample_Name)) +
    xlab(paste0("PC3: ", percentVar[3], "% variance")) +
    ylab(paste0("PC4: ", percentVar[4], "% variance")) +
    ggtitle(glue("PCA: {Comparison}")) +
    scale_colour_manual(values = color_values) +
    theme.my.own)
saveplot(PCAplot_pc34, plotname = "PCA_PC3vsPC4")

# ************************FactoExtra************************
# calculate the variance for top 500 gene

rv <- rowVars(assay(vsd))
ntop <- 500
# select the ntop genes by variance
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
res.pca <- PCA(t(assay(vsd)[select, ]), graph = FALSE, scale.unit = FALSE)
summary.PCA(res.pca)

# Visualize eigenvalues/variances
fviz_screeplot(res.pca, addlabels = TRUE)
eig.val <- get_eigenvalue(res.pca)
var <- get_pca_var(res.pca)

## Genes + PCA Biplots
heat.colors <- brewer.pal(6, "RdYlBu")
## Genes + PCA Biplots
(Genes_Biplot <- fviz_pca_biplot(res.pca, repel = TRUE, labelsize = 6))
(Genes_Biplot <- Genes_Biplot + theme(text = element_text(size = 17),
                                      axis.title = element_text(size = 17),
                                      axis.text = element_text(size = 17),
                                      legend.key.size = unit(3, 'cm'),
                                      legend.text = element_text(size = 20))
)
saveplot(Genes_Biplot, plotname = "Genes_Biplot")

(Genes_contributions_Biplot <- fviz_pca_var(res.pca, col.var = "contrib", repel = TRUE,
                                            labelsize = 6,
                                            gradient.cols = c("Gray", "blue", "pink", "yellow", "orange",
                                                              "green", "red", "black")
                                            ))
(Genes_contributions_Biplot <- (Genes_contributions_Biplot + theme(text = element_text(size = 17),
                                                       axis.title = element_text(size = 17),
                                                       axis.text = element_text(size = 17),
                                                       legend.key.size = unit(1, 'cm'),
                                                       legend.text = element_text(size = 20)
                                                        )))
saveplot(Genes_contributions_Biplot, plotname = "Genes_contributions_Biplot")
# Contributions of variables to PC2
(top25_genes_dim2 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 25))
saveplot(top25_genes_dim2, plotname = "top25_genes_dim2")
# # Contributions of variables to PC1
(top25_genes_dim1 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 25))
saveplot(top25_genes_dim1, plotname = "top25_genes_dim1")

# ********************************DGE Results********************************
### Running the differential expression pipeline
dds <- DESeq(dds)

### Building the results table
### 2nd term will be the Nr.(Infected)
switch(Comparison,
"BL6_InfectedVsControl" = {res <- results(dds, cooksCutoff = FALSE, independentFiltering = FALSE,
                                          contrast = c("condition", "Infected", "control"))},#BL6_InfectedVsControl
# "adult_GFVsd7_GF"   = {res <- results(dds, contrast = c("condition", "d7", "adult"))}, #adult_GFVsd7_GF
# "adult_WTVsd7_WT"   = {res <- results(dds, contrast = c("condition", "d7", "adult"))}, #adult_WTVsd7_WT
# "d7_WTVsd7_spF"     = {res <- results(dds, contrast = c("MouseType", "SPF", "WLD"))}, # d7_WTVsd7_spF
# "d7_GFVsd7_SPF"     = {res <- results(dds, contrast = c("MouseType", "SPF", "GF"))},# d7_GFVsd7_SPF
# "adult_GFVsadult_WT" = {res <- results(dds, contrast = c("MouseType", "WLD", "GF"))},# adult_GFVsadult_WT
# "adult_GFVsadult_SPF" = {res <- results(dds, contrast = c("MouseType", "SPF", "GF"))},# adult_GFVsadult_SPF
)
write.csv(as.data.frame(res), file = file.path(Comparison_path , glue("DGE_Results_{Comparison}.csv")))

### Histogram of p-values
hist(res$pvalue, breaks = 100, col = "grey50", border = "blue")

# Map Gene symbols to Ensembl and Entrez ID
resdf <- tibble::rownames_to_column(as.data.frame(res), var = "symbol")
gn <- resdf$symbol
# Mapping the Symbol to ENTREZ ID
entrez <- mapIds(org.Mm.eg.db, keys = gn, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
ensembl_id <- mapIds(org.Mm.eg.db, keys = gn, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
resdf$entrez <- entrez
resdf$ensemblID <- ensembl_id
resdf <- resdf %>% filter(!is.na(symbol) & !is.na(entrez))

## Differentially Expressed Genes that are statistically Significant
resdf$diffexpressed <- "NS"
# if log2Foldchange > 1.0 and pvalue < 0.05, set as "UP"
resdf$diffexpressed[resdf$log2FoldChange > 1.0 & resdf$pvalue < 0.05] <- "UP"
# if log2Foldchange < -1.0 and pvalue < 0.05, set as "DOWN"
resdf$diffexpressed[resdf$log2FoldChange < -1.0 & resdf$pvalue < 0.05] <- "DOWN"
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
# resdf$delabel <- NA
# resdf$delabel[resdf$diffexpressed != "NS"] <- resdf$symbol[resdf$diffexpressed != "NS"]

# Significant DE Genes Table
significantgenes_df <- resdf[(abs(resdf$log2FoldChange) > 1) & (resdf$pvalue < 0.05), ]
write.csv(significantgenes_df, file = file.path(Comparison_path , glue("SignificantDEgenes_{Comparison}.csv")))

## ********************************Volcano Plots based on Enhanced Volcano********************************
# Volcano Plot with pvalue
volcano1 <- EnhancedVolcano(resdf,
                            lab = resdf$symbol,
                            x = "log2FoldChange",
                            y = "pvalue",
                            xlab = bquote(~ Log[2] ~ "FoldChange"),
                            pCutoff = 0.05,
                            FCcutoff = 1.0,
                            title = glue("DE genes: Log2FoldChange Vs -Log10 Pvalue: {Comparison}"),
                            subtitle = bquote(~ Log[2] ~ "|FoldChange| = 1, pvalue < 0.05"),
                            pointSize = 2.0,
                            labSize = 6,
                            boxedLabels = FALSE,
                            gridlines.major = FALSE,
                            gridlines.minor = FALSE,
                            colAlpha = 0.5,
                            legendPosition = "bottom",
                            legendLabSize = 12,
                            legendIconSize = 4.0,
                            drawConnectors = T,
                            widthConnectors = 0.75,
                            max.overlaps = 12,
                            axisLabSize = 22
)
(volcano1 <- volcano1 + scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 1),
                                           sec.axis = sec_axis(~ . * 1, labels = NULL, breaks = NULL)) +
    scale_x_continuous(limits = c(-5, 10), breaks = seq(-5, 10, 2),
                       sec.axis = sec_axis(~ . * 1, labels = NULL, breaks = NULL)) )
#xlab(expression(DownRegulated %<->% UpRegulated)))
saveplot(volcano1, plotname = "Volcano_pvalue")

# Volcano Plot with padj.
volcano2 <- EnhancedVolcano(resdf,
                            lab = resdf$symbol,
                            x = "log2FoldChange",
                            y = "padj",
                            xlab = bquote(~ Log[2] ~ "FoldChange"),
                            ylab = bquote(~ Log[10] ~ "Padj"),
                            pCutoff = 0.05,
                            FCcutoff = 1.0,
                            title = glue("DE genes: Log2FoldChange Vs -Log10 Padj: {Comparison}"),
                            subtitle = bquote(~ Log[2] ~ "|FoldChange| = 1, pvalue < 0.05"),
                            pointSize = 2.0,
                            labSize = 6.0,
                            boxedLabels = FALSE,
                            gridlines.major = FALSE,
                            gridlines.minor = FALSE,
                            colAlpha = 0.5,
                            legendPosition = "bottom",
                            legendLabSize = 12,
                            legendIconSize = 4.0,
                            drawConnectors = T,
                            widthConnectors = 0.75,
                            max.overlaps = 12,
                            axisLabSize = 22
)
(volcano2 <- volcano2 + scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1),
                                           sec.axis = sec_axis(~ . * 1, labels = NULL, breaks = NULL)) +
    scale_x_continuous(limits = c(-5, 10), breaks = seq(-5, 10, 2),
                       sec.axis = sec_axis(~ . * 1, labels = NULL, breaks = NULL)))
#xlab(expression(DownRegulated %<->% UpRegulated))
saveplot(volcano2, plotname = "Volcano_padj")

### Number of Genes from different strains that are contributing to UP/DOWN regulation.
significantgenes_df_UP <- significantgenes_df[(significantgenes_df$log2FoldChange) > 1, ]# UP Regulation Table
significantgenes_df_DOWN <- significantgenes_df[(significantgenes_df$log2FoldChange) < -1, ]# DOWN Regulation Table
nrow(significantgenes_df_UP)
nrow(significantgenes_df_DOWN)
write.csv(significantgenes_df_UP, file.path(Comparison_path, glue("SignificantDE_UPgenes_{Comparison}.csv")))
write.csv(significantgenes_df_DOWN, file.path(Comparison_path, glue("SignificantDE_DOWNgenes_{Comparison}.csv")))

## ********************************Z-score based Gene Heatmaps********************************
# Determining the significant Genes based on Log2FC and pvalue thresholds
sigs2df <- resdf[(abs(resdf$log2FoldChange) > 1) & (resdf$pvalue < 0.05), ]

# mat <- counts(dds1, normalized = TRUE)[rownames(sigsdf),]
mat <- counts(dds, normalized = TRUE)[(significantgenes_df$symbol) %in% rownames(counts(dds)), ]
mat.zs <- t(apply(mat, 1, scale)) # Calculating the zscore for each row
colnames(mat.zs) <- coldata$Sample_Name# need to provide correct sample names for each of the columns

(AllGenes_Heatmap <- Heatmap(mat.zs,
                            cluster_columns = TRUE,
                            cluster_rows = TRUE,
                            column_labels = colnames(mat.zs),
                            name = glue("DE Genes- {Comparison}"),
                            show_row_names = FALSE,
                            use_raster = TRUE,
                            raster_quality = 10,
                            column_names_gp = grid::gpar(fontsize = 22),
                            #row_labels = sigs2df[rownames(mat2.zs), ]$symbol
                            heatmap_legend_param = list(legend_direction = "horizontal",
                                                        legend_width = unit(x= 5, units = "cm"))))

jpeg(file = file.path(Comparison_path, glue("/DEGenes_heatmap1_{Comparison}.jpeg") ),
    width = 1000, height = 1000, units = "px", pointsize = 12,
    bg = "white", res = NA, family = "", restoreConsole = TRUE,
    type = "windows",
    symbolfamily = "default"
)
draw(AllGenes_Heatmap, heatmap_legend_side = "bottom")
while (!is.null(dev.list()))  dev.off()
# Long heatmap
LongHeatMap_Allgenes <- Heatmap(mat.zs,
                                cluster_columns = TRUE,
                                cluster_rows = TRUE,
                                column_labels = colnames(mat.zs),
                                row_labels = rownames(mat.zs),#sigsdf[rownames(mat2.zs), ]$symbol
                                name = "All Genes",
                                show_row_names = TRUE,
                                use_raster = TRUE,
                                raster_quality = 5,
                                column_names_gp = grid::gpar(fontsize = 18),
                                #row_labels = sigsdf[rownames(mat2.zs), ]$symbol
)
jpeg(file = file.path(Comparison_path, glue("/DEGenes_heatmap2_{Comparison}.jpeg") ),
     width = 1000, height = 2000, units = "px", pointsize = 12, bg = "white", res = NA,
     family = "", restoreConsole = TRUE, type = "windows", symbolfamily = "default")
draw(LongHeatMap_Allgenes, heatmap_legend_side = "bottom")
while (!is.null(dev.list()))  dev.off() #dev.off()


### Heatmap with tighter constraints (all genes together!)
sigs1df <- resdf[(resdf$baseMean > 250) & (abs(resdf$log2FoldChange) > 2) & (resdf$pvalue < 0.05), ]
mat1 <- counts(dds, normalized = TRUE)[(sigs1df$symbol), ]
mat1.zs <- t(apply(mat1, MARGIN = 1, scale)) # Calculating the zscore for each row
colnames(mat1.zs) <- coldata$Sample_Name # need to provide correct sample names for each of the columns

Tightconstraints_Heatmap <- Heatmap(mat1.zs,
                                    cluster_columns = TRUE,
                                    cluster_rows = TRUE,
                                    column_labels = colnames(mat1.zs),
                                    name = glue("DE Genes - {Comparison}"),
                                    row_labels = rownames(mat1.zs),
                                    column_names_gp = grid::gpar(fontsize = 20),
                                    row_names_gp = grid::gpar(fontsize = 20),
                                    heatmap_legend_param = list(legend_direction = "horizontal",
                                                                legend_width = unit(x = 5, units = "cm")))
jpeg(file = file.path(Comparison_path, glue("/DEGenes_heatmap3_{Comparison}.jpeg") ),
     width = 1000, height = 1500, units = "px", pointsize = 12, bg = "white", res = NA,
     family = "", restoreConsole = TRUE, type = "windows", symbolfamily = "default")
draw(Tightconstraints_Heatmap, heatmap_legend_side = "bottom")
while (!is.null(dev.list()))  dev.off() #dev.off()

# ********************************Functional Analysis using Cluster Profiler********************************
## GO over-representation analysis
### GO Terms for UP Regulated Genes
UPgene_ENS_ID <- (significantgenes_df_UP$ensemblID)
GO_UPRegResults <- enrichGO(gene = UPgene_ENS_ID, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL",
                            ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                            readable = TRUE)
GO_UPRegResults_df <- as.data.frame(GO_UPRegResults)
write.csv(GO_UPRegResults_df, file.path(Comparison_path , glue("GO_UPRegResults_df_{Comparison}.csv")))

# Functional Analysis Plots
if (nrow(GO_UPRegResults_df > 0)) {
  GO_UPReg_Barplot <- plot(barplot(GO_UPRegResults, showCategory = 25, font.size = 15,
                                   title = "UpRegulated", label_format = 45))
  saveplot(plot = GO_UPReg_Barplot, plotname = "GO_UPReg_Barplot")
  GO_UPReg_Dotplot <- plot(dotplot(GO_UPRegResults, showCategory = 25, font.size = 15,
                                   title = "UpRegulated", label_format = 45))
  saveplot(plot = GO_UPReg_Dotplot, plotname = "GO_UPReg_Dotplot")
  GO_UPReg_Cnetplot <- plot(cnetplot(GO_UPRegResults, showCategory = 15, font.size = 20))
  saveplot(plot = GO_UPReg_Cnetplot, plotname = "GO_UPReg_Cnetplot")
  GO_UPReg_Heatplot <- plot(heatplot(GO_UPRegResults, foldChange = 1))
  saveplot(plot = GO_UPReg_Heatplot, plotname = "GO_UPReg_Heatplot")
  edox2 <- pairwise_termsim(GO_UPRegResults)
  GO_UPReg_enrichtreeplot <- plot(treeplot(edox2))
  saveplot(plot = GO_UPReg_enrichtreeplot, plotname = "GO_UPReg_enrichtreeplot")
  (GO_UPReg_emapplot <- emapplot(edox2, showCategory = 25, repel = TRUE))
  saveplot(plot = GO_UPReg_emapplot, plotname = "GO_UPReg_emapplot")
} else {
  print("There were 0 rows in ORA analysis results!")
}

### GO Terms for DOWN Regulated Genes
DOWNgene_ENS_ID <- (significantgenes_df_DOWN$ensemblID)
GO_DOWNRegResults <- enrichGO(gene = DOWNgene_ENS_ID, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "BP",
                              pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
GO_DOWNRegResults_df <- as.data.frame(GO_DOWNRegResults)
write.csv(GO_DOWNRegResults_df, file.path(Comparison_path , glue("GO_DOWNRegResults_df_{Comparison}.csv")))

if (nrow(GO_DOWNRegResults_df > 0)) {
  GO_DOWNReg_Barplot <- plot(barplot(GO_DOWNRegResults, showCategory = 25, font.size = 15,
                                     title = "DOWNRegulated", label_format = 45))
  saveplot(plot = GO_DOWNReg_Barplot, plotname = "GO_DOWNReg_Barplot")
  GO_DOWNReg_Dotplot <- plot(dotplot(GO_DOWNRegResults, showCategory = 25, font.size = 15,
                                     title = "DOWNRegulated", label_format = 45))
  saveplot(plot = GO_DOWNReg_Dotplot, plotname = "GO_DOWNReg_Dotplot")
  GO_DOWNReg_Cnetplot <- plot(cnetplot(GO_DOWNRegResults, showCategory = 15, font.size = 20))
  saveplot(plot = GO_DOWNReg_Cnetplot, plotname = "GO_DOWNReg_Cnetplot")
  GO_DOWNReg_Heatplot <- plot(heatplot(GO_DOWNRegResults, foldChange = 1))
  saveplot(plot = GO_DOWNReg_Heatplot, plotname = "GO_DOWNReg_Heatplot")
  edox2 <- pairwise_termsim(GO_DOWNRegResults)
  GO_DOWNReg_enrichtreeplot <- plot(treeplot(edox2))
  saveplot(plot = GO_DOWNReg_enrichtreeplot, plotname = "GO_DOWNReg_enrichtreeplot")

  GO_DOWNReg_emapplot <- emapplot(edox2, showCategory = 25, repel = TRUE)
  saveplot(plot = GO_DOWNReg_emapplot, plotname = "GO_DOWNReg_emapplot")

}else {
  print("There were 0 rows in ORA analysis results!")
}


# Log2Fold Change Table

# Creating a list of samplenmes that can be provided
samplelist <- as.list(coldata$Sample_Name)

# Design a generic function
run_individualDESeq <- function(namelist, cd, cm) {
  # cd is coldata, cm is countsmatrix
  # Create an empty list for the results
  resultslist <- list()
  for (i in 1:14) {
    name <- namelist[[i]]
    print(glue("Working with Strain:{name}"))
    cm_input <- cm[, c(i, 15:16)]
    coldata_input <- cd[c(i, 15:16), ]
    rownames(coldata_input)
    colnames(cm_input)
    all(rownames(coldata_input) %in% colnames(cm_input))
    ncol(cm_input) == nrow(coldata_input)
    dds <- DESeqDataSetFromMatrix(countData = cm_input, colData = coldata_input, design = ~condition)
    dds <- DESeq(dds)
    dds1 <- DESeq(dds, minReplicatesForReplace = Inf)
    res <- results(dds1, cooksCutoff = FALSE, independentFiltering = FALSE,
                   contrast = c("condition", "Infected", "control"))
    res_df <- as.data.frame(res) # convert the results table to a df
    res_df <- tibble::rownames_to_column(res_df, "symbol")
    res_df <- res_df %>% filter(!is.na(log2FoldChange))
    colnames(res_df)[which(names(res_df) == "log2FoldChange")] <- glue("log2FoldChange_{name}")
    res_df <- as.data.frame(res_df)
    print(glue("completed {name} strain DESeq calculations and produced result tables!!"))
    # Appending each new results table into a resultslist
    resultslist[[length(resultslist) + 1]] <- res_df
  }
  return(resultslist)
}

# Calling the function to run DESeq for each Strain
resulttable <- run_individualDESeq(namelist = samplelist, cd = coldata, cm = countsmatrix)
# Mapping -> converting the list of objects to list of dataframes
list_df <- Map(as.data.frame, resulttable)
# merge all data frames in list based on symbol column
table <- list_df %>% reduce(full_join, by = "symbol")

l2fctable <- table %>% select("symbol",
  "log2FoldChange_T", "log2FoldChange_S54", "log2FoldChange_S55", "log2FoldChange_L57",
  "log2FoldChange_L57","log2FoldChange_A58","log2FoldChange_L60","log2FoldChange_S61",
  "log2FoldChange_A62","log2FoldChange_L63","log2FoldChange_A64","log2FoldChange_S65",
  "log2FoldChange_L66","log2FoldChange_A68","log2FoldChange_L69")

write.csv(l2fctable, file = file.path(Comparison_path , glue("LogFCTable_{Comparison}.csv")))

# Calculating CPM Values

run_individualDESeq2 <- function(namelist, cd, cm) {
  # cd is coldata, cm is countsmatrix
  # Create an empty list for the results
  cpmlist <- list()
  for (i in 1:14) {
    name <- namelist[[i]]
    print(glue("Working with Strain:{name}"))
    cm_input <- cm[, i]
    coldata_input <- cd[i, ]
    # print(head(cm_input, 3))
    # print(head(coldata_input, 3))
    # CPM Values
    # as DGEList
    dge_er <- DGEList(counts = cm_input)
    dim(dge_er)
    colnames(dge_er)
    # dge_er$samples
    ## calculate norm. factors
    nr <- calcNormFactors(dge_er)
    ## get normalized counts
    cpm_df <- as.data.frame(cpm(nr))
    colnames(cpm_df) <- glue("CPM_{name}")
    cpm_df <- tibble::rownames_to_column(cpm_df, "symbol")
    # Appending each new result table into cpmlist
    cpmlist[[length(cpmlist) + 1]] <- cpm_df
  }
  return(cpmlist)
}

cpmT <- run_individualDESeq2(namelist = samplelist, cd = coldata, cm = countsmatrix)
# Mapping -> converting the list of objects to list of dataframes
list_df2 <- Map(as.data.frame, cpmT)
# merge all data frames in list based on symbol column
cpmtable <- list_df2 %>% reduce(full_join, by = "symbol")
nrow(cpmtable)
rownames(cpmtable) <- cpmtable[, 1]
cpmtable <- subset(cpmtable, select = -symbol)
keep <- rowSums(cpmtable) > 0
cpmtable <- cpmtable[keep, ]
nrow(cpmtable)

write.csv(cpmtable, file = file.path(Comparison_path , glue("CPMTable_{Comparison}.csv")))

sessionInfo()
