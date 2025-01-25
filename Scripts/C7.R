library(Seurat)
library(SeuratExtend)
library(ggplotify)
library(ggplot2)
library(here)
#installing and loading required packages----
pacman::p_load(char = c('lubridate',"gtsummary", 'tidyverse', "dplyr", "here", "rio", "scales", "boot", 
                        "magrittr",  "mvtnorm", "zoo", "patchwork", "mgcv", "PropCIs", "writexl", 
                        "ggsignif", "ggpubr", "ggeasy", "cowplot","ggExtra", "PupillometryR","hrbrthemes", "ggstance",
                        "survival","survminer","sysfonts","showtext","nlme"))

# Load object ----
load("data/all_merged_subset_labelled.RData")
# Renname idents
all_merged_subset_labelled <- RenameIdents(object = all_merged_subset_labelled,
                                           "Secretory cells" = "Secretory cells",
                                           "Goblet cells" = "Goblet cells",
                                           "CD3+ T cells" = "CD3+ T cells",
                                           "Phagocytes" = "Monocytes/Macrophages",
                                           "B cells" = "B cells",
                                           "Developing ciliated cells" = "Ciliated cells",
                                           "Neutrophils" = "Neutrophils",
                                           "Squamous cells" = 'Squamous cells',
                                           "FOXJ1++ Ciliated cells" = "Ciliated cells",
                                           "Stressed cells" = "Stressed cells",
                                           "B cells" = "B cells",
                                           "Club cells" = "Club cells",
                                           "Deuterosomal cells" = "Ciliated cells",
                                           "BEST++ Cilia++ Ciliated cells" = "Ciliated cells",
                                           "Ionocytes" = "Ionocytes",
                                           "Dendritic cells" = "Dendritic cells",
                                           "Neurons" = "Neurons")

all_merged_subset_labelled <- subset(all_merged_subset_labelled,
                                     idents = c("Stressed cells","Neurons","Club cells"),
                                     invert = T)
Neutrophils <- subset(all_merged_subset_labelled,
                      idents = c("CD3+ T cells"),
                      invert = F)
Idents(Neutrophils)<-"HIV_Status"

# generate a waterfall plot
# Downloading the C7 Immunologic Signature gene sets
library(msigdbr)
C7_human <- msigdbr(species = "Homo sapiens", category = "C7")
genesets <- split(C7_human$gene_symbol, C7_human$gs_name)

# Run GeneSetAnalysis with C7 gene sets
Neutrophils <- GeneSetAnalysis(Neutrophils,
                               genesets = genesets)
matr <- Neutrophils@misc$AUCell$genesets
# Calculate row means across cells for each gene set
row_means <- rowMeans(matr)

# Order gene sets by row mean in decreasing order and select the top 20
top20_indices <- order(row_means, decreasing = TRUE)[1:50]
top20 <- matr[top20_indices, ]


WaterfallPlot(top20,
              f=Neutrophils$HIV_Status,
              ident.1 = "HIV+ ART<3 Months",
              ident.2 = "HIV+ ART>1 Year")









