
#installing and loading required packages----
pacman::p_load(char = c('lubridate',"gtsummary", 'tidyverse', "dplyr", "here", "rio", "scales",
"boot","magrittr",  "mvtnorm", "zoo", "patchwork", "mgcv", "PropCIs", "writexl","ggsignif",
"ggpubr", "ggeasy", "cowplot","ggExtra", "PupillometryR","hrbrthemes", "ggstance","survival",
"survminer","sysfonts","showtext","nlme"))

install.packages("remotes")
remotes::install_github("aertslab/SCENIC")
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
library(nichenetr)
library(Seurat)
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(SCopeLoomR)


dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")


dbFiles <- download.file(url = "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather", destfile = "scenic_database",mode = "wb")

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
all_merged_subset_labelled$cell_cluster <- paste0(all_merged_subset_labelled@active.ident) 
DimPlot(all_merged_subset_labelled,reduction = 'umap.harmony',label = T,repel = T)+NoLegend()
# subset for T cells ----
T_cells <- subset(all_merged_subset_labelled,
                  idents=c("CD3+ T cells"),
                  invert=F)
DimPlot(T_cells,
        reduction = 'umap.harmony',
        label = T,
        repel = T,
        split.by = "HIV_Status")+
  NoLegend()

# Perform differential expression
de_genes3MvsNeg <- FindMarkers(T_cells, 
                               group.by = "HIV_Status",
                               ident.1 = "HIV+ ART<3 Months",
                               ident.2 = "HIV-",
                               min.pct = 0.25)

de_genes1YvsNeg <- FindMarkers(T_cells, 
                              group.by = "HIV_Status",
                              ident.1 = "HIV+ ART>1 Year",
                              ident.2 = "HIV-",
                              min.pct = 0.25)
# Extract expression Matrix ----

sce <- as.SingleCellExperiment(all_merged_subset_labelled)
experMat <-counts(sce)
cellInfo <- colData(sce)
cellInfo <- data.frame(seuratClusters=Idents(all_merged_subset_labelled))

# Saving into loom
dir.create("scenic_database/loom")
loom <- build_loom("scenic_database/loom/all_cells_seurat.loom",dgem=experMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

# Input
# Expression matrix
loom <- open_loom("scenic_database/loom/all_cells_seurat.loom")

experMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)
dim(experMat)

# Cell info / phenodata
cellInfo$nGene <- colSums(experMat>0)
head(cellInfo)

cellInfo <- data.frame(cellInfo)
cbind(table(cellInfo$seuratClusters))

dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")


# Initialise scenic settings
library(feather)
hg19_500bp <- read_feather(file.path(dbDir, "hg19-500bp-upstream-7species.mc9nr.feather"))
hg19_tss <- read_feather(file.path(dbDir, "hg19-tss-centered-10kb-7species.mc9nr.feather"))
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
org <- "hgnc"
dbDir <- "scenic_database"
myDatasetTitle <- "all_cells_seurat"
data("defaultDbNames")
dbs <- defaultDbNames[[org]]
dbs <- c("hg19-500bp-upstream-7species.mc9nr.feather", 
         "hg19-tss-centered-10kb-7species.mc9nr.feather")
scenicOptions <- initializeScenic(org = org,
                                  dbDir = dbDir,
                                  dbs = dbs,
                                  datasetTitle = myDatasetTitle,
                                  nCores = 4)
list.files(dbDir)






