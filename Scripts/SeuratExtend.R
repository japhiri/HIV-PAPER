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
load("Data/Single_Cell_Data/all_merged_subset_labelled.RData")
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
DimPlot(all_merged_subset_labelled,reduction = 'umap.harmony',label = TRUE,label.size = 6)+NoLegend()+
  labs(x='UMAP-1',y='UMAP-2')
DimPlot2(all_merged_subset_labelled,label = TRUE, box = TRUE, label.color = 'black',repel = TRUE, theme = NoLegend())
FeaturePlot3(all_merged_subset_labelled,
             reduction = "umap.harmony",
             color = "ryb", 
             feature.1 = "CD3D",
             feature.2 = "CD19",
             feature.3 = "CSF3R")
feature_plot <- FeaturePlot3.grid(all_merged_subset_labelled,
                  reduction = "umap.harmony",
                  features = c("CD3D","CD19","CSF3R","FCGR3A","CD14","CD79A","CD3E","CD8A","CD4","NKG7","PRF1","GZMB"))+
  labs(x="UMAP-1",y="UMAP-2")
feature_plot
ggsave(here('HIV PAPER FIGURES/Figure 1','feature_plot.png'),
       plot=((feature_plot|plot_layout(ncol = 1,nrow = 2,widths = c(1,1,1,1)))),
       width = 8,height = 8,unit='in',dpi=300)
# Headmap for differential expressed gened for each cluster
all_merged_subset_labelled <- FindVariableFeatures(all_merged_subset_labelled,
                     selection.method = 'vst',
                     verbose = TRUE)
genes <- VariableFeatures(all_merged_subset_labelled)
toplot <- CalcStats(all_merged_subset_labelled,
                    features = genes,
                    method = 'zscore',
                    order = 'p',
                    n = 5)
cluster_heatmap <- Heatmap(toplot,lab_fill='zscore',plot.margin = margin(l = 30))    
cluster_heatmap <- as.ggplot(cluster_heatmap)
ggsave(here('HIV PAPER FIGURES/Figure 1','Fig1j.png'),
       plot=((cluster_heatmap|plot_layout(ncol = 2,nrow = 2,widths = c(1,1,1,1)))),
       width = 8,height = 10,unit='in',dpi=300)

ggsave(here('HIV PAPER FIGURES/Figure 1','Fig1j.pdf'),
       plot=((cluster_heatmap|plot_layout(ncol = 2,nrow = 2,widths = c(1,1,1,1)))),
       width = 8,height = 10,unit='in',dpi=300)

# generate a waterfall plot
# Downloading the C7 Immunologic Signature gene sets
library(msigdbr)
C7_human <- msigdbr(species = "Homo sapiens", category = "C7")
genesets <- split(C7_human$gene_symbol, C7_human$gs_name)

# Run GeneSetAnalysis with C7 gene sets
all_merged_subset_labelled <- GeneSetAnalysis(all_merged_subset_labelled,
                                              genesets = genesets)


all_merged_subset_labelled <- GeneSetAnalysis(all_merged_subset_labelled,
                                              genesets = hall50$human)


matr <- all_merged_subset_labelled@misc$AUCell$genesets
all_merged_subset_labelled$cell_clusters <- Idents(all_merged_subset_labelled)
Idents(all_merged_subset_labelled)<-"HIV_Status"
WaterfallPlot(matr,
              f=all_merged_subset_labelled$HIV_Status,
              ident.1 = "HIV+ ART>1 Year",
              ident.2 = "HIV-")
GSEAplot(all_merged_subset_labelled,
         ident.1 = "Neutrophils",
         geneset = hall50$human$HALLMARK_TNFA_SIGNALING_VIA_NFKB)
# Barplot
ClusterDistrBar(origin = all_merged_subset_labelled$sample,
                cluster = all_merged_subset_labelled$cell_clusters,
                flip = TRUE, percent = TRUE)
# Create an enhanced violin plot
genes <- c("CD3D","CD79A","CSF3R")
cells <- colnames(all_merged_subset_labelled)[all_merged_subset_labelled$cell_clusters %in% c("B cells","Neutrophils","CD3+ T cells")]
VlnPlot2(all_merged_subset_labelled,
         features = genes,
         ncol = 1, pt=FALSE,
         #group.by = "cell_clusters",
         #split.by = "HIV_Status",
         #stat.method = "wilcox.test",
         hide.ns = TRUE,
         cells = cells)

# Gene set enrichment analysis using seuratextend
library(SeuratExtend)
library(dplyr)
options(max.print = 12,spe="human")

all_merged_subset_labelled <- GeneSetAnalysisGO(all_merged_subset_labelled,
                                                #parent = "immune_system_process",
                                                nCores = 16)
matr <- all_merged_subset_labelled@misc$AUCell$immune_system_process
matr <- RenameGO(matr)
head(matr, 2,3)

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
# Subsetting for neutrophils----
Neutrophils <- subset(all_merged_subset_labelled,
                      idents = "Neutrophils",
                      invert = F)
Neutrophils$Condition <- paste0(Neutrophils$sample,Neutrophils$HIV_Status,Neutrophils$Carriage_Status)
# Headmap for differential expressed genes neutrophils in each HIV status
Idents(Neutrophils) <- "HIV_Status"
Neutrophils <- FindVariableFeatures(Neutrophils,
                                    selection.method = 'vst',
                                    verbose = TRUE)
genes <- VariableFeatures(Neutrophils)
toplot <- CalcStats(Neutrophils,
                    features = genes,
                    method = 'zscore',
                    order = 'p',n = 10
                    )
cluster_heatmap <- Heatmap(toplot,lab_fill='zscore',plot.margin = margin(l = 30))    
cluster_heatmap <- as.ggplot(cluster_heatmap)
cluster_heatmap
ggsave(here('HIV PAPER FIGURES/Figure 1','NeutHeatMap.png'),
       plot=((cluster_heatmap|plot_layout(ncol = 2,nrow = 2,widths = c(1,1,1,1)))),
       width = 8,height = 10,unit='in',dpi=300)

ggsave(here('HIV PAPER FIGURES/Figure 1','NeutHeatMap.pdf'),
       plot=((cluster_heatmap|plot_layout(ncol = 2,nrow = 2,widths = c(1,1,1,1)))),
       width = 10,height = 10,unit='in',dpi=300)

