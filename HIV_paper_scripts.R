# Load required packages----
pacman::p_load(char = c("lubridate","gtsummary", "tidyverse", "dplyr", "here", "rio", "scales", "boot", "Matrix","ggpubr",
                        "magrittr",  "mvtnorm", "zoo", "patchwork", "mgcv", "PropCIs", "writexl","DropletUtils","SeuratWrappers",
                        "Seurat","rjson","R2HTML","DT","cowplot","RCurl","glmGamPoi","DESeq2","ggrepel","rPanglaoDB","ComplexHeatmap",
                        "EnhancedVolcano","RColorBrewer","circlize","rmarkdown","biomaRt","biomartr","clusterProfiler","multinichenetr",
                        "AnnotationDbi","org.Hs.eg.db","CEMiTool","enrichplot","pathview","scmap","SingleR","S4Vectors","TMB","muscat",
                        "SingleCellExperiment","apeglm","edgeR","purrr","tibble","png","RColorBrewer","scran","microViz","CommPath",
                        "reshape2","scater","Azimuth","scCATCH","CellChat","SoupX","knitr","DoubletFinder","ggpubr","viridis","GSVA",
                        "ggmin","cluster","foreach","doParallel","BPCells","ggimage","ggbeeswarm","grid","data.table","scriabin",
                        "clusterExperiment","destiny","gam","corrplot","ggthemes","base64enc","Biobase","CATALYST","dittoSeq","viridis",
                        "DelayedArray","DelayedMatrixStats","limma","lme4","batchelor","HDF5Array","terra","ggrastr","nloptr","ggsignif",
                        'lubridate',"gtsummary", 'tidyverse', "dplyr", "here", "rio", "scales", "boot", 
                        "magrittr",  "mvtnorm", "zoo", "patchwork", "mgcv", "PropCIs", "writexl", 
                        "ggsignif", "ggpubr", "ggeasy", "cowplot","ggExtra", "PupillometryR","hrbrthemes", "ggstance",
                        "survival","survminer","sysfonts","showtext","nlme",'glue'))


#installing and loading required packages----
pacman::p_load(char = c('lubridate',"gtsummary", 'tidyverse', "dplyr", "here", "rio", "scales", "boot", 
                        "magrittr",  "mvtnorm", "zoo", "patchwork", "mgcv", "PropCIs", "writexl", 
                        "ggsignif", "ggpubr", "ggeasy", "cowplot","ggExtra", "PupillometryR","hrbrthemes", "ggstance",
                        "survival","survminer","sysfonts","showtext","nlme"))
library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(here)

# Figure 1 ----
# Load all cell cluster
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
# Define a consistent color palette for 17 clusters
num_clusters <- 17
base_palette <- brewer.pal(min(num_clusters,12),'Paired')
additional_colors <- colorRampPalette(brewer.pal(8,'Set3'))(num_clusters-12)
color_palette <- c(base_palette,additional_colors)
names(color_palette) <- clusters
# Plot UMAP with consistent colors
# Fig1d ----
main_umap <- DimPlot(all_merged_subset_labelled,
        reduction = 'umap.harmony',
        label = TRUE,
        repel = TRUE,
        label.size = 5)+
  scale_color_manual(values = color_palette)+
  #NoLegend()+
  labs(x='UMAP-1',
       y='UMAP-2')+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())
main_umap <- as.ggplot(main_umap)
main_umap


# Save Figure 1d
ggsave(here('HIV PAPER FIGURES/Figure 1','Fig1d.png'),
       plot=((main_umap|plot_layout(ncol = 2,nrow = 2,widths = c(1,1,1,1)))),
       width = 20,height = 8,unit='in',dpi=300)
ggsave(here('HIV PAPER FIGURES/Figure 1','Fig1d.pdf'),
       plot=((main_umap|plot_layout(ncol = 2,nrow = 2,widths = c(1,1,1,1)))),
       width = 20,height = 8,unit='in',dpi=700)
# Figure 1e ----
main_umap_sex <- DimPlot(all_merged_subset_labelled,
                     reduction = 'umap.harmony',
                     group.by = "sex",
                     label = FALSE,
                     repel = TRUE,
                     label.size = 5)+
  #scale_color_manual(values = color_palette)+
  #NoLegend()+
  labs(x='UMAP-1',
       y='UMAP-2')+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_blank())
main_umap_sex <- as.ggplot(main_umap_sex)
main_umap_sex


# Save Figure 1e
ggsave(here('HIV PAPER FIGURES/Figure 1','Fig1e.png'),
       plot=((main_umap_sex|plot_layout(ncol = 2,nrow = 2,widths = c(1,1,1,1)))),
       width = 20,height = 8,unit='in',dpi=300)
ggsave(here('HIV PAPER FIGURES/Figure 1','Fig1e.pdf'),
       plot=((main_umap_sex|plot_layout(ncol = 2,nrow = 2,widths = c(1,1,1,1)))),
       width = 20,height = 8,unit='in',dpi=700)
# Figure 1f ----
main_umap_HIV_Status <- DimPlot(all_merged_subset_labelled,
                         reduction = 'umap.harmony',
                         group.by = "HIV_Status",
                         label = FALSE,
                         repel = TRUE,
                         label.size = 5)+
  #scale_color_manual(values = color_palette)+
  #NoLegend()+
  labs(x='UMAP-1',
       y='UMAP-2')+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_blank())
main_umap_HIV_Status <- as.ggplot(main_umap_HIV_Status)
main_umap_HIV_Status


# Save Figure 1f
ggsave(here('HIV PAPER FIGURES/Figure 1','Fig1f.png'),
       plot=((main_umap_HIV_Status|plot_layout(ncol = 2,nrow = 2,widths = c(1,1,1,1)))),
       width = 20,height = 8,unit='in',dpi=300)
ggsave(here('HIV PAPER FIGURES/Figure 1','Fig1f.pdf'),
       plot=((main_umap_HIV_Status|plot_layout(ncol = 2,nrow = 2,widths = c(1,1,1,1)))),
       width = 20,height = 8,unit='in',dpi=700)

# Figure 1g----
# Extract cluster identities
cluster_info <- Idents(all_merged_subset_labelled)
clusters <- unique(cluster_info)
# Calculate proportion for each cluster
cluster_counts <- table(cluster_info)
cluster_df <- as.data.frame(cluster_counts)
colnames(cluster_df) <- c('Cluster','Count')
# Compute percentages
cluster_df$fraction <- cluster_df$Count / sum(cluster_df$Count)

# Compute the cumulative percentages (top of each rectangle)
cluster_df$ymax <- cumsum(cluster_df$fraction)

# Compute the bottom of each rectangle
cluster_df$ymin <- c(0, head(cluster_df$ymax, n=-1))

# Compute label position
cluster_df$labelPosition <- (cluster_df$ymax + cluster_df$ymin) / 2

# Compute a good label
cluster_df$label <- paste0(cluster_df$Cluster, "\n value: ", cluster_df$Count)

# Make the plot
all_data_pie <- ggplot(cluster_df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Cluster
                                       )) +
  geom_rect() +
  geom_point(aes(y = (ymax + ymin) / 2, x = 3.5, fill = Cluster), 
             shape = 21, size = 0, stroke = 0)+
  #geom_label(x=4,
             #aes(y=labelPosition, label=label), 
             #size=3) +
  scale_fill_manual(values = color_palette)+
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void()+
  theme(#legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = 'bold'))+
  guides(fill=guide_legend(override.aes = list(shape=21, size=5,color=NA)))
all_data_pie

# Save Figure 1g
ggsave(here('HIV PAPER FIGURES/Figure 1','Fig1g.png'),
       plot=((all_data_pie|plot_layout(ncol = 3,nrow = 2,widths = c(1,1,1,1),heights = c(1,1,1,1)))),
       width = 20,height = 7,unit='in',dpi=300)
ggsave(here('HIV PAPER FIGURES/Figure 1','Fig1g.pdf'),
       plot=((all_data_pie|plot_layout(ncol = 3,nrow = 2,widths = c(1,1,1,1),heights = c(1,1,1,1)))),
       width = 20,height = 7,unit='in',dpi=700)

# Figure 1h----
# Dotplot for each cluster
#load the csv file with marker_classes
Markers <- read.csv('HIV PAPER FIGURES/Markes.csv')
desired_order <- c("CD3+ T cells","B cells","Neutrophils","Dendritic cells",
                   "Phagocytes","Goblet cells","Secretory cells","Basal cells",
                   "Ciliated cells","Squamous cells","Ionocytes")
all_merged_subset_labelled <- SetIdent(all_merged_subset_labelled,
                                       value = factor(Idents(all_merged_subset_labelled),
                                                      levels = desired_order))
main_dotplot <- DotPlot(all_merged_subset_labelled,
        features = unique(Markers$Marker),
        cols = c("yellow","blue"),
        #split.by = "HIV_Status",
        dot.scale = 9,
        scale = TRUE,
        scale.by = "radius")+
  RotatedAxis()+
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   size = 10),
        axis.text.y = element_text(size = 10,
                                   face = 'bold'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10,
                                    face = 'bold'),
        axis.title = element_blank())
main_dotplot <- as.ggplot(main_dotplot)
main_dotplot

# Save Figure 1h
ggsave(here('HIV PAPER FIGURES/Figure 1','Fig1h.png'),
       plot=((main_dotplot|plot_layout(ncol = 1,nrow = 1,widths = c(1,1,1,1),heights = c(1,1,1,1)))),
       width = 13,height = 4.5,unit='in',dpi=300)
ggsave(here('HIV PAPER FIGURES/Figure 1','Fig1h.pdf'),
       plot=((main_dotplot|plot_layout(ncol = 1,nrow = 1,widths = c(1,1,1,1),heights = c(1,1,1,1)))),
       width = 13,height = 4.5,unit='in',dpi=700)

# Figure 1i ----
all_merged_subset_labelled$cell_clusters <- Idents(all_merged_subset_labelled)
# Calculate cluster proportion 
cluster_df <- as.data.frame(table(all_merged_subset_labelled$cell_clusters,
                                  all_merged_subset_labelled$HIV_Status))
colnames(cluster_df)<-c('Cluster','HIV Status','Count')

# Fig 1i
bargraph_main <- cluster_df %>%
  ggplot(aes(x=`HIV Status`,
             y=Count, fill=Cluster))+
  geom_col(position = "fill")+
  #scale_fill_manual(values = color_palette)+
  scale_fill_manual(values = c(
    "CD3+ T cells" = "#A3D981",
    "B cells" = "#F3C7D4",
    "Neutrophils" = "#C3A1A9",
    "Dendritic cells" = "#FFEB5D",
    "Phagocytes" = "#8BD1D1",
    "Goblet cells" = "#3582C4",
    "Secretory cells" = "#C3DBEF",
    "Basal cells" = "#47A248",
    "Squamous cells" = "#B5603B",
    "Ionocytes" = "#DDA92C"
  )) +
  labs(x='',
       y='Frequency of cells')+
  theme_classic2()+
  theme(legend.position = "right",
        legend.key.size = unit(0.5,'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(angle = 0,hjust = 1,face = 'bold'))+
  coord_flip()
bargraph_main

# Save Figure 1i
ggsave(here('HIV PAPER FIGURES/Figure 1','Fig1i.png'),
       plot=((bargraph_main|plot_layout(ncol = 1,nrow = 1,widths = c(1,1,1,1),heights = c(1,1,1,1)))),
       width = 13,height = 3,unit='in',dpi=300)
ggsave(here('HIV PAPER FIGURES/Figure 1','Fig1i.pdf'),
       plot=((bargraph_main|plot_layout(ncol = 1,nrow = 1,widths = c(1,1,1,1),heights = c(1,1,1,1)))),
       width = 13,height = 3,unit='in',dpi=700)



# Differential gene expression analysis in neutrophils ----
# Perform defferential expression analysis within t he same cell type across conditions
# Perform pseudobulking
all_merged_subset_labelled$condition <- paste(all_merged_subset_labelled$cell_clusters,
                                              all_merged_subset_labelled$HIV_Status,sep = '_')
Idents(all_merged_subset_labelled) <- 'condition'
view(all_merged_subset_labelled@meta.data)
pseudo_exp <- AggregateExpression(all_merged_subset_labelled,
                                  assays = 'RNA',
                                  slot = 'counts',
                                  return.seurat = TRUE,
                                  group.by = c('sample',
                                               'HIV_Status',
                                               'cell_clusters'))
tail(Cells(pseudo_exp))
view(pseudo_exp@meta.data)

pseudo_exp$cell_clusters <- sapply(strsplit(Cells(pseudo_exp),
                                            split = '_'),'[',3)
pseudo_exp$HIV_Status <- sapply(strsplit(Cells(pseudo_exp),
                                     split = '_'),'[',2)
pseudo_exp$sample <- sapply(strsplit(Cells(pseudo_exp),
                                         split = '_'),'[',1)
pseudo_exp$condition <- paste(pseudo_exp$cell_clusters,
                              pseudo_exp$HIV_Status,sep = '-')

# Set the Idents to condition 
Idents(pseudo_exp)<-'condition'
# Neutrophils differential expression
Neut_3MvsNeg <- FindMarkers(object = pseudo_exp,
                            ident.1 = 'Neutrophils-HIV+ ART<3 Months',
                            ident.2 = 'Neutrophils-HIV-',
                            #min.pct = 0.25,
                            test.use = 'DESeq2')
Neut_3MvsNeg <- Neut_3MvsNeg %>%
  filter(p_val<0.05, abs(avg_log2FC)>1.5)
Neut_1YvsNeg <- FindMarkers(object = pseudo_exp,
                            ident.1 = 'Neutrophils-HIV+ ART>1 Year',
                            ident.2 = 'Neutrophils-HIV-',
                            #min.pct = 0.25,
                            test.use = 'DESeq2')
Neut_1YvsNeg <- Neut_1YvsNeg %>%
  filter(p_val<0.05, abs(avg_log2FC)>1.5)
Neut_3Mvs1Y <- FindMarkers(object = pseudo_exp,
                            ident.1 = 'Neutrophils-HIV+ ART<3 Months',
                            ident.2 = 'Neutrophils-HIV+ ART>1 Year',
                            #min.pct = 0.25,
                            test.use = 'DESeq2')
Neut_3Mvs1Y <- Neut_3Mvs1Y %>%
  filter(p_val<0.05, abs(avg_log2FC)>1.5)
# Gene set enrichment analysis to determining the important genes in pathways
Neut_diff_exp <- unique(rbind(Neut_3MvsNeg,
                              Neut_1YvsNeg,
                              Neut_3Mvs1Y))
Neut_diff_exp <- Neut_diff_exp %>% 
  tibble::rownames_to_column(var = "gene")
# Rank differentially expressed genes
Neut_ranked_genes <- Neut_diff_exp %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::mutate(ranked=avg_log2FC) %>%
  dplyr::select(gene,ranked) %>%
  dplyr::pull(ranked,gene) 
str(Neut_ranked_genes)
# Get hallmark gene sets
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  as.data.frame()
# Prepare gene sets for GSEA
gene_sets <- split(m_t2g$gene_symbol,m_t2g$gs_name)
str(gene_sets)
Neut_ranked_genes <- as.numeric(Neut_ranked_genes)
names(Neut_ranked_genes) <- Neut_diff_exp$gene
common_genes <- intersect(names(Neut_ranked_genes),unlist(gene_sets))
Neut_ranked_genes <- Neut_ranked_genes[common_genes]
# Perform GSEA with gene symbols
neut_gsea_results <- GSEA(Neut_ranked_genes,
                          TERM2GENE = gene_sets,
                          pvalueCutoff = 0.05)
