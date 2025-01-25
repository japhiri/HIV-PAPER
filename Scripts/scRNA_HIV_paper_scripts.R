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
                        "magrittr",  "mvtnorm", "zoo", "patchwork", "mgcv", "PropCIs", "writexl", "presto",
                        "ggsignif", "ggpubr", "ggeasy", "cowplot","ggExtra", "PupillometryR","hrbrthemes", "ggstance",
                        "survival","survminer","sysfonts","showtext","nlme",'glue'))

# install.packages("devtools")
devtools::install_github("saeyslab/nichenetr")
devtools::install_github("saeyslab/multinichenetr")
devtools::install_github('immunogenomics/presto')


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


# Load all cell cluster ----
load("data/Single_Cell_Data/all_merged_subset_labelled.RData")
DimPlot(all_merged_subset_labelled,
        reduction = "umap.harmony",label = T,repel = T)+NoLegend()+
  labs(x="UMAP1",y="UMAP2")+
  theme_minimal()
# Renname idents
all_merged_subset_labelled <- RenameIdents(object = all_merged_subset_labelled,
                                           "Secretory cells" = "Secretory cells",
                                           "Goblet cells" = "Goblet cells",
                                           "CD3+ T cells" = "CD3+ T cells",
                                           "Phagocytes" = "Mono/Mac",
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
all_merged_subset_labelled$Clusters <- paste0(all_merged_subset_labelled@active.ident)
DimPlot(all_merged_subset_labelled,
        reduction = "umap.harmony",
        label = T,repel = T,
        #split.by = "HIV_Status"
        )+NoLegend()+
  theme_minimal()+
  labs(x="UMAP-1",y="UMAP-2")+
  theme(legend.position = 'none',
        axis.title = element_text(face = 'bold'),
        axis.ticks = element_blank(),
        axis.text = element_blank())

# Define a consistent color palette for 17 clusters
Clusters <- unique(all_merged_subset_labelled$Clusters)
num_clusters <- 11
base_palette <- brewer.pal(min(num_clusters,11),'Paired')
#additional_colors <- colorRampPalette(brewer.pal(8,'Set3'))(num_clusters-12)
#color_palette <- c(base_palette,additional_colors)
#names(color_palette) <- Clusters
names(base_palette) <- Clusters
# Plot UMAP with consistent colors
# Figure 3a (Main UMAP) ----
main_umap <- DimPlot(all_merged_subset_labelled,
        reduction = 'umap.harmony',
        label = TRUE,
        repel = TRUE,
        label.size = 5)+
  scale_color_manual(values = base_palette)+
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
# Save Figure 3a
ggsave(main_umap,filename="Figures/Figure 3/Fig3a.png",
       width = 12,height = 9,dpi = 1080,units = "in")

ggsave(main_umap,filename="Figures/Figure 3/Fig3a.pdf",
       width = 12,height = 9,dpi = 1080, units = "in")


# Figure 3b (Pie chart)----
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
  scale_fill_manual(values = base_palette)+
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void()+
  theme(legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = 'bold'))+
  guides(fill=guide_legend(override.aes = list(shape=21, size=5,color=NA)))
all_data_pie

# Save Figure 3b
ggsave(all_data_pie,filename="Figures/Figure 3/Fig3b.png",
       width = 20,height = 7,dpi = 1080,units = "in")

ggsave(all_data_pie,filename="Figures/Figure 3/Fig3b.pdf",
       width = 20,height = 7,dpi = 1080, units = "in")

# Figure 3c (Heatmap) -----
# Heatmap for differential expressed gened for each cluster 
all_merged_subset_labelled <- FindVariableFeatures(all_merged_subset_labelled,
                                                   selection.method = 'vst',
                                                   verbose = TRUE)
genes <- VariableFeatures(all_merged_subset_labelled)
toplot <- SeuratExtend::CalcStats(all_merged_subset_labelled,
                    features = genes,
                    method = 'zscore',
                    order = 'p',
                    n = 5)
cluster_heatmap <-SeuratExtend::Heatmap(toplot,
                                        lab_fill='zscore',
                                        plot.margin = margin(l = 30),
                                        angle = 90)    
cluster_heatmap <- as.ggplot(cluster_heatmap)
cluster_heatmap
# Save Figure 3c
ggsave(cluster_heatmap,filename="Figures/Figure 3/Fig3c.png",
       width = 4,height = 8,dpi = 1080,units = "in")

ggsave(cluster_heatmap,filename="Figures/Figure 3/Fig3c.pdf",
       width = 4,height = 8,dpi = 1080, units = "in")

# Figure 3d (Main Bargraph)----
# Calculate cluster proportion 
cluster_df <- as.data.frame(table(all_merged_subset_labelled$Clusters,
                                  all_merged_subset_labelled$HIV_Status))
colnames(cluster_df)<-c('Cluster','HIV Status','Count')

bargraph_main <- cluster_df %>%
  ggplot(aes(x=`HIV Status`,
             y=Count, fill=factor(Cluster,levels=c("Secretory cells","Goblet cells",
                                           "CD3+ T cells","Mono/Mac","B cells",
                                           "Ciliated cells","Neutrophils","Squamous cells",
                                           "Ionocytes","Dendritic cells","Basal cells"))))+
  geom_col(position = "fill")+
  scale_fill_manual(values = base_palette)+
  labs(x='',
       y='Frequency of cells')+
  theme_classic()+
  theme(legend.position = "right",
        legend.key.size = unit(0.7,'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(angle = 0,hjust = 1,face = 'bold'))+
  coord_flip()
bargraph_main

# Save Figure 3d
ggsave(bargraph_main,filename="Figures/Figure 3/Fig3d.png",
       width = 10,height = 5,dpi = 1080,units = "in")

ggsave(bargraph_main,filename="Figures/Figure 3/Fig3d.pdf",
       width = 10,height = 5,dpi = 1080, units = "in")

# Perform GO analysis for each cluster -----
# Load DEG markers for each cluster
Markers <- FindAllMarkers(
  all_merged_subset_labelled,
  assay = "RNA",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.01,
  verbose = TRUE,
  only.pos = TRUE,
  min.cells.feature = 3,
  min.cells.group = 3)

write.csv(Markers,'scRNAseq_Results/Markers.csv',row.names = T)
#Markers <- read_csv('scRNAseq_Results/Markers.csv',row.name=T)
  #column_to_rownames(Markers)

# Filter for significant DEGs (adjusted p-value <0.05)
significant_markers <- Markers %>%
dplyr::filter(p_val_adj <= 0.05)


# Split DEGs by cluster
DEGs_by_cluster <- split(significant_markers$gene,
significant_markers$cluster)

# Perform actual GO analysis
run_go_analysis <- function(genes, cluster_name){
    enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
    )
}
# Perform GO analysis for all clusters
go_results <- lapply(names(DEGs_by_cluster),function(cluster){
    cat("Running GO analysis for cluster:", cluster, "\n")
    run_go_analysis(DEGs_by_cluster[[cluster]], cluster)
})

names(go_results) <- names(DEGs_by_cluster)

# Extract the @result slot and process data
zscore_data <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::select(Description, zScore)
      rownames(data) <- data$Description
      data <- data %>%
        dplyr::select(-Description)
      return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(go_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(go_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(go_results)

# Subset the go_results object
go_results_top20 <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::arrange(p.adjust) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 5)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
go_results_top20 <- Filter(Negate(is.null), go_results_top20)

# Check the first cluster's top 20 for example
head(go_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(go_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

GO_Terms_Heatmap <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
    #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
    #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
GO_Terms_Heatmap
GO_Terms_Heatmap <- as.ggplot(GO_Terms_Heatmap)

# Save Figure 3e
ggsave(GO_Terms_Heatmap,filename="Figures/Figure 3/Fig3e.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(GO_Terms_Heatmap,filename="Figures/Figure 3/Fig3e.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
GO_Terms_Heatmap <- SeuratExtend::Heatmap(heatmap_matrix,
                                          lab_fill='zscore',
                                          plot.margin = margin(l = 30),
                                          y_text_position = "right",
                                          angle = 90,
                                          hjust = 0.5,
                                          color_scheme = color_palette,
                                          vjust = 0.5)
GO_Terms_Heatmap <- as.ggplot(GO_Terms_Heatmap)
GO_Terms_Heatmap
# Save Figure 3e
ggsave(GO_Terms_Heatmap,filename="Figures/Figure 3/Fig3e2.png",
       width = 11.5,height = 12,dpi = 1080,units = "in")

ggsave(GO_Terms_Heatmap,filename="Figures/Figure 3/Fig3e2.pdf",
       width = 11.5,height = 12,dpi = 1080, units = "in")

# Dotplot of genes ------
# Interested in immune genes from top20 GO terms
# Extracting immune genes from top20 GO terms

top20_go_terms <- rownames(heatmap_matrix) %>%
  as.data.frame()
colnames(top20_go_terms)<-c("Go_terms")

rownames(top20_go_terms) <- top20_go_terms$Go_terms %>%
  as.matrix()

# Extract the top20 Go terms
# Extract the @result slot and process data
geneID_data <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::select(Description,geneID)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
geneID_matrix <- geneID_data[[1]]
colnames(geneID_matrix) <- names(geneID_data)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(geneID_data)) {
  current_cluster <- geneID_data[[i]]
  colname <- names(go_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  geneID_matrix <- merge(geneID_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(geneID_matrix) <- geneID_matrix$Row.names
  geneID_matrix <- geneID_matrix[, -1]  # Drop redundant Row.names column
}
geneID_matrix <- as.matrix(geneID_matrix)

# Match rownames and subset
top20_geneID_matrix <- geneID_matrix[match(rownames(top20_go_terms), rownames(geneID_matrix), nomatch = 0), ]


top20_geneID_matrix %>% as.data.frame(top20_geneID_matrix)

top20_genes <- lapply(top20_geneID_matrix[,1:11], as.vector)
top20_genes <- lapply(top20_genes,unique)

top20_genes <- unique(unlist(top20_genes)) %>%
  as.vector()
top20_genes <- strsplit(top20_genes,split = "/")
top20_genes <- unlist(top20_genes) %>%
  as.data.frame()

colnames(top20_genes)<-"Genes"

DotPlot(all_merged_subset_labelled,
        features = unique(top20_genes$Genes),
        #split.by = "HIV_Status",
        #cols = c("red","white","blue")
        )+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_blank())

# Running GSEA on all clusters -----
# Split DEGs by cluster
DEGs_by_cluster <- split(significant_markers$gene,
                         significant_markers$cluster)
# Perform actual GO analysis
organism <- org.Hs.eg.db
run_gsea_analysis <- function(genes, cluster_name){
  gseGO(
    geneList = genes,
    ont ="BP", 
    keyType = "SYMBOL", 
    nPerm = 10000, 
    minGSSize = 3, 
    maxGSSize = 800, 
    pvalueCutoff = 0.05, 
    verbose = TRUE, 
    OrgDb = organism, 
    pAdjustMethod = "none"
  )
}
# Perform GO analysis for all clusters
gsea_results <- lapply(names(DEGs_by_cluster),function(cluster){
  cat("Running gsea analysis for cluster:", cluster, "\n")
  run_go_analysis(DEGs_by_cluster[[cluster]], cluster)
})

names(gsea_results) <- names(DEGs_by_cluster)

# Extract the @result slot and process data
foldEnrichment_data <- lapply(gsea_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::select(Description, FoldEnrichment)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
foldEnrichment_matrix <- foldEnrichment_data[[1]]
colnames(foldEnrichment_matrix) <- names(gsea_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(foldEnrichment_data)) {
  current_cluster <- foldEnrichment_data[[i]]
  colname <- names(gsea_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  foldEnrichment_matrix <- merge(foldEnrichment_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(foldEnrichment_matrix) <- foldEnrichment_matrix$Row.names
  foldEnrichment_matrix <- foldEnrichment_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
foldEnrichment_matrix[is.na(foldEnrichment_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(foldEnrichment_matrix) <- names(gsea_results)

# Subset the go_results object
gsea_results_top20 <- lapply(gsea_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::arrange(p.adjust) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 5)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
gsea_results_top20 <- Filter(Negate(is.null), gsea_results_top20)

# Check the first cluster's top 20 for example
head(gsea_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_gsea_results <- dplyr::bind_rows(gsea_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,FoldEnrichment)
rownames(merged_top20_gsea_results)<-NULL

merged_top20_gsea_results <- merged_top20_gsea_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = FoldEnrichment)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_gsea_results[is.na(merged_top20_gsea_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_gsea_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_gsea_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

GO_Terms_Heatmap <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
GO_Terms_Heatmap
GO_Terms_Heatmap <- as.ggplot(GO_Terms_Heatmap)

# Save Figure 3e
ggsave(GO_Terms_Heatmap,filename="Figures/Figure 3/Fig3e.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(GO_Terms_Heatmap,filename="Figures/Figure 3/Fig3e.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")

# Running GO analysis on all clusters pre ranking the most important GO terms -----
# Split DEGs by cluster
DEGs_by_cluster <- split(significant_markers$gene,
                         significant_markers$cluster)
# Perform actual GO analysis
organism <- org.Hs.eg.db
run_go_analysis <- function(genes, cluster_name){
  enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
}
# Perform GO analysis for all clusters
go_results <- lapply(names(DEGs_by_cluster),function(cluster){
  cat("Running gsea analysis for cluster:", cluster, "\n")
  run_go_analysis(DEGs_by_cluster[[cluster]], cluster)
})

names(go_results) <- names(DEGs_by_cluster)

# Extract the @result slot and process data
zscore_data <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(go_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(go_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                                 by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(go_results)

# Extract the @result slot and process data
pvalue_data <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::select(Description, p.adjust,FoldEnrichment)
    rownames(data) <- data$Description
    data <- data %>%
      #dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
pvalue_matrix <- pvalue_data[[1]]
colnames(pvalue_matrix) <- names(go_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(pvalue_data)) {
  current_cluster <- pvalue_data[[i]]
  colname <- names(go_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  pvalue_matrix <- merge(pvalue_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(pvalue_matrix) <- pvalue_matrix$Row.names
  pvalue_matrix <- pvalue_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
pvalue_matrix[is.na(pvalue_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(pvalue_matrix) <- names(go_results)
# Pivot longer
pvalue_matrix_longer <- pvalue_matrix %>%
  mutate(GO_terms = rownames(pvalue_matrix)) %>%
  pivot_longer(cols = c("Secretory cells","Goblet cells","CD3+ T cells","Mono/Mac",
                        "B cells","Ciliated cells","Neutrophils","Squamous cells",
                        "Ionocytes","Dendritic cells","Basal cells"),
               names_to = "cluster",
               values_to = "p.adjust") %>%
  as.data.frame()

pvalue_matrix_longer <- pvalue_matrix_longer %>%
  dplyr::arrange(p.adjust)
top50pvalue_matrix_longer <- head(pvalue_matrix_longer,1000)
top50pvalue_matrix_longer <- top50pvalue_matrix_longer[!duplicated(top50pvalue_matrix_longer$GO_terms), ]  

heatmap_matrix <- pvalue_matrix %>%
  filter(rownames(pvalue_matrix) %in% top50pvalue_matrix_longer$GO_terms) %>%
  as.matrix() %>%
  head(50)


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

GO_Terms_Heatmap <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  #col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("blue", "red")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
GO_Terms_Heatmap
GO_Terms_Heatmap <- as.ggplot(GO_Terms_Heatmap)

# Save Figure 3e
ggsave(GO_Terms_Heatmap,filename="Figures/Figure 3/Fig3e.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(GO_Terms_Heatmap,filename="Figures/Figure 3/Fig3e.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")


GO_Terms_Heatmap <-SeuratExtend::Heatmap(heatmap_matrix,
                                        lab_fill='zscore',
                                        plot.margin = margin(l = 30),
                                        angle = 90,
                                        cluster_rows=T)  
GO_Terms_Heatmap
GO_Terms_Heatmap <- as.ggplot(GO_Terms_Heatmap)

# Save Figure 3e
ggsave(GO_Terms_Heatmap,filename="Figures/Figure 3/Fig3e3.png",
       width = 7,height = 7,dpi = 1080,units = "in")

ggsave(GO_Terms_Heatmap,filename="Figures/Figure 3/Fig3e3.pdf",
       width = 7,height = 7,dpi = 1080, units = "in")

# Filetering based on p adjusted and returning the top-50 based on FoldEnrichment values -----
# Example: Assuming `enrichGO_list` is your list of enrichGO results for each cell cluster
filtered_results <- lapply(go_results, function(cluster_results) {
  cluster_results@result %>% 
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::select(Description, FoldEnrichment)
})

# Combine results into a single dataframe for creating a matrix
combined_results <- bind_rows(filtered_results, .id = "Cluster") %>%
  arrange(desc(FoldEnrichment)) %>%
  head(50)
  
  # Create a matrix with GO terms as rows and clusters as columns
  heatmap_matrix <- combined_results %>%
    pivot_wider(names_from = Cluster, values_from = FoldEnrichment, values_fill = NA) %>%
    column_to_rownames(var = "Description") %>%
    as.matrix()
  
  heatmap_matrix[is.na(heatmap_matrix)] <- 0
  
  # Plot the heatmap
  pheatmap(heatmap_matrix,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           cluster_rows = T,
           cluster_cols = F,
           main = "GO Term Fold Enrichment")  
  
# GO-terms in Neutrophils grouped by HIV status----
Neutrophils <- subset(all_merged_subset_labelled,
                      idents="Neutrophils",
                      invert=F)
# Make HIV status as the Ident
Idents(Neutrophils)<-"HIV_Status"
# Find differentially expressed markers
Neutrophil_markers <- FindAllMarkers(
  Neutrophils,
  assay = "RNA",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.01,
  verbose = TRUE,
  only.pos = TRUE,
  min.cells.feature = 3,
  min.cells.group = 3)
# Filter for significant markers
DEG_Neutrophils <- Neutrophil_markers %>%
  dplyr::filter(p_val_adj<=0.05)

# Split DEGs by HIV status
Neut_DEGs_by_HIV_status <- split(DEG_Neutrophils$gene,
                                 DEG_Neutrophils$cluster)

# Perform actual GO analysis
run_Neut_go_analysis <- function(genes, cluster_name){
  enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
}
# Perform GO analysis for all clusters
Neut_go_results <- lapply(names(Neut_DEGs_by_HIV_status),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_Neut_go_analysis(Neut_DEGs_by_HIV_status[[cluster]], cluster)
})

names(Neut_go_results) <- names(Neut_DEGs_by_HIV_status)

# Extract the @result slot and process data
zscore_data <- lapply(Neut_go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(Neut_go_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(Neut_go_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(Neut_go_results)

# Subset the go_results object
Neut_go_results_top20 <- lapply(Neut_go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::arrange(p.adjust) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 20)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
Neut_go_results_top20 <- Filter(Negate(is.null), Neut_go_results_top20)

# Check the first cluster's top 20 for example
head(Neut_go_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(Neut_go_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

GO_Terms_Heatmap <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
GO_Terms_Heatmap
GO_Terms_Heatmap <- as.ggplot(GO_Terms_Heatmap)

# Save Figure 3e
ggsave(GO_Terms_Heatmap,filename="Figures/Figure 3/Fig3e4.png",
       width = 4,height = 8,dpi = 1080,units = "in")

ggsave(GO_Terms_Heatmap,filename="Figures/Figure 3/Fig3e4.pdf",
       width = 4,height = 8,dpi = 1080, units = "in")

GO_Terms_Heatmap <- SeuratExtend::Heatmap(
  heatmap_matrix,
  plot.margin = margin(l = 30),
  angle=90
)
GO_Terms_Heatmap
# Save Figure 3e
ggsave(GO_Terms_Heatmap,filename="Figures/Figure 3/Fig3e4.png",
       width = 8,height = 12,dpi = 1080,units = "in")

ggsave(GO_Terms_Heatmap,filename="Figures/Figure 3/Fig3e4.pdf",
       width = 8,height = 12,dpi = 1080, units = "in")
# CEll TO CELL COMMUNICATION NETWORK USING MULTINICHENETR -----
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
library(nichenetr)
library(Seurat)
organism="human"
# Make gene names syntactically valid using make.names to avoid loss of genes and load the organism name
lr_network = readRDS("NicheNet/lr_network_human_allInfo_30112033.rds")%>%
  mutate(ligand = convert_alias_to_symbols(ligand, organism = organism),
         receptor = convert_alias_to_symbols(receptor, organism = organism)) %>% 
  mutate(ligand = make.names(ligand), receptor = make.names(receptor)) %>% 
  distinct(ligand, receptor)

ligand_target_matrix = readRDS("NicheNet/ligand_target_matrix_nsga2r_final.rds")
colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
  convert_alias_to_symbols(organism = organism) %>% make.names()
rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
  convert_alias_to_symbols(organism = organism) %>% make.names()

lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]

# Prepare the seurat object for conversion to single cell experiment object
all_merged_subset_labelled$Clusters <- paste0(all_merged_subset_labelled@active.ident)
all_merged_subset_labelled$HIV_Status <- paste0(all_merged_subset_labelled$HIV_Status)
all_merged_subset_labelled$sample <- paste0(all_merged_subset_labelled$sample)
# Make names syntactically valid
all_merged_subset_labelled$HIV_Status <- make.names(all_merged_subset_labelled$HIV_Status)
all_merged_subset_labelled$sample <- make.names(all_merged_subset_labelled$sample)
all_merged_subset_labelled$Clusters <- make.names(all_merged_subset_labelled$Clusters)
# Convert to a single cell experiment
sce <- Seurat::as.SingleCellExperiment(all_merged_subset_labelled)
SingleCellExperiment::colData(sce)
# Prepare cell-cell communication analysis
sample_id="sample"
group_id="HIV_Status"
celltype_id="Clusters"
covariates=NA
batches=NA
# Define contracts
contrasts_oi = c("'HIV..ART.3.Months-(HIV..ART.1.Year+HIV.)/2','HIV..ART.1.Year-(HIV..ART.3.Months+HIV.)/2','HIV.-(HIV..ART.1.Year+HIV..ART.3.Months)/2'")
contrast_tbl = tibble(contrast =
                        c("HIV..ART.3.Months-(HIV..ART.1.Year+HIV.)/2","HIV..ART.1.Year-(HIV..ART.3.Months+HIV.)/2","HIV.-(HIV..ART.1.Year+HIV..ART.3.Months)/2"),
                      group = c("HIV..ART.3.Months","HIV..ART.1.Year","HIV."))

senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in%
            c(senders_oi,receivers_oi)]

conditions_keep = c("HIV..ART.3.Months","HIV..ART.1.Year","HIV.")
sce = sce[,SummarizedExperiment::colData(sce)[,group_id] %in% 
            conditions_keep]
# Running Multinichenet core analysis
# Cell-type filtering
min_cells = 5
abundance_info =multinichenetr::get_abundance_info(
  sce = sce, 
  sample_id = sample_id, 
  group_id = group_id, 
  celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, 
  receivers_oi = receivers_oi, 
  batches = batches)

abundance_info$abund_plot_sample

abundance_df_summarized = abundance_info$abundance_data %>%
  dplyr::mutate(keep = as.logical(keep)) %>%
  dplyr::group_by(group_id,celltype_id) %>%
  dplyr::summarise(samples_present = sum((keep)))

celltypes_absent_one_condition = abundance_df_summarized %>%
  dplyr::filter(samples_present == 0) %>% 
  dplyr::pull(celltype_id) %>%
  unique()

celltypes_present_one_condition = abundance_df_summarized %>%
  dplyr::filter(samples_present >= 2) %>%
  dplyr::pull(celltype_id) %>%
  unique()

total_nr_conditions = SummarizedExperiment::colData(sce)[,group_id] %>%
  unique() %>%
  length()

absent_celltypes = abundance_df_summarized %>% 
  filter(samples_present < 2) %>% 
  group_by(celltype_id) %>% 
  count() %>% 
  #filter(n == total_nr_conditions) %>% 
  pull(celltype_id)

analyse_condition_specific_celltypes = FALSE
if(analyse_condition_specific_celltypes == TRUE){
  senders_oi = senders_oi %>% setdiff(absent_celltypes)
  receivers_oi = receivers_oi %>% setdiff(absent_celltypes)
} else {
  senders_oi = senders_oi %>% 
    setdiff(union(absent_celltypes, condition_specific_celltypes))
  receivers_oi = receivers_oi %>% 
    setdiff(union(absent_celltypes, condition_specific_celltypes))
}

sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)]
# Gene filtering
min_sample_prop = 0.50
fraction_cutoff = 0.05


frq_list = multinichenetr::get_frac_exprs(
  sce = sce,
  sample_id = sample_id,
  celltype_id = celltype_id,
  group_id = group_id,
  batches = batches,
  min_cells = min_cells,
  fraction_cutoff = fraction_cutoff,
  min_sample_prop = min_sample_prop)

genes_oi = frq_list$expressed_df %>%
  dplyr::filter(expressed == TRUE) %>%
  dplyr::pull(gene) %>% unique()

sce = sce[genes_oi,]

# Pseudobulk expression calculation
abundance_expression_info = multinichenetr::process_abundance_expression_info(
  sce = sce,
  sample_id = sample_id,
  group_id = group_id,
  celltype_id = celltype_id,
  min_cells = min_cells,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network,
  batches = batches,
  frq_list = frq_list,
  abundance_info = abundance_info)

abundance_expression_info$celltype_info$pb_df %>% head()
abundance_expression_info$celltype_info$pb_df_group %>% head()
abundance_expression_info$sender_receiver_info$pb_df %>% head()
abundance_expression_info$sender_receiver_info$pb_df_group %>% head()

# Differential expression (DE) analysis
DE_info = multinichenetr::get_DE_info(
  sce = sce,
  sample_id = sample_id,
  group_id = group_id,
  celltype_id = celltype_id,
  batches = batches,
  covariates = covariates,
  contrasts_oi = contrasts_oi,
  min_cells = min_cells,
  expressed_df = frq_list$expressed_df)
# Check DE results
DE_info$celltype_de$de_output_tidy %>% head()
DE_info$hist_pvals

empirical_pval = FALSE
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
  celltype_de = DE_info_emp$de_output_tidy_emp %>% select(-p_val, -p_adj) %>% 
    rename(p_val = p_emp, p_adj = p_adj_emp)
} else {
  celltype_de = DE_info$celltype_de$de_output_tidy
} 

# Combine DE information for ligand-senders and receivers-receivers
sender_receiver_de = multinichenetr::combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network)
sender_receiver_de %>% head(20)

# ligand activity prediction
logFC_threshold = 0.50
p_val_threshold = 0.05
p_val_adj = FALSE


geneset_assessment = contrast_tbl$contrast %>%
  lapply(
    multinichenetr::process_geneset_data,
    celltype_de,
    logFC_threshold,
    p_val_adj,
    p_val_threshold) %>%
  bind_rows()
geneset_assessment

# In case we want to use the adjusted p-values as threshold
geneset_assessment_adjustedPval = contrast_tbl$contrast %>%
  lapply(
    multinichenetr::process_geneset_data,
    celltype_de,
    logFC_threshold,
    p_val_adj = TRUE,
    p_val_threshold) %>%
  bind_rows()
geneset_assessment_adjustedPval
  
  
# Ligand activity analysis and ligand target inference
top_n_target = 250
verbose = TRUE
cores_system = 8
n.cores = min(cores_system, celltype_de$cluster_id %>% unique() %>% length())

ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  multinichenetr::get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = verbose,
    n.cores = n.cores)))  
  
  
ligand_activities_targets_DEgenes$ligand_activities %>% head(20)

# Prioritization: Ranking cell-cell communication patters through multi criteria prioritization
ligand_activity_down = FALSE

sender_receiver_tb
  
  
  
  
  
  
  
  
  
  
  




























