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
                        "DelayedArray","DelayedMatrixStats","limma","lme4","batchelor","HDF5Array","terra","ggrastr","nloptr","ggsignif"))





remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)
remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE)


install.packages("azimuth")

options(Seurat.object.assay.version = "v5")


# Create a Seurat object for each sample using the for loop command ----

for (file in c("CUH124_raw_feature_bc_matrix",
               "CUG11X_raw_feature_bc_matrix",
               "CUF131_raw_feature_bc_matrix",
               "CUF130_raw_feature_bc_matrix",
               "CUF134_raw_feature_bc_matrix",
               "CUF13I_raw_feature_bc_matrix",
               "CUF135_raw_feature_bc_matrix",
               "CUF136_raw_feature_bc_matrix",
               "CUF13J_raw_feature_bc_matrix",
               "CUF13K_raw_feature_bc_matrix",
               "CUF137_raw_feature_bc_matrix",
               "CUF12Y_raw_feature_bc_matrix"
               )){
  seurat_data <- Read10X(data.dir = paste0("data/",# Directory containing the _feature_bc_matrix folders
  file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.features = 30,# only load droplets that express atleast 30 genes because they are likely cells
                                   project = file)
  assign(file, seurat_obj)
}

save(CUF12Y_raw_feature_bc_matrix,file = "data/CUF12Y_raw_feature_bc_matrix.RData")
saveRDS(CUF12Y_raw_feature_bc_matrix,file = "data/CUF12Y_raw_feature_bc_matrix.rds")

# Create a merged Seurat object ----
all_merged_seurat <- merge(x=CUF131_raw_feature_bc_matrix,
                       y=c(CUH124_raw_feature_bc_matrix,
                           CUG11X_raw_feature_bc_matrix,
                           CUF130_raw_feature_bc_matrix,
                           CUF134_raw_feature_bc_matrix,
                           CUF13I_raw_feature_bc_matrix,
                           CUF135_raw_feature_bc_matrix,
                           CUF136_raw_feature_bc_matrix,
                           CUF13J_raw_feature_bc_matrix,
                           CUF13K_raw_feature_bc_matrix,
                           CUF137_raw_feature_bc_matrix,
                           CUF12Y_raw_feature_bc_matrix
                           ),
                       add.cell.id = c("CUF131",
                                       "CUH124",
                                       "CUG11X",
                                       "CUF130",
                                       "CUF134",
                                       "CUF13I",
                                       "CUF135",
                                       "CUF136",
                                       "CUF13J",
                                       "CUF13K",
                                       "CUF137",
                                       "CUF12Y"
                                       ))



VlnPlot(all_merged_seurat, features = c("mitoRatio","nFeature_RNA"),
        pt.size = 0)+NoLegend()

view(all_merged_seurat@meta.data)


# Add number of genes per UMI for each cell to metadata----
all_merged_seurat$log10GenesPerUMI <- log10(all_merged_seurat$nFeature_RNA) / log10(all_merged_seurat$nCount_RNA)

# Compute percent mito ratio ----
all_merged_seurat$mitoRatio <- PercentageFeatureSet(object = all_merged_seurat, pattern = "^MT-")
all_merged_seurat$mitoRatio <- all_merged_seurat@meta.data$mitoRatio / 100
# Compute percent ribosomal protein----
all_merged_seurat$Percent_ribo<-PercentageFeatureSet(object = all_merged_seurat, pattern = "^RP[SL]")
all_merged_seurat$Percent_ribo <- all_merged_seurat@meta.data$Percent_ribo / 100

# compute percent of hemoglobin genes----
all_merged_seurat$Percent_hemoglobin<-PercentageFeatureSet(object = all_merged_seurat, pattern = "^HB[^(P)]")
all_merged_seurat$Percent_hemoglobin <- all_merged_seurat@meta.data$Percent_hemoglobin/ 100


# Create metadata dataframe ----
metadata <- all_merged_seurat@meta.data

# Add cell IDs to metadata---
metadata$cells <- rownames(metadata)




# Create sample column ----
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^CUH124"))] <- "CUH124"
metadata$sample[which(str_detect(metadata$cells, "^CUG11X"))] <- "CUG11X"
metadata$sample[which(str_detect(metadata$cells, "^CUF131"))] <- "CUF131"
metadata$sample[which(str_detect(metadata$cells, "^CUF130"))] <- "CUF130"
metadata$sample[which(str_detect(metadata$cells, "^CUF134"))] <- "CUF134"
metadata$sample[which(str_detect(metadata$cells, "^CUF13I"))] <- "CUF13I"
metadata$sample[which(str_detect(metadata$cells, "^CUF135"))] <- "CUF135"
metadata$sample[which(str_detect(metadata$cells, "^CUF136"))] <- "CUF136"
metadata$sample[which(str_detect(metadata$cells, "^CUF13J"))] <- "CUF13J"
metadata$sample[which(str_detect(metadata$cells, "^CUF13K"))] <- "CUF13K"
metadata$sample[which(str_detect(metadata$cells, "^CUF137"))] <- "CUF137"
metadata$sample[which(str_detect(metadata$cells, "^CUF12Y"))] <- "CUF12Y"






# Create HIV status column----
metadata$HIV_Status <- NA
metadata$HIV_Status[which(str_detect(metadata$cells, "^CUH124"))] <- "HIV+ ART>1 Year"
metadata$HIV_Status[which(str_detect(metadata$cells, "^CUG11X"))] <- "HIV+ ART>1 Year"
metadata$HIV_Status[which(str_detect(metadata$cells, "^CUF131"))] <- "HIV-"
metadata$HIV_Status[which(str_detect(metadata$cells, "^CUF130"))] <- "HIV+ ART>1 Year"
metadata$HIV_Status[which(str_detect(metadata$cells, "^CUF134"))] <- "HIV+ ART<3 Months"
metadata$HIV_Status[which(str_detect(metadata$cells, "^CUF13I"))] <- "HIV+ ART>1 Year"
metadata$HIV_Status[which(str_detect(metadata$cells, "^CUF135"))] <- "HIV+ ART>1 Year"
metadata$HIV_Status[which(str_detect(metadata$cells, "^CUF136"))] <- "HIV-"
metadata$HIV_Status[which(str_detect(metadata$cells, "^CUF13J"))] <- "HIV+ ART<3 Months"
metadata$HIV_Status[which(str_detect(metadata$cells, "^CUF13K"))] <- "HIV+ ART<3 Months"
metadata$HIV_Status[which(str_detect(metadata$cells, "^CUF137"))] <- "HIV-"
metadata$HIV_Status[which(str_detect(metadata$cells, "^CUF12Y"))] <- "HIV+ ART>1 Year"


# Create Carriage status column----
metadata$Carriage_Status <- NA
metadata$Carriage_Status[which(str_detect(metadata$cells, "^CUH124"))] <- "SPN+"
metadata$Carriage_Status[which(str_detect(metadata$cells, "^CUG11X"))] <- "SPN-"
metadata$Carriage_Status[which(str_detect(metadata$cells, "^CUF131"))] <- "SPN+"
metadata$Carriage_Status[which(str_detect(metadata$cells, "^CUF130"))] <- "SPN-"
metadata$Carriage_Status[which(str_detect(metadata$cells, "^CUF134"))] <- "SPN+"
metadata$Carriage_Status[which(str_detect(metadata$cells, "^CUF13I"))] <- "SPN+"
metadata$Carriage_Status[which(str_detect(metadata$cells, "^CUF135"))] <- "SPN+"
metadata$Carriage_Status[which(str_detect(metadata$cells, "^CUF136"))] <- "SPN+"
metadata$Carriage_Status[which(str_detect(metadata$cells, "^CUF13J"))] <- "SPN+"
metadata$Carriage_Status[which(str_detect(metadata$cells, "^CUF13K"))] <- "SPN-"
metadata$Carriage_Status[which(str_detect(metadata$cells, "^CUF137"))] <- "SPN+"
metadata$Carriage_Status[which(str_detect(metadata$cells, "^CUF12Y"))] <- "SPN+"

# Create Carriage density column ----
metadata$Carriage_density <- NA
metadata$Carriage_density[which(str_detect(metadata$cells, "^CUH124"))] <- 5862500
metadata$Carriage_density[which(str_detect(metadata$cells, "^CUG11X"))] <- NA
metadata$Carriage_density[which(str_detect(metadata$cells, "^CUF131"))] <- 26800000
metadata$Carriage_density[which(str_detect(metadata$cells, "^CUF130"))] <- NA
metadata$Carriage_density[which(str_detect(metadata$cells, "^CUF134"))] <- 686750
metadata$Carriage_density[which(str_detect(metadata$cells, "^CUF13I"))] <- 670
metadata$Carriage_density[which(str_detect(metadata$cells, "^CUF135"))] <- 33500
metadata$Carriage_density[which(str_detect(metadata$cells, "^CUF136"))] <- 67
metadata$Carriage_density[which(str_detect(metadata$cells, "^CUF13J"))] <- 26800
metadata$Carriage_density[which(str_detect(metadata$cells, "^CUF13K"))] <- NA
metadata$Carriage_density[which(str_detect(metadata$cells, "^CUF137"))] <- 10720
metadata$Carriage_density[which(str_detect(metadata$cells, "^CUF12Y"))] <- 2345



# Create Carriage serotype column ----
metadata$serotype <- NA
metadata$serotype[which(str_detect(metadata$cells, "^CUH124"))] <- "NVT"
metadata$serotype[which(str_detect(metadata$cells, "^CUG11X"))] <- NA
metadata$serotype[which(str_detect(metadata$cells, "^CUF131"))] <- "3"
metadata$serotype[which(str_detect(metadata$cells, "^CUF130"))] <- NA
metadata$serotype[which(str_detect(metadata$cells, "^CUF134"))] <- "9"
metadata$serotype[which(str_detect(metadata$cells, "^CUF13I"))] <- "1"
metadata$serotype[which(str_detect(metadata$cells, "^CUF135"))] <- "10"
metadata$serotype[which(str_detect(metadata$cells, "^CUF136"))] <- "NVT"
metadata$serotype[which(str_detect(metadata$cells, "^CUF13J"))] <- "10"
metadata$serotype[which(str_detect(metadata$cells, "^CUF13K"))] <- NA
metadata$serotype[which(str_detect(metadata$cells, "^CUF137"))] <- "NVT"
metadata$serotype[which(str_detect(metadata$cells, "^CUF12Y"))] <- "NVT"

# Create AGE column ----
metadata$age <- NA
metadata$age[which(str_detect(metadata$cells, "^CUH124"))] <- 35
metadata$age[which(str_detect(metadata$cells, "^CUG11X"))] <- 21
metadata$age[which(str_detect(metadata$cells, "^CUF131"))] <- 28
metadata$age[which(str_detect(metadata$cells, "^CUF130"))] <- 34
metadata$age[which(str_detect(metadata$cells, "^CUF134"))] <- 26
metadata$age[which(str_detect(metadata$cells, "^CUF13I"))] <- 27
metadata$age[which(str_detect(metadata$cells, "^CUF135"))] <- 34
metadata$age[which(str_detect(metadata$cells, "^CUF136"))] <- 36
metadata$age[which(str_detect(metadata$cells, "^CUF13J"))] <- 24
metadata$age[which(str_detect(metadata$cells, "^CUF13K"))] <- 27
metadata$age[which(str_detect(metadata$cells, "^CUF137"))] <- 28
metadata$age[which(str_detect(metadata$cells, "^CUF12Y"))] <- 25



# Create SEX column ----
metadata$sex <- NA
metadata$sex[which(str_detect(metadata$cells, "^CUH124"))] <- "Female"
metadata$sex[which(str_detect(metadata$cells, "^CUG11X"))] <- "Female"
metadata$sex[which(str_detect(metadata$cells, "^CUF131"))] <- "Female"
metadata$sex[which(str_detect(metadata$cells, "^CUF130"))] <- "Female"
metadata$sex[which(str_detect(metadata$cells, "^CUF134"))] <- "Male"
metadata$sex[which(str_detect(metadata$cells, "^CUF13I"))] <- "Male"
metadata$sex[which(str_detect(metadata$cells, "^CUF135"))] <- "Male"
metadata$sex[which(str_detect(metadata$cells, "^CUF136"))] <- "Female"
metadata$sex[which(str_detect(metadata$cells, "^CUF13J"))] <- "Female"
metadata$sex[which(str_detect(metadata$cells, "^CUF13K"))] <- "Female"
metadata$sex[which(str_detect(metadata$cells, "^CUF137"))] <- "Male"
metadata$sex[which(str_detect(metadata$cells, "^CUF12Y"))] <- "Male"


# Create viral load column ----
metadata$Viral_Load <- NA
metadata$Viral_Load[which(str_detect(metadata$cells, "^CUH124"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$cells, "^CUG11X"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$cells, "^CUF131"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$cells, "^CUF130"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$cells, "^CUF134"))] <- 312500
metadata$Viral_Load[which(str_detect(metadata$cells, "^CUF13I"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$cells, "^CUF135"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$cells, "^CUF136"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$cells, "^CUF13J"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$cells, "^CUF13K"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$cells, "^CUF137"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$cells, "^CUF12Y"))] <- "Not Detected"


# Create absolute CD4 count column ----
metadata$abs_CD4_count <- NA
metadata$abs_CD4_count[which(str_detect(metadata$cells, "^CUH124"))] <- 805
metadata$abs_CD4_count[which(str_detect(metadata$cells, "^CUG11X"))] <- 488
metadata$abs_CD4_count[which(str_detect(metadata$cells, "^CUF131"))] <- 381
metadata$abs_CD4_count[which(str_detect(metadata$cells, "^CUF130"))] <- 535
metadata$abs_CD4_count[which(str_detect(metadata$cells, "^CUF134"))] <- 14
metadata$abs_CD4_count[which(str_detect(metadata$cells, "^CUF13I"))] <- 452
metadata$abs_CD4_count[which(str_detect(metadata$cells, "^CUF135"))] <- 274
metadata$abs_CD4_count[which(str_detect(metadata$cells, "^CUF136"))] <- 1179
metadata$abs_CD4_count[which(str_detect(metadata$cells, "^CUF13J"))] <- 689
metadata$abs_CD4_count[which(str_detect(metadata$cells, "^CUF13K"))] <- 907
metadata$abs_CD4_count[which(str_detect(metadata$cells, "^CUF137"))] <- 431
metadata$abs_CD4_count[which(str_detect(metadata$cells, "^CUF12Y"))] <- 311



# rename columns ----
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
# Add metadata back to Seurat object ----
all_merged_seurat@meta.data <- metadata


View(all_merged_seurat@meta.data)

# Create .RData object to load at any time ----
save(all_merged_seurat, file = "data/all_merged_seurat.RData")
load("data/all_merged_seurat.RData")


all_merged_seurat<-JoinLayers(all_merged_seurat)

view(all_merged_seurat@meta.data)


# Filter out low quality cells using selected thresholds - these will change with experiment ----
all_filtered_seurat <- subset(x = all_merged_seurat,
                          subset= (nUMI >= 100) & 
                            (nUMI <= 50000) &
                            (nGene >= 50) & 
                            (log10GenesPerUMI > 0.70) & 
                            (Percent_hemoglobin<0.0001) &
                            (mitoRatio < 0.30))

length(all_filtered_seurat@meta.data$cells)



# Gene-level filtering ----

all_filtered_seurat<-JoinLayers(all_filtered_seurat)

# Extract counts
counts <- LayerData(object = all_filtered_seurat, layer = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 3

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
all_filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = all_filtered_seurat@meta.data)

length(all_filtered_seurat@meta.data$cells)
# PLOT THE QUALITY MATRIX----
# remove CUF12Y
new_all_filtered_seurat <- subset(all_filtered_seurat,
                                  cells = which(all_filtered_seurat$sample!='CUF12Y'))
items <- c('nGene','nUMI','mitoRatio')

vln_plots <- list()

for (features in items) {
  vln <- VlnPlot(new_all_filtered_seurat,
                 features = features,
                 pt.size = 0.0001,
                 alpha = 0.5)+
    NoLegend()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 8))
  vln_plots[[features]] <- vln
  
}

#arranged_plots <- grid.arrange(grobs = vln_plots,ncol=2,nrow=2)
arranged_plots <- grid.arrange(grobs = vln_plots,ncol=length(items))
# Save the quality matrix plots
ggsave('Thesis Figures/Quality.Matrix.png',arranged_plots,width = 10,height = 4,units = 'in')
ggsave('Thesis Figures/Quality.Matrix.pdf',arranged_plots,width = 10,height = 4,units = 'in')

nGenenUMI <- new_all_filtered_seurat@meta.data %>%
  ggplot(aes(nUMI,nGene, color=mitoRatio))+
  geom_point(size=0.1)+
  geom_smooth(method = 'nlm')+
  stat_cor(method = 'spearman',size=8,face='bold',label.y = 9000)+
  labs(x='nUMI',
       y='nGene')+
  scale_color_gradient()+
  theme_bw()+
  theme(axis.title = element_text(size = 16,face = 'bold'),
        panel.border = element_rect(linewidth = 0.1),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.1))
nGenenUMI
ggsave('Thesis Figures/nGenenUMI.png',nGenenUMI,width = 6,height = 6,units = 'in')
ggsave('Thesis Figures/nGenenUMI.pdf',nGenenUMI,width = 6,height = 6,units = 'in')


VlnPlot(all_filtered_seurat, features = c("mitoRatio"),
        pt.size = 0)+NoLegend()+
  theme(axis.title.x = element_blank())


VlnPlot(all_filtered_seurat, features = c("nGene"),
        pt.size = 0)+NoLegend()+
  theme(axis.title.x = element_blank())

VlnPlot(all_filtered_seurat, features = c("nUMI"),
        pt.size = 0)+NoLegend()+
  theme(axis.title.x = element_blank())





VlnPlot(all_filtered_seurat, features = c("nGene","mitoRatio"),
        pt.size = 0)

# save all filtered seurat object ----
save(all_filtered_seurat, file = "data/all_filtered_seurat.RData")
# Load all filtered seurat data
load("data/all_filtered_seurat.RData")

# Annotate cells using projection to Shalek lab dataset----
nasal_reference <- readRDS("Olympia/nasal_reference_REAL.rds")
load("data/all_filtered_seurat.RData")

view(nasal_reference@meta.data) 

# Mapping----
anchors <- FindTransferAnchors(reference = nasal_reference,
                               query = all_filtered_seurat,
                               normalization.method = "SCT",
                               reference.reduction = "pca",
                               dims = 1:30)



all_filtered_seurat <- MapQuery(anchorset = anchors,
                            reference = nasal_reference,
                            query = all_filtered_seurat,
                            refdata = list(celltype="Coarse_Cell_Annotations"),
                            reference.reduction = "pca",
                            reduction.model = "umap")

# Save filtered seurat----
save(all_filtered_seurat, file = "data/all_filtered_seurat.RData")
load("data/all_filtered_seurat.RData")


# SPLIT seurat object----
all_filtered_seurat[["RNA"]]<-split(all_filtered_seurat[["RNA"]],
                                f=all_filtered_seurat$sample)


all_filtered_seurat<-NormalizeData(all_filtered_seurat, 
                                   assay = "RNA",
                                   normalization.method = "LogNormalize",)
all_filtered_seurat<-FindVariableFeatures(all_filtered_seurat,
                                          selection.method = "vst",
                                          nfeatures = 2000)
all_filtered_seurat<-ScaleData(all_filtered_seurat,
                               vars.to.regress = c("mitoRatio","nGene","nUMI","Carriage_Status"))
all_filtered_seurat<-RunPCA(all_filtered_seurat, npcs = 50, ndims.print = 1:10)

ElbowPlot(all_filtered_seurat, ndims = 50)+
  labs(title="Elbow plot",x="Principle Component", y="Standard Deviation")

# ----
save(all_filtered_seurat, file = "data/all_filtered_seurat.RData")
load("data/all_filtered_seurat.RData")

DimPlot(all_filtered_seurat, 
        reduction = "pca"
        #, group.by = c("sex","Carriage_Status")
        )

all_filtered_seurat<-FindNeighbors(all_filtered_seurat, dims=1:20, reduction = "pca")
all_filtered_seurat<-FindClusters(all_filtered_seurat, resolution = 0.3)
all_filtered_seurat<-RunUMAP(all_filtered_seurat, dims = 1:20, reduction = "pca")
all_filtered_seurat<-RunTSNE(all_filtered_seurat, dims = 1:20, reduction = "pca")

DimPlot(all_filtered_seurat, reduction = "umap",
        #group.by = "sample",
        #split.by = "sample",
        label = TRUE)+
  theme_classic()+
  theme(axis.line = element_line())


DimPlot(all_filtered_seurat, reduction = "tsne",
        #split.by = c("sample"),
        #group.by = c("sample"),
        label = TRUE)+
  labs(x="tSNE-1", y="tSNE-2")+
  theme_classic()+
  theme(axis.line = element_line())

# -----
save(all_filtered_seurat, file = "data/all_filtered_seurat.RData")
load("data/all_filtered_seurat.RData")
# Harmony integration----
all_filtered_seurat<-IntegrateLayers(object=all_filtered_seurat,
                                 method=CCAIntegration,
                                 group.by = "HIV_Status",
                                 dims = 1:20,
                                 orig.reduction ="pca",
                                 new.reduction="CCA",
                                 verbose = FALSE)


all_filtered_seurat<-FindNeighbors(all_filtered_seurat, 
                                   reduction = "CCA",
                                   dims = 1:30)

all_filtered_seurat<-FindClusters(all_filtered_seurat, resolution = 0.3,
                              cluster.name = "CCA_clusters",
                              method = "igraph")

all_filtered_seurat<-RunUMAP(all_filtered_seurat, reduction = "CCA",
                         dims = 1:30,reduction.name = "umap.CCA")

all_filtered_seurat<-RunTSNE(all_filtered_seurat, reduction = "CCA",
                         dims = 1:30,reduction.name = "umap.CCA")


DimPlot(all_filtered_seurat,
        reduction = "umap.CCA",
        label = T,
        repel = T, 
        #pt.size = 1.2,
        #group.by = c("predicted.celltype"),
        #alpha = 0.9,
        #split.by=c("sample")
)+
  theme_bw()+NoLegend()


# Save all filtered seurat object----
save(all_filtered_seurat, file = "data/all_filtered_seurat.RData")
load("data/all_filtered_seurat.RData")

all_filtered_seurat_subset <- subset(all_filtered_seurat,
                         idents = c("0"),
                         invert = T)



all_filtered_seurat_subset<-IntegrateLayers(object=all_filtered_seurat_subset,
                                     method=HarmonyIntegration,
                                     group.by = "HIV_Status",
                                     dims = 1:20,
                                     orig.reduction ="pca",
                                     new.reduction="harmony",
                                     verbose = FALSE)


all_filtered_seurat_subset<-FindNeighbors(all_filtered_seurat_subset, 
                                          reduction = "harmony",
                                          dims = 1:20)

all_filtered_seurat_subset<-FindClusters(all_filtered_seurat_subset, resolution = 0.3,
                                  cluster.name = "harmony_clusters",
                                  method = "igraph",
                                  group.singletons=T)

all_filtered_seurat_subset<-RunUMAP(all_filtered_seurat_subset, reduction = "harmony",
                             dims = 1:20,reduction.name = "umap.harmony")

all_filtered_seurat_subset<-RunTSNE(all_filtered_seurat_subset, reduction = "harmony",
                             dims = 1:20,reduction.name = "tsne.harmony")


DimPlot(all_filtered_seurat_subset,
        reduction = "umap.harmony",
        label = T,
        repel = F, 
        #pt.size = 1.2,
        #group.by = c("predicted.celltype"),
        #alpha = 0.9,
        split.by=c("HIV_Status")
)+
  theme_bw()#+NoLegend()


# Save all_filtered seurat subset object ----
save(all_filtered_seurat_subset, file = "data/all_filtered_seurat_subset.RData")
load("data/all_filtered_seurat_subset.RData")


DimPlot(all_filtered_seurat_subset,
        reduction = "umap.harmony",
        label = T)+
  theme_bw()


FeaturePlot(all_filtered_seurat_subset,reduction = "umap.harmony",
            features = c(#"KRT15","SPINK4","FCGBP","KRT7","CXCL17","F3","MUC5AC","BPIFA1",
                         "FOXJ1","SCEL","MUC5AC","PTPRC","CSF3R","CD3D"),
            #split.by = c("sample")
            )

 
# Cell Annotation using singleR package ----
ref<-celldex::HumanPrimaryCellAtlasData()
view(ref@colData)

Immune_cells_results<-SingleR(test = as.SingleCellExperiment(immune_cells),
                              ref=ref, labels = ref$label.main)
immune_cells$singleR_labels<-Immune_cells_results$labels
immune_cells[[]]


DimPlot(immune_cells, 
        reduction = "umap.harmony",
        label = T,
        pt.size = 0.9,
        alpha = 0.9,
        #group.by = "singleR_labels",
        #split.by = "HIV_Status"
        )+
  theme_bw()



FeaturePlot(immune_cells,reduction = "umap.harmony",
            features = c("CD3D","CD8A","CD19","MZB1"#,
                         #"CD19","MZB1","NCAM1"#,"KLRB1","CD68","MS4A7",
              #"PTPRC"
              ),
            pt.size = 0.1,
            #split.by = c("HIV_Status")
)

meta <- immune_cells@meta.data
umap <- immune_cells@reductions$umap.harmony@cell.embeddings
soup_channel <- setClusters()
immune_cells_sc <- autoEstCont(immune_cells)


immune_cell_markers <- FindAllMarkers(
  immune_cells,
  assay = "RNA",
  logfc.threshold = 0.25,
  test.use = "bimod",
  slot = "data",
  min.pct = 0.1,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = TRUE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  min.cells.feature = 3,
  min.cells.group = 3,
  mean.fxn = NULL,
  fc.name = NULL,
  base = 2,
  return.thresh = 0.01,
  densify = FALSE)

write.csv(immune_cell_markers, file = "results/immune_cell_markers.csv")

immune_cells <- subset(immune_cells,
                         idents = c("4","6","9"),
                         invert = T)

immune_cells_labelled <- RenameIdents(object = immune_cells,
                             "0"="T cells",
                             "1"="T cells",
                             "2"="Neutrophils",
                             "3"="T cells",
                             "5"="B cells",
                             "7"="Macrophages",
                             "8"="Monocytes")



DimPlot(immune_cells_labelled,
        reduction = "umap.harmony",
        label = T,
        label.size = 5,
        #split.by = "HIV_Status"
        )+
  theme_bw()+NoLegend()


save(immune_cells_labelled, file = "data/immune_cells_labelled.RData")
load("data/immune_cells_labelled.RData")


# Rename idents of epithelial cell object clusters ----
all_filtered_seurat_subset <- RenameIdents(object = all_filtered_seurat_subset,
                                          "0" = "",
                                          "1" = "Monocytes/Macrophages",
                                          "2" = "",
                                          "3" = "",
                                          "4" = "Squamous cells",
                                          "5" = "Ciliated cells",
                                          "6" = "",
                                          "7" = "",
                                          "8" = "",
                                          "9" = " T cells",
                                          "10" = "Neutrophils",
                                          "11" = "Ciliated cells",
                                          "12" = "",
                                          "13" = "Ciliated cells",
                                          "14" = "B cells",
                                          "15" = "",
                                          "16" = "Ionocytes")



Idents(all_filtered_seurat)<-"predicted.celltype"
view(all_filtered_seurat@meta.data)




