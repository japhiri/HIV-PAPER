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




devtools::install_github("yingyonghui/CommPath")
library(CommPath)

remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)
remotes::install_github("mojaveazure/seurat-object", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)
remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE)
devtools::install_github("BlishLab/scriabin", ref = "main")
install.packages("devtools")
devtools::install_github("saeyslab/nichenetr")
devtools::install_github("saeyslab/multinichenetr")

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
  seurat_data <- Read10X(data.dir = paste0("Data/Single_Cell_Data/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.features = 30,
                                   project = file)
  assign(file, seurat_obj)
}

# -----CUH124----
CUH124 <- CUH124_raw_feature_bc_matrix

# Add number of genes per UMI for each cell to metadata
CUH124$log10GenesPerUMI <- log10(CUH124$nFeature_RNA) / log10(CUH124$nCount_RNA)

# Compute percent mito ratio 
CUH124$mitoRatio <- PercentageFeatureSet(object = CUH124, pattern = "^MT-")
CUH124$mitoRatio <- CUH124@meta.data$mitoRatio / 100

# Rename metadata columns
CUH124@meta.data <- CUH124@meta.data %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)


length(CUH124@meta.data$seq_folder)
view(CUH124@meta.data)


# Standard seurat workflow 
CUH124 <- NormalizeData(CUH124)
CUH124 <- FindVariableFeatures(CUH124)
CUH124 <- ScaleData(CUH124)
CUH124 <- RunPCA(CUH124)
CUH124 <- FindNeighbors(CUH124)
CUH124 <- FindClusters(CUH124)
CUH124 <- RunUMAP(CUH124,dims = 1:30)
DimPlot(CUH124,reduction = "umap")

# Featureplot
VlnPlot(CUH124,
        features = c("mitoRatio","nGene","nUMI"))


CUH124_filt <- subset(x = CUH124,
                 subset= (nUMI >= 100) & 
                   (nUMI <= 50000) &
                   (nGene >= 100) &
                   (mitoRatio < 0.25))
VlnPlot(CUH124_filt,
        features = c("mitoRatio","nGene","nUMI"))


# Standard seurat workflow 
CUH124_filt <- NormalizeData(CUH124_filt)
CUH124_filt <- FindVariableFeatures(CUH124_filt)
CUH124_filt <- ScaleData(CUH124_filt)
CUH124_filt <- RunPCA(CUH124_filt)
CUH124_filt <- FindNeighbors(CUH124_filt)
CUH124_filt <- FindClusters(CUH124_filt)
CUH124_filt <- RunUMAP(CUH124_filt,dims = 1:30)
DimPlot(CUH124_filt,reduction = "umap")

FeaturePlot(CUH124_filt,
            reduction = "umap",
            features = c("CD3D","CSF3R","CD19","FOXJ1"))


length(CUH124_filt@meta.data$nUMI)
# View metadata 
View(CUH124_filt@meta.data)

# -----CUG11X----
CUG11X <- CUG11X_raw_feature_bc_matrix

# Add number of genes per UMI for each cell to metadata
CUG11X$log10GenesPerUMI <- log10(CUG11X$nFeature_RNA) / log10(CUG11X$nCount_RNA)

# Compute percent mito ratio 
CUG11X$mitoRatio <- PercentageFeatureSet(object = CUG11X, pattern = "^MT-")
CUG11X$mitoRatio <- CUG11X@meta.data$mitoRatio / 100

# Rename metadata columns
CUG11X@meta.data <- CUG11X@meta.data %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)


length(CUG11X@meta.data$seq_folder)
view(CUG11X@meta.data)

# Standard seurat workflow 
CUG11X <- NormalizeData(CUG11X)
CUG11X <- FindVariableFeatures(CUG11X)
CUG11X <- ScaleData(CUG11X)
CUG11X <- RunPCA(CUG11X)
CUG11X <- FindNeighbors(CUG11X)
CUG11X <- FindClusters(CUG11X)
CUG11X <- RunUMAP(CUG11X,dims = 1:30)
DimPlot(CUG11X,reduction = "umap")

# Featureplot
VlnPlot(CUG11X,
        features = c("mitoRatio","nGene","nUMI"))


CUG11X_filt <- subset(x = CUG11X,
                      subset= (nUMI >= 100) & 
                        (nUMI <= 50000) &
                        (nGene >= 100) &
                        (mitoRatio < 0.20))


length(CUG11X_filt@meta.data$seq_folder)

CUG11X_filt <- NormalizeData(CUG11X_filt)
CUG11X_filt <- FindVariableFeatures(CUG11X_filt)
CUG11X_filt <- ScaleData(CUG11X_filt)
CUG11X_filt <- RunPCA(CUG11X_filt)
CUG11X_filt <- FindNeighbors(CUG11X_filt)
CUG11X_filt <- FindClusters(CUG11X_filt)
CUG11X_filt <- RunUMAP(CUG11X_filt,dims = 1:30)
DimPlot(CUG11X_filt,reduction = "umap")

VlnPlot(CUG11X_filt,
        features = c("mitoRatio","nGene","nUMI"))

FeaturePlot(CUG11X_filt,
            reduction = "umap",
            features = c("CD3D","CSF3R","CD19","FOXJ1"))

DimPlot(CUG11X_filt,reduction = "umap")

length(CUG11X_filt@meta.data$nUMI)
# View metadata 
View(CUG11X_filt@meta.data)

# -----CUF130----
CUF130 <- CUF130_raw_feature_bc_matrix

# Add number of genes per UMI for each cell to metadata
CUF130$log10GenesPerUMI <- log10(CUF130$nFeature_RNA) / log10(CUF130$nCount_RNA)

# Compute percent mito ratio 
CUF130$mitoRatio <- PercentageFeatureSet(object = CUF130, pattern = "^MT-")
CUF130$mitoRatio <- CUF130@meta.data$mitoRatio / 100

# Rename metadata columns
CUF130@meta.data <- CUF130@meta.data %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Standard seurat workflow 
CUF130 <- NormalizeData(CUF130)
CUF130 <- FindVariableFeatures(CUF130)
CUF130 <- ScaleData(CUF130)
CUF130 <- RunPCA(CUF130)
CUF130 <- FindNeighbors(CUF130)
CUF130 <- FindClusters(CUF130)
CUF130 <- RunUMAP(CUF130,dims = 1:30)
DimPlot(CUF130,reduction = "umap")

# Featureplot
VlnPlot(CUF130,
        features = c("mitoRatio","nGene","nUMI"))


CUF130_filt <- subset(x = CUF130,
                      subset= (nUMI >= 100) & 
                        (nUMI <= 50000) &
                        (nGene >= 150) &
                        (mitoRatio < 0.30))

CUF130_filt <- NormalizeData(CUF130_filt)
CUF130_filt <- FindVariableFeatures(CUF130_filt)
CUF130_filt <- ScaleData(CUF130_filt)
CUF130_filt <- RunPCA(CUF130_filt)
CUF130_filt <- FindNeighbors(CUF130_filt)
CUF130_filt <- FindClusters(CUF130_filt)
CUF130_filt <- RunUMAP(CUF130_filt,dims = 1:30)
DimPlot(CUF130_filt,reduction = "umap")

VlnPlot(CUF130_filt,
        features = c("mitoRatio","nGene","nUMI"))

FeaturePlot(CUF130_filt,
            reduction = "umap",
            features = c("CD3D","CSF3R","CD19","FOXJ1"))


length(CUF130_filt@meta.data$nUMI)
# View metadata 
View(CUF130_filt@meta.data)

# -----CUF131----
CUF131 <- CUF131_raw_feature_bc_matrix

# Add number of genes per UMI for each cell to metadata
CUF131$log10GenesPerUMI <- log10(CUF131$nFeature_RNA) / log10(CUF131$nCount_RNA)

# Compute percent mito ratio 
CUF131$mitoRatio <- PercentageFeatureSet(object = CUF131, pattern = "^MT-")
CUF131$mitoRatio <- CUF131@meta.data$mitoRatio / 100

# Rename metadata columns
CUF131@meta.data <- CUF131@meta.data %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Standard seurat workflow 
CUF131 <- NormalizeData(CUF131)
CUF131 <- FindVariableFeatures(CUF131)
CUF131 <- ScaleData(CUF131)
CUF131 <- RunPCA(CUF131)
CUF131 <- FindNeighbors(CUF131)
CUF131 <- FindClusters(CUF131)
CUF131 <- RunUMAP(CUF131,dims = 1:30)
DimPlot(CUF131,reduction = "umap")

# Featureplot
VlnPlot(CUF131,
        features = c("mitoRatio","nGene","nUMI"))


CUF131_filt <- subset(x = CUF131,
                      subset= (nUMI >= 100) & 
                        (nUMI <= 50000) &
                        (nGene >= 100) &
                        (mitoRatio < 0.30))

CUF131_filt <- NormalizeData(CUF131_filt)
CUF131_filt <- FindVariableFeatures(CUF131_filt)
CUF131_filt <- ScaleData(CUF131_filt)
CUF131_filt <- RunPCA(CUF131_filt)
CUF131_filt <- FindNeighbors(CUF131_filt)
CUF131_filt <- FindClusters(CUF131_filt)
CUF131_filt <- RunUMAP(CUF131_filt,dims = 1:30)
DimPlot(CUF131_filt,reduction = "umap")

VlnPlot(CUF131_filt,
        features = c("mitoRatio","nGene","nUMI"))

FeaturePlot(CUF131_filt,
            reduction = "umap",
            features = c("CD3D","CSF3R","CD19","FOXJ1"))

DimPlot(CUF131_filt,reduction = "umap")

length(CUF131_filt@meta.data$nUMI)
# View metadata 
View(CUF131_filt@meta.data)

# -----CUF134----
CUF134 <- CUF134_raw_feature_bc_matrix

# Add number of genes per UMI for each cell to metadata
CUF134$log10GenesPerUMI <- log10(CUF134$nFeature_RNA) / log10(CUF134$nCount_RNA)

# Compute percent mito ratio 
CUF134$mitoRatio <- PercentageFeatureSet(object = CUF134, pattern = "^MT-")
CUF134$mitoRatio <- CUF134@meta.data$mitoRatio / 100

# Rename metadata columns
CUF134@meta.data <- CUF134@meta.data %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Standard seurat workflow 
CUF134 <- NormalizeData(CUF134)
CUF134 <- FindVariableFeatures(CUF134)
CUF134 <- ScaleData(CUF134)
CUF134 <- RunPCA(CUF134)
CUF134 <- FindNeighbors(CUF134)
CUF134 <- FindClusters(CUF134)
CUF134 <- RunUMAP(CUF134,dims = 1:30)
DimPlot(CUF134,reduction = "umap")

# Featureplot
VlnPlot(CUF134,
        features = c("mitoRatio","nGene","nUMI"))


CUF134_filt <- subset(x = CUF134,
                      subset= (nUMI >= 100) & 
                        (nUMI <= 50000) &
                        (nGene >= 100) &
                        (mitoRatio < 0.25))

CUF134_filt <- NormalizeData(CUF134_filt)
CUF134_filt <- FindVariableFeatures(CUF134_filt)
CUF134_filt <- ScaleData(CUF134_filt)
CUF134_filt <- RunPCA(CUF134_filt)
CUF134_filt <- FindNeighbors(CUF134_filt)
CUF134_filt <- FindClusters(CUF134_filt)
CUF134_filt <- RunUMAP(CUF134_filt,dims = 1:30)
DimPlot(CUF134_filt,reduction = "umap")

VlnPlot(CUF134_filt,
        features = c("mitoRatio","nGene","nUMI"))

FeaturePlot(CUF134_filt,
            reduction = "umap",
            features = c("CD3D","CSF3R","CD19","KRT7"))

DimPlot(CUF134_filt,reduction = "umap")

length(CUF134_filt@meta.data$nUMI)
# View metadata 
View(CUF134_filt@meta.data)

# -----CUF135----
CUF135 <- CUF135_raw_feature_bc_matrix

# Add number of genes per UMI for each cell to metadata
CUF135$log10GenesPerUMI <- log10(CUF135$nFeature_RNA) / log10(CUF135$nCount_RNA)

# Compute percent mito ratio 
CUF135$mitoRatio <- PercentageFeatureSet(object = CUF135, pattern = "^MT-")
CUF135$mitoRatio <- CUF135@meta.data$mitoRatio / 100

# Rename metadata columns
CUF135@meta.data <- CUF135@meta.data %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Standard seurat workflow 
CUF135 <- NormalizeData(CUF135)
CUF135 <- FindVariableFeatures(CUF135)
CUF135 <- ScaleData(CUF135)
CUF135 <- RunPCA(CUF135)
CUF135 <- FindNeighbors(CUF135)
CUF135 <- FindClusters(CUF135)
CUF135 <- RunUMAP(CUF135,dims = 1:30)
DimPlot(CUF135,reduction = "umap")

# Featureplot
VlnPlot(CUF135,
        features = c("mitoRatio","nGene","nUMI"))


CUF135_filt <- subset(x = CUF135,
                      subset= (nUMI >= 100) & 
                        (nUMI <= 50000) &
                        (nGene >= 100) &
                        (mitoRatio < 0.25))

CUF135_filt <- NormalizeData(CUF135_filt)
CUF135_filt <- FindVariableFeatures(CUF135_filt)
CUF135_filt <- ScaleData(CUF135_filt)
CUF135_filt <- RunPCA(CUF135_filt)
CUF135_filt <- FindNeighbors(CUF135_filt)
CUF135_filt <- FindClusters(CUF135_filt)
CUF135_filt <- RunUMAP(CUF135_filt,dims = 1:30)
DimPlot(CUF135_filt,reduction = "umap")

VlnPlot(CUF135_filt,
        features = c("mitoRatio","nGene","nUMI"))

FeaturePlot(CUF135_filt,
            reduction = "umap",
            features = c("CD3D","CSF3R","CD19","KRT7"))

DimPlot(CUF135_filt,reduction = "umap")

length(CUF135_filt@meta.data$nUMI)
# View metadata 
View(CUF135_filt@meta.data)

# -----CUF136----
CUF136 <- CUF136_raw_feature_bc_matrix

# Add number of genes per UMI for each cell to metadata
CUF136$log10GenesPerUMI <- log10(CUF136$nFeature_RNA) / log10(CUF136$nCount_RNA)

# Compute percent mito ratio 
CUF136$mitoRatio <- PercentageFeatureSet(object = CUF136, pattern = "^MT-")
CUF136$mitoRatio <- CUF136@meta.data$mitoRatio / 100

# Rename metadata columns
CUF136@meta.data <- CUF136@meta.data %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Standard seurat workflow 
CUF136 <- NormalizeData(CUF136)
CUF136 <- FindVariableFeatures(CUF136)
CUF136 <- ScaleData(CUF136)
CUF136 <- RunPCA(CUF136)
CUF136 <- FindNeighbors(CUF136)
CUF136 <- FindClusters(CUF136)
CUF136 <- RunUMAP(CUF136,dims = 1:30)
DimPlot(CUF136,reduction = "umap")

# Featureplot
VlnPlot(CUF136,
        features = c("mitoRatio","nGene","nUMI"))


CUF136_filt <- subset(x = CUF136,
                      subset= (nUMI >= 100) & 
                        (nUMI <= 50000) &
                        (nGene >= 100) &
                        (mitoRatio < 0.30))

CUF136_filt <- NormalizeData(CUF136_filt)
CUF136_filt <- FindVariableFeatures(CUF136_filt)
CUF136_filt <- ScaleData(CUF136_filt)
CUF136_filt <- RunPCA(CUF136_filt)
CUF136_filt <- FindNeighbors(CUF136_filt)
CUF136_filt <- FindClusters(CUF136_filt)
CUF136_filt <- RunUMAP(CUF136_filt,dims = 1:30)
DimPlot(CUF136_filt,reduction = "umap")

VlnPlot(CUF136_filt,
        features = c("mitoRatio","nGene","nUMI"))

FeaturePlot(CUF136_filt,
            reduction = "umap",
            features = c("CD3D","CSF3R","CD19","KRT7"))

DimPlot(CUF136_filt,reduction = "umap")

length(CUF136_filt@meta.data$nUMI)
# View metadata 
View(CUF136_filt@meta.data)



# -----CUF137----
CUF137 <- CUF137_raw_feature_bc_matrix

# Add number of genes per UMI for each cell to metadata
CUF137$log10GenesPerUMI <- log10(CUF137$nFeature_RNA) / log10(CUF137$nCount_RNA)

# Compute percent mito ratio 
CUF137$mitoRatio <- PercentageFeatureSet(object = CUF137, pattern = "^MT-")
CUF137$mitoRatio <- CUF137@meta.data$mitoRatio / 100

# Rename metadata columns
CUF137@meta.data <- CUF137@meta.data %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Standard seurat workflow 
CUF137 <- NormalizeData(CUF137)
CUF137 <- FindVariableFeatures(CUF137)
CUF137 <- ScaleData(CUF137)
CUF137 <- RunPCA(CUF137)
CUF137 <- FindNeighbors(CUF137)
CUF137 <- FindClusters(CUF137)
CUF137 <- RunUMAP(CUF137,dims = 1:30)
DimPlot(CUF137,reduction = "umap")

# Featureplot
VlnPlot(CUF137,
        features = c("mitoRatio","nGene","nUMI"))


CUF137_filt <- subset(x = CUF137,
                      subset= (nUMI >= 100) & 
                        (nUMI <= 50000) &
                        (nGene >= 100) &
                        (mitoRatio < 0.20))

CUF137_filt <- NormalizeData(CUF137_filt)
CUF137_filt <- FindVariableFeatures(CUF137_filt)
CUF137_filt <- ScaleData(CUF137_filt)
CUF137_filt <- RunPCA(CUF137_filt)
CUF137_filt <- FindNeighbors(CUF137_filt)
CUF137_filt <- FindClusters(CUF137_filt)
CUF137_filt <- RunUMAP(CUF137_filt,dims = 1:30)
DimPlot(CUF137_filt,reduction = "umap")

VlnPlot(CUF137_filt,
        features = c("mitoRatio","nGene","nUMI"))

FeaturePlot(CUF137_filt,
            reduction = "umap",
            features = c("CD3D","CSF3R","CD19","KRT7"))

DimPlot(CUF137_filt,reduction = "umap")

length(CUF137_filt@meta.data$nUMI)
# View metadata 
View(CUF137_filt@meta.data)



# -----CUF13J----
CUF13J <- CUF13J_raw_feature_bc_matrix

# Add number of genes per UMI for each cell to metadata
CUF13J$log10GenesPerUMI <- log10(CUF13J$nFeature_RNA) / log10(CUF13J$nCount_RNA)

# Compute percent mito ratio 
CUF13J$mitoRatio <- PercentageFeatureSet(object = CUF13J, pattern = "^MT-")
CUF13J$mitoRatio <- CUF13J@meta.data$mitoRatio / 100

# Rename metadata columns
CUF13J@meta.data <- CUF13J@meta.data %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Standard seurat workflow 
CUF13J <- NormalizeData(CUF13J)
CUF13J <- FindVariableFeatures(CUF13J)
CUF13J <- ScaleData(CUF13J)
CUF13J <- RunPCA(CUF13J)
CUF13J <- FindNeighbors(CUF13J)
CUF13J <- FindClusters(CUF13J)
CUF13J <- RunUMAP(CUF13J,dims = 1:30)
DimPlot(CUF13J,reduction = "umap")

# Featureplot
VlnPlot(CUF13J,
        features = c("mitoRatio","nGene","nUMI"))


CUF13J_filt <- subset(x = CUF13J,
                      subset= (nUMI >= 100) & 
                        (nUMI <= 50000) &
                        (nGene >= 100) &
                        (mitoRatio < 0.30))

CUF13J_filt <- NormalizeData(CUF13J_filt)
CUF13J_filt <- FindVariableFeatures(CUF13J_filt)
CUF13J_filt <- ScaleData(CUF13J_filt)
CUF13J_filt <- RunPCA(CUF13J_filt)
CUF13J_filt <- FindNeighbors(CUF13J_filt)
CUF13J_filt <- FindClusters(CUF13J_filt)
CUF13J_filt <- RunUMAP(CUF13J_filt,dims = 1:30)
DimPlot(CUF13J_filt,reduction = "umap")

VlnPlot(CUF13J_filt,
        features = c("mitoRatio","nGene","nUMI"))

FeaturePlot(CUF13J_filt,
            reduction = "umap",
            features = c("CD3D","CSF3R","CD19","KRT7"))

DimPlot(CUF13J_filt,reduction = "umap")

length(CUF13J_filt@meta.data$nUMI)
# View metadata 
View(CUF13J_filt@meta.data)



# -----CUF13K----
CUF13K <- CUF13K_raw_feature_bc_matrix

# Add number of genes per UMI for each cell to metadata
CUF13K$log10GenesPerUMI <- log10(CUF13K$nFeature_RNA) / log10(CUF13K$nCount_RNA)

# Compute percent mito ratio 
CUF13K$mitoRatio <- PercentageFeatureSet(object = CUF13K, pattern = "^MT-")
CUF13K$mitoRatio <- CUF13K@meta.data$mitoRatio / 100

# Rename metadata columns
CUF13K@meta.data <- CUF13K@meta.data %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Standard seurat workflow 
CUF13K <- NormalizeData(CUF13K)
CUF13K <- FindVariableFeatures(CUF13K)
CUF13K <- ScaleData(CUF13K)
CUF13K <- RunPCA(CUF13K)
CUF13K <- FindNeighbors(CUF13K)
CUF13K <- FindClusters(CUF13K)
CUF13K <- RunUMAP(CUF13K,dims = 1:30)
DimPlot(CUF13K,reduction = "umap")

# Featureplot
VlnPlot(CUF13K,
        features = c("mitoRatio","nGene","nUMI"))


CUF13K_filt <- subset(x = CUF13K,
                      subset= (nUMI >= 100) & 
                        (nUMI <= 50000) &
                        (nGene >= 100) &
                        (mitoRatio < 0.30))

CUF13K_filt <- NormalizeData(CUF13K_filt)
CUF13K_filt <- FindVariableFeatures(CUF13K_filt)
CUF13K_filt <- ScaleData(CUF13K_filt)
CUF13K_filt <- RunPCA(CUF13K_filt)
CUF13K_filt <- FindNeighbors(CUF13K_filt)
CUF13K_filt <- FindClusters(CUF13K_filt)
CUF13K_filt <- RunUMAP(CUF13K_filt,dims = 1:30)
DimPlot(CUF13K_filt,reduction = "umap")

VlnPlot(CUF13K_filt,
        features = c("mitoRatio","nGene","nUMI"))

FeaturePlot(CUF13K_filt,
            reduction = "umap",
            features = c("CD3D","CSF3R","CD19","KRT7"))

DimPlot(CUF13K_filt,reduction = "umap")

length(CUF13K_filt@meta.data$nUMI)
# View metadata 
View(CUF13K_filt@meta.data)


# -----CUF12Y----
CUF12Y <- CUF12Y_raw_feature_bc_matrix

# Add number of genes per UMI for each cell to metadata
CUF12Y$log10GenesPerUMI <- log10(CUF12Y$nFeature_RNA) / log10(CUF12Y$nCount_RNA)

# Compute percent mito ratio 
CUF12Y$mitoRatio <- PercentageFeatureSet(object = CUF12Y, pattern = "^MT-")
CUF12Y$mitoRatio <- CUF12Y@meta.data$mitoRatio / 100

# Rename metadata columns
CUF12Y@meta.data <- CUF12Y@meta.data %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Standard seurat workflow 
CUF12Y <- NormalizeData(CUF12Y)
CUF12Y <- FindVariableFeatures(CUF12Y)
CUF12Y <- ScaleData(CUF12Y)
CUF12Y <- RunPCA(CUF12Y)
CUF12Y <- FindNeighbors(CUF12Y)
CUF12Y <- FindClusters(CUF12Y)
CUF12Y <- RunUMAP(CUF12Y,dims = 1:30)
DimPlot(CUF12Y,reduction = "umap")

# Featureplot
VlnPlot(CUF12Y,
        features = c("mitoRatio","nGene","nUMI"))


CUF12Y_filt <- subset(x = CUF12Y,
                      subset= (nUMI >= 100) & 
                        (nUMI <= 50000) &
                        (nGene >= 100) &
                        (mitoRatio < 0.30))

CUF12Y_filt <- NormalizeData(CUF12Y_filt)
CUF12Y_filt <- FindVariableFeatures(CUF12Y_filt)
CUF12Y_filt <- ScaleData(CUF12Y_filt)
CUF12Y_filt <- RunPCA(CUF12Y_filt)
CUF12Y_filt <- FindNeighbors(CUF12Y_filt)
CUF12Y_filt <- FindClusters(CUF12Y_filt)
CUF12Y_filt <- RunUMAP(CUF12Y_filt,dims = 1:30)
DimPlot(CUF12Y_filt,reduction = "umap")

VlnPlot(CUF12Y_filt,
        features = c("mitoRatio","nGene","nUMI"))

FeaturePlot(CUF12Y_filt,
            reduction = "umap",
            features = c("CD3D","CSF3R","CD19","KRT7"))

DimPlot(CUF12Y_filt,reduction = "umap")

length(CUF12Y_filt@meta.data$nUMI)
# View metadata 
View(CUF12Y_filt@meta.data)





# -----CUF13I----
CUF13I <- CUF13I_raw_feature_bc_matrix

# Add number of genes per UMI for each cell to metadata
CUF13I$log10GenesPerUMI <- log10(CUF13I$nFeature_RNA) / log10(CUF13I$nCount_RNA)

# Compute percent mito ratio 
CUF13I$mitoRatio <- PercentageFeatureSet(object = CUF13I, pattern = "^MT-")
CUF13I$mitoRatio <- CUF13I@meta.data$mitoRatio / 100

# Rename metadata columns
CUF13I@meta.data <- CUF13I@meta.data %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Standard seurat workflow 
CUF13I <- NormalizeData(CUF13I)
CUF13I <- FindVariableFeatures(CUF13I)
CUF13I <- ScaleData(CUF13I)
CUF13I <- RunPCA(CUF13I)
CUF13I <- FindNeighbors(CUF13I)
CUF13I <- FindClusters(CUF13I)
CUF13I <- RunUMAP(CUF13I,dims = 1:30)
DimPlot(CUF13I,reduction = "umap")

# Featureplot
VlnPlot(CUF13I,
        features = c("mitoRatio","nGene","nUMI"))


CUF13I_filt <- subset(x = CUF13I,
                      subset= (nUMI >= 100) & 
                        (nUMI <= 50000) &
                        (nGene >= 100) &
                        (mitoRatio < 0.30))

CUF13I_filt <- NormalizeData(CUF13I_filt)
CUF13I_filt <- FindVariableFeatures(CUF13I_filt)
CUF13I_filt <- ScaleData(CUF13I_filt)
CUF13I_filt <- RunPCA(CUF13I_filt)
CUF13I_filt <- FindNeighbors(CUF13I_filt)
CUF13I_filt <- FindClusters(CUF13I_filt)
CUF13I_filt <- RunUMAP(CUF13I_filt,dims = 1:30)
DimPlot(CUF13I_filt,reduction = "umap")

VlnPlot(CUF13K_filt,
        features = c("mitoRatio","nGene","nUMI"))

FeaturePlot(CUF13K_filt,
            reduction = "umap",
            features = c("CD3D","CSF3R","CD19","KRT7"))

DimPlot(CUF13I_filt,reduction = "umap")

length(CUF13I_filt@meta.data$nUMI)
# View metadata 
View(CUF13I_filt@meta.data)





# Create a merged Seurat object ----
all_merged <- merge(x=CUF131_filt,
                    y=c(CUH124_filt,
                        CUG11X_filt,
                        CUF130_filt,
                        CUF134_filt,
                        CUF13I_filt,
                        CUF135_filt,
                        CUF136_filt,
                        CUF13J_filt,
                        CUF13K_filt,
                        #CUF12Y_filt,
                        CUF137_filt),
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
                                    #"CUF12Y",
                                    "CUF137"))


# save all_merged object as .RData file -----
save(all_merged, file = "Data/Single_Cell_Data/all_merged.RData")
load("Data/Single_Cell_Data/all_merged.RData")
view(all_merged@meta.data)
length(all_merged@meta.data$seq_folder)

all_merged[["RNA"]]<-split(all_merged[["RNA"]],f=all_merged$seq_folder)



# Standard seurat workflow 
all_merged <- NormalizeData(all_merged,assay = "RNA",normalization.method = "LogNormalize")
all_merged <- FindVariableFeatures(all_merged,selection.method = "vst",nfeatures = 2000)
all_merged <- ScaleData(all_merged,vars.to.regress = c("mitoRatio","nGene","nUMI"))
all_merged <- RunPCA(all_merged,npcs = 50, ndims.print = 1:10)

# Harmony integration----
all_merged<-IntegrateLayers(object=all_merged,
                                     method=HarmonyIntegration,
                                     group.by = "seq_folder",
                                     dims = 1:20,
                                     orig.reduction ="pca",
                                     new.reduction="harmony",
                                     verbose = FALSE)


all_merged<-FindNeighbors(all_merged, 
                                   reduction = "harmony",
                                   dims = 1:30)

all_merged<-FindClusters(all_merged, resolution = 0.3,
                                  cluster.name = "harmony_clusters",
                                  method = "igraph")

all_merged<-RunUMAP(all_merged, reduction = "harmony",
                             dims = 1:30,reduction.name = "umap.harmony")

all_merged<-RunTSNE(all_merged, reduction = "harmony",
                             dims = 1:30,reduction.name = "tsne.harmony")


DimPlot(all_merged,
        reduction = "umap.harmony",
        label = T,
        repel = T, 
        #pt.size = 1.2,
        #group.by = c("predicted.celltype"),
        #alpha = 0.9,
        #split.by=c("seq_folder")
)+
  theme_bw()+NoLegend()





# creating a metadata file----
metadata <- all_merged@meta.data
# Adding metadata to seurat object----
# ----Sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$seq_folder, "^CUH124"))] <- "CUH124"
metadata$sample[which(str_detect(metadata$seq_folder, "^CUG11X"))] <- "CUG11X"
metadata$sample[which(str_detect(metadata$seq_folder, "^CUF131"))] <- "CUF131"
metadata$sample[which(str_detect(metadata$seq_folder, "^CUF130"))] <- "CUF130"
metadata$sample[which(str_detect(metadata$seq_folder, "^CUF134"))] <- "CUF134"
metadata$sample[which(str_detect(metadata$seq_folder, "^CUF13I"))] <- "CUF13I"
metadata$sample[which(str_detect(metadata$seq_folder, "^CUF135"))] <- "CUF135"
metadata$sample[which(str_detect(metadata$seq_folder, "^CUF136"))] <- "CUF136"
metadata$sample[which(str_detect(metadata$seq_folder, "^CUF13J"))] <- "CUF13J"
metadata$sample[which(str_detect(metadata$seq_folder, "^CUF13K"))] <- "CUF13K"
metadata$sample[which(str_detect(metadata$seq_folder, "^CUF137"))] <- "CUF137"
#metadata$sample[which(str_detect(metadata$seq_folder, "^CUF12Y"))] <- "CUF12Y"

# ---- Create HIV status column
metadata$HIV_Status <- NA
metadata$HIV_Status[which(str_detect(metadata$seq_folder, "^CUH124"))] <- "HIV+ ART>1 Year"
metadata$HIV_Status[which(str_detect(metadata$seq_folder, "^CUG11X"))] <- "HIV+ ART>1 Year"
metadata$HIV_Status[which(str_detect(metadata$seq_folder, "^CUF131"))] <- "HIV-"
metadata$HIV_Status[which(str_detect(metadata$seq_folder, "^CUF130"))] <- "HIV+ ART>1 Year"
metadata$HIV_Status[which(str_detect(metadata$seq_folder, "^CUF134"))] <- "HIV+ ART<3 Months"
metadata$HIV_Status[which(str_detect(metadata$seq_folder, "^CUF13I"))] <- "HIV+ ART>1 Year"
metadata$HIV_Status[which(str_detect(metadata$seq_folder, "^CUF135"))] <- "HIV+ ART>1 Year"
metadata$HIV_Status[which(str_detect(metadata$seq_folder, "^CUF136"))] <- "HIV-"
metadata$HIV_Status[which(str_detect(metadata$seq_folder, "^CUF13J"))] <- "HIV+ ART<3 Months"
metadata$HIV_Status[which(str_detect(metadata$seq_folder, "^CUF13K"))] <- "HIV+ ART<3 Months"
metadata$HIV_Status[which(str_detect(metadata$seq_folder, "^CUF137"))] <- "HIV-"
#metadata$HIV_Status[which(str_detect(metadata$cells, "^CUF12Y"))] <- "HIV+ ART>1 Year"


# ---- Create Carriage status column
metadata$Carriage_Status <- NA
metadata$Carriage_Status[which(str_detect(metadata$seq_folder, "^CUH124"))] <- "SPN+"
metadata$Carriage_Status[which(str_detect(metadata$seq_folder, "^CUG11X"))] <- "SPN-"
metadata$Carriage_Status[which(str_detect(metadata$seq_folder, "^CUF131"))] <- "SPN+"
metadata$Carriage_Status[which(str_detect(metadata$seq_folder, "^CUF130"))] <- "SPN-"
metadata$Carriage_Status[which(str_detect(metadata$seq_folder, "^CUF134"))] <- "SPN+"
metadata$Carriage_Status[which(str_detect(metadata$seq_folder, "^CUF13I"))] <- "SPN+"
metadata$Carriage_Status[which(str_detect(metadata$seq_folder, "^CUF135"))] <- "SPN+"
metadata$Carriage_Status[which(str_detect(metadata$seq_folder, "^CUF136"))] <- "SPN+"
metadata$Carriage_Status[which(str_detect(metadata$seq_folder, "^CUF13J"))] <- "SPN+"
metadata$Carriage_Status[which(str_detect(metadata$seq_folder, "^CUF13K"))] <- "SPN-"
metadata$Carriage_Status[which(str_detect(metadata$seq_folder, "^CUF137"))] <- "SPN+"
#metadata$Carriage_Status[which(str_detect(metadata$seq_folder, "^CUF12Y"))] <- "SPN+"

# ---- Create Carriage density column
metadata$Carriage_density <- NA
metadata$Carriage_density[which(str_detect(metadata$seq_folder, "^CUH124"))] <- 5862500
metadata$Carriage_density[which(str_detect(metadata$seq_folder, "^CUG11X"))] <- NA
metadata$Carriage_density[which(str_detect(metadata$seq_folder, "^CUF131"))] <- 26800000
metadata$Carriage_density[which(str_detect(metadata$seq_folder, "^CUF130"))] <- NA
metadata$Carriage_density[which(str_detect(metadata$seq_folder, "^CUF134"))] <- 686750
metadata$Carriage_density[which(str_detect(metadata$seq_folder, "^CUF13I"))] <- 670
metadata$Carriage_density[which(str_detect(metadata$seq_folder, "^CUF135"))] <- 33500
metadata$Carriage_density[which(str_detect(metadata$seq_folder, "^CUF136"))] <- 67
metadata$Carriage_density[which(str_detect(metadata$seq_folder, "^CUF13J"))] <- 26800
metadata$Carriage_density[which(str_detect(metadata$seq_folder, "^CUF13K"))] <- NA
metadata$Carriage_density[which(str_detect(metadata$seq_folder, "^CUF137"))] <- 10720
#metadata$Carriage_density[which(str_detect(metadata$seq_folder, "^CUF12Y"))] <- 2345



# ---- Create Carriage Serotype column
metadata$serotype <- NA
metadata$serotype[which(str_detect(metadata$seq_folder, "^CUH124"))] <- "NVT"
metadata$serotype[which(str_detect(metadata$seq_folder, "^CUG11X"))] <- NA
metadata$serotype[which(str_detect(metadata$seq_folder, "^CUF131"))] <- "3"
metadata$serotype[which(str_detect(metadata$seq_folder, "^CUF130"))] <- NA
metadata$serotype[which(str_detect(metadata$seq_folder, "^CUF134"))] <- "9"
metadata$serotype[which(str_detect(metadata$seq_folder, "^CUF13I"))] <- "1"
metadata$serotype[which(str_detect(metadata$seq_folder, "^CUF135"))] <- "10"
metadata$serotype[which(str_detect(metadata$seq_folder, "^CUF136"))] <- "NVT"
metadata$serotype[which(str_detect(metadata$seq_folder, "^CUF13J"))] <- "10"
metadata$serotype[which(str_detect(metadata$seq_folder, "^CUF13K"))] <- NA
metadata$serotype[which(str_detect(metadata$seq_folder, "^CUF137"))] <- "NVT"
#metadata$serotype[which(str_detect(metadata$seq_folder, "^CUF12Y"))] <- "NVT"

# ---- Create AGE column
metadata$age <- NA
metadata$age[which(str_detect(metadata$seq_folder, "^CUH124"))] <- 35
metadata$age[which(str_detect(metadata$seq_folder, "^CUG11X"))] <- 21
metadata$age[which(str_detect(metadata$seq_folder, "^CUF131"))] <- 28
metadata$age[which(str_detect(metadata$seq_folder, "^CUF130"))] <- 34
metadata$age[which(str_detect(metadata$seq_folder, "^CUF134"))] <- 26
metadata$age[which(str_detect(metadata$seq_folder, "^CUF13I"))] <- 27
metadata$age[which(str_detect(metadata$seq_folder, "^CUF135"))] <- 34
metadata$age[which(str_detect(metadata$seq_folder, "^CUF136"))] <- 36
metadata$age[which(str_detect(metadata$seq_folder, "^CUF13J"))] <- 24
metadata$age[which(str_detect(metadata$seq_folder, "^CUF13K"))] <- 27
metadata$age[which(str_detect(metadata$seq_folder, "^CUF137"))] <- 28
#metadata$age[which(str_detect(metadata$seq_folder, "^CUF12Y"))] <- 25



# ---- Create SEX column
metadata$sex <- NA
metadata$sex[which(str_detect(metadata$seq_folder, "^CUH124"))] <- "Female"
metadata$sex[which(str_detect(metadata$seq_folder, "^CUG11X"))] <- "Female"
metadata$sex[which(str_detect(metadata$seq_folder, "^CUF131"))] <- "Female"
metadata$sex[which(str_detect(metadata$seq_folder, "^CUF130"))] <- "Female"
metadata$sex[which(str_detect(metadata$seq_folder, "^CUF134"))] <- "Male"
metadata$sex[which(str_detect(metadata$seq_folder, "^CUF13I"))] <- "Male"
metadata$sex[which(str_detect(metadata$seq_folder, "^CUF135"))] <- "Male"
metadata$sex[which(str_detect(metadata$seq_folder, "^CUF136"))] <- "Female"
metadata$sex[which(str_detect(metadata$seq_folder, "^CUF13J"))] <- "Female"
metadata$sex[which(str_detect(metadata$seq_folder, "^CUF13K"))] <- "Female"
metadata$sex[which(str_detect(metadata$seq_folder, "^CUF137"))] <- "Male"
#metadata$sex[which(str_detect(metadata$seq_folder, "^CUF12Y"))] <- "Male"


# ----Create viral load column
metadata$Viral_Load <- NA
metadata$Viral_Load[which(str_detect(metadata$seq_folder, "^CUH124"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$seq_folder, "^CUG11X"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$seq_folder, "^CUF131"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$seq_folder, "^CUF130"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$seq_folder, "^CUF134"))] <- 312500
metadata$Viral_Load[which(str_detect(metadata$seq_folder, "^CUF13I"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$seq_folder, "^CUF135"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$seq_folder, "^CUF136"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$seq_folder, "^CUF13J"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$seq_folder, "^CUF13K"))] <- "Not Detected"
metadata$Viral_Load[which(str_detect(metadata$seq_folder, "^CUF137"))] <- "Not Detected"
#metadata$Viral_Load[which(str_detect(metadata$cells, "^CUF12Y"))] <- "Not Detected"


# ----Create absolute CD4 count column
metadata$abs_CD4_count <- NA
metadata$abs_CD4_count[which(str_detect(metadata$seq_folder, "^CUH124"))] <- 805
metadata$abs_CD4_count[which(str_detect(metadata$seq_folder, "^CUG11X"))] <- 488
metadata$abs_CD4_count[which(str_detect(metadata$seq_folder, "^CUF131"))] <- 381
metadata$abs_CD4_count[which(str_detect(metadata$seq_folder, "^CUF130"))] <- 535
metadata$abs_CD4_count[which(str_detect(metadata$seq_folder, "^CUF134"))] <- 14
metadata$abs_CD4_count[which(str_detect(metadata$seq_folder, "^CUF13I"))] <- 452
metadata$abs_CD4_count[which(str_detect(metadata$seq_folder, "^CUF135"))] <- 274
metadata$abs_CD4_count[which(str_detect(metadata$seq_folder, "^CUF136"))] <- 1179
metadata$abs_CD4_count[which(str_detect(metadata$seq_folder, "^CUF13J"))] <- 689
metadata$abs_CD4_count[which(str_detect(metadata$seq_folder, "^CUF13K"))] <- 907
metadata$abs_CD4_count[which(str_detect(metadata$seq_folder, "^CUF137"))] <- 431
#metadata$abs_CD4_count[which(str_detect(metadata$seq_folder, "^CUF12Y"))] <- 311

# Add metadata back to Seurat object ----
all_merged@meta.data <- metadata


View(all_merged@meta.data)

# Dimplot----
DimPlot(all_merged,
        reduction = "umap.harmony",
        split.by = "HIV_Status",
        label = T)+
  theme_bw()

# Save all-merged seurat with metadata ----
save(all_merged,file = "Data/Single_Cell_Data/all_merged.RData")
load("Data/Single_Cell_Data//all_merged.RData")

# Remove clusters that only appear in 1 or 2 cells ----
all_merged_subset <- subset(all_merged,idents = c("0","1"), invert = T)

DimPlot(all_merged_subset,
        reduction = "umap.harmony",
        #split.by = "HIV_Status",
        label = T)+
  theme_bw()
# Harmony integration of subsetted samples----
all_merged_subset<-IntegrateLayers(object=all_merged_subset,
                            method=HarmonyIntegration,
                            group.by = "sample",
                            dims = 1:20,
                            orig.reduction ="pca",
                            new.reduction="harmony",
                            verbose = FALSE)


all_merged_subset<-FindNeighbors(all_merged_subset, 
                          reduction = "harmony",
                          dims = 1:30)

all_merged_subset<-FindClusters(all_merged_subset, resolution = 0.3,
                         cluster.name = "harmony_clusters",group.singletons = TRUE,
                         method = "igraph")

all_merged_subset<-RunUMAP(all_merged_subset, reduction = "harmony",
                    dims = 1:30,reduction.name = "umap.harmony")

all_merged_subset<-RunTSNE(all_merged_subset, reduction = "harmony",
                    dims = 1:30,reduction.name = "tsne.harmony")


DimPlot(all_merged_subset,
        reduction = "umap.harmony",
        label = T,
        repel = T, 
        #pt.size = 1.2,
        #group.by = c("predicted.celltype"),
        #alpha = 0.9,
        #split.by=c("HIV_Status")
)+
  theme_bw()#+NoLegend()

VlnPlot(all_merged_subset,
        features = c("CD3D","CSF3R","CD19","SCEL","FOXJ1","CD79A","MUC5AC","KRT15","KRT5"))
FeaturePlot(all_merged_subset,
            reduction = "umap.harmony",
            features = c("CD3D","CSF3R","CD19","SCEL","FOXJ1","CD79A","MUC5AC","KRT15","KRT5"))

# Save all_merged_subset -----
save(all_merged_subset,file = "Data/Single_Cell_Data/all_merged_subset.RData")
load("Data/Single_Cell_Data//all_merged_subset.RData")

# Finding all markers -----
all_merged_subset<-JoinLayers(all_merged_subset)
all_merged_subset_markers <- FindAllMarkers(
  all_merged_subset,
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

# Save all_merged subset markers ----
write.csv(all_merged_subset_markers,"scRNAseq_Results/all_merged_subset_markers.csv")

# Rename idents of epithelial cell object clusters ----
all_merged_subset_labelled <- RenameIdents(object = all_merged_subset,
                                          "0" = "Secretory cells",
                                          "1" = "Goblet cells",
                                          "2" = "CD3+ T cells",
                                          "3" = "Neurons",
                                          "4" = "Goblet cells",
                                          "5" = "Phagocytes",
                                          "6" = "Basal cells",
                                          "7" = "Developing ciliated cells",
                                          "8" = "Neutrophils",
                                          "9" = "Goblet cells",
                                          "10" = "Squamous cells",
                                          "11" = "FOXJ1++ Ciliated cells",
                                          "12" = "Stressed cells",
                                          "13" = "B cells",
                                          "14" = "Club cells",
                                          "15" = "Deuterosomal cells",
                                          "16" = "BEST++ Cilia++ Ciliated cells",
                                          "17" = "Ionocytes",
                                          "18" = "Dendritic cells")

# finding markers for all_merged subset labelled object
all_merged_subset_labelled_markers <- FindAllMarkers(
  all_merged_subset_labelled,
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



# Save all_merged_subset_labelled ----
save(all_merged_subset_labelled,file = "data/all_merged_subset_labelled.RData")
load("Data/Single_Cell_Data/all_merged_subset_labelled.RData")

DimPlot(all_merged_subset_labelled,
        reduction = "umap.harmony")
# Save all_merged subset labelled markers ----
write.csv(all_merged_subset_labelled_markers,"results/all_merged_subset_labelled_markers.csv")

Markers <- FindAllMarkers(
  all_merged_subset_labelled,
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

write.csv(Markers,'scRNAseq_Results/Markers.csv')

# To get number of cells in the seurat object -----

ncells <- ncol(all_merged_subset_labelled)

# To get number of genes in the seurat object ----

ngenes <- nrow(all_merged_subset_labelled)


# To get average number of genes in the seurat object -----

average_genes_per_cell <- mean(Matrix::rowSums(all_merged_subset_labelled@assays$RNA@layers$counts))



# Dimplot----
load("data/all_merged_subset_labelled.RData")
saveRDS(all_merged_subset_labelled, file = 'seurat.rds')
main_umap <- DimPlot(all_merged_subset_labelled,
        reduction = "umap.harmony",
        label = T,
        repel = T, 
        label.size = 5,
        #pt.size = 1.2,
        #group.by = c("predicted.celltype"),
        #alpha = 0.9,
        #split.by=c("sample")
)+
  theme_bw()+NoLegend()+
  labs(x="UMAP-1",y="UMAP-2",title = "")+
  theme(#axis.title = element_text(size = 10,face = "bold"),
        axis.title = element_blank(),
        #panel.border = element_rect(linewidth = 2),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())#+
  #geom_hline(yintercept = 0,color="lightblue",linetype='dashed')+
  #geom_vline(xintercept = 0,color="lightblue",linetype='dashed')
main_umap

main_umap <- as.ggplot(main_umap)
# Save the ggplot object
ggsave(here("outputs", "main_umap.pdf"),
       plot = ((main_umap|plot_layout(ncol = 3, nrow = 2, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 7, unit="in", dpi = 700,limitsize = T)


ggsave(here("outputs", "Fig1a.pdf"),
       plot = ((main_umap|plot_layout(ncol = 3, nrow = 2, width = c(1,1,1,1)))),
       width = 20, height =7, unit="in", dpi = 700)


# Dimplot of grouped by HIVstatus
main_umap_HIV <- DimPlot(all_merged_subset_labelled,
                              reduction = "umap.harmony",
                              label = F,
                              repel = T, 
                              #pt.size = 1.2,
                              group.by = c("HIV_Status"),
                              #alpha = 0.9,
                              #split.by=c("HIV_Status")
)+
  theme_bw()+NoLegend()+
  labs(x="UMAP-1",y="UMAP-2",title = "")+
  theme(#axis.title = element_text(size = 10,face = "bold"),
    axis.title = element_blank(),
    #panel.border = element_rect(linewidth = 2),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank())#+
#geom_hline(yintercept = 0,color="lightblue",linetype='dashed')+
#geom_vline(xintercept = 0,color="lightblue",linetype='dashed')
main_umap_HIV

main_umap_HIV <- as.ggplot(main_umap_HIV)
# Save the ggplot object
ggsave(here("Thesis Figures", "main_umap_HIV.pdf"),
       plot = ((main_umap_HIV|plot_layout(ncol = 3, nrow = 2, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 7, unit="in", dpi = 700,limitsize = F)

ggsave(here("Outputs", "main_umap_HIV.pdf"),
       plot = ((main_umap_HIV|plot_layout(ncol = 3, nrow = 2, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 7, unit="in", dpi = 700,limitsize = F)

# Dimplot of grouped by sex
main_umap_sex <- DimPlot(all_merged_subset_labelled,
                     reduction = "umap.harmony",
                     label = F,
                     repel = T, 
                     #pt.size = 1.2,
                     group.by = c("sex"),
                     #alpha = 0.9,
                     #split.by=c("HIV_Status")
)+
  theme_bw()+NoLegend()+
  labs(x="UMAP-1",y="UMAP-2",title = "")+
  theme(#axis.title = element_text(size = 10,face = "bold"),
    axis.title = element_blank(),
    #panel.border = element_rect(linewidth = 2),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank())#+
#geom_hline(yintercept = 0,color="lightblue",linetype='dashed')+
#geom_vline(xintercept = 0,color="lightblue",linetype='dashed')
main_umap_sex

main_umap_sex <- as.ggplot(main_umap_sex)
# Save the ggplot object
ggsave(here("Outputs", "main_umap_sex.pdf"),
       plot = ((main_umap_sex|plot_layout(ncol = 3, nrow = 2, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 7, unit="in", dpi = 700,limitsize = F)

# Dimplot of grouped by Carriage status
main_umap_Carriage <- DimPlot(all_merged_subset_labelled,
                         reduction = "umap.harmony",
                         label = F,
                         repel = T, 
                         #pt.size = 1.2,
                         group.by = c("Carriage_Status"),
                         #alpha = 0.9,
                         #split.by=c("HIV_Status")
)+
  theme_bw()+NoLegend()+
  labs(x="UMAP-1",y="UMAP-2",title = "")+
  theme(#axis.title = element_text(size = 10,face = "bold"),
    axis.title = element_blank(),
    #panel.border = element_rect(linewidth = 2),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank())#+
#geom_hline(yintercept = 0,color="lightblue",linetype='dashed')+
#geom_vline(xintercept = 0,color="lightblue",linetype='dashed')
main_umap_Carriage

main_umap_Carriage <- as.ggplot(main_umap_Carriage)
# Save the ggplot object
ggsave(here("Outputs", "main_umap_Carriage.png"),
       plot = ((main_umap_Carriage|plot_layout(ncol = 3, nrow = 2, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 7, unit="in", dpi = 700,limitsize = F)

# All umaps in one
ggsave(here("Outputs", "umap.pdf"),
       plot = ((main_umap|main_umap_HIV|main_umap_sex|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 7, unit="in", dpi = 700,limitsize = F)


FeaturePlot(all_merged_subset_labelled,
            reduction = "umap.harmony",
            features = c("KRT7","CXCL17","F3","AQP5","CP"))

VlnPlot(all_merged_subset_labelled,
        features = c("CAPS","CETN2","MORN2",
                     "C9orf24","TPPP3","C20orf85",
                     "RSPH1","C1orf194","PTAFR","PIGR"),
        ncol = 5)


VlnPlot(all_merged_subset_labelled,
        features = c("SCEL",""),
        #idents = c("CD3+ T cells"),
        #pt.size = 0.0,
        #split.plot = F,
        log = F,
        #y.max = 2.05,
        #cols = c("green","red","blue"),
        #split.by = c("HIV_Status")
        )+
  NoLegend()#+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
  
all_merged_subset_labelled$cell_clusters <- paste0(all_merged_subset_labelled@active.ident)  

all_merged_subset_labelled_sce <- as.SingleCellExperiment(all_merged_subset_labelled)

cell_counts <- as.data.frame(table(all_merged_subset_labelled_sce$cell_clusters))
  
cluster_heatmap <- DoHeatmap(subset(all_merged_subset_labelled,
                 downsample=120),
          features = c("CST1","ID1","STATH","CYP2F1","HSD11B2","MUC5AC",
                       "SERPINB3","ALOX15","SERPINB4","MT1X","SLC26A2",
                       "CD3G","CCL5","GZMA","CD96","CD7",
                       "AC005062.1","CATSPERB","AC083837.1","LINC00342","GRAMD1B",
                       "IL32","AQP4-AS1","BUB3","COA1","FAM13A",
                       "C1QC","HLA-DPB1","LGALS1","HMOX2","LYZ",
                       "KRT5","CSRP2","TP63","FGFR3","PCP4L1",
                       "OMG","TMEM190","IGFBP7","C20orf85","TPPP3",
                       "CSF3R","G0S2","CXCL8","IL1R2","PROK2",
                       "AC092691.1","PTPRG","LSAMP","CSMD1","SHISA9",
                       "SPRR3","CNFN","MAL","SPRR2D","SPRR2A",
                       "CFAP157","RP1","DNAH3","CDHR3","VWA3A",
                       "MS4A1","CD79A","IGKC","CD79B","BANK1",
                       "BPIFA1","DUOXA2","C3","AZGP1","BPIFB1",
                       "CDC20B","CCNO","DEUP1","E2F7","CCDC74A",
                       "AC008415.1","AC068587.4","CFAP65","SHANK2","AQP4-AS1",
                       "DGKI","PDE1C","ASCL3","TMEM61","CFTR",
                       "SLC8A1","THBS1","CCL2","PLXDC2","CXCL10"),
          slot = "scale.data",
          group.by = "ident",
          group.bar = T,
          angle = 90,
          size = 4)+
  scale_fill_viridis()+
  NoLegend()
cluster_heatmap

cluster_heatmap <- as.ggplot(cluster_heatmap)
# Save the ggplot object
ggsave(here("Outputs", "cluster_heatmap.png"),
       plot = ((cluster_heatmap|plot_layout(ncol = 1, nrow = 2, width = c(1,1,1,1,1,1,1,1)))),
       width = 15, height = 14, unit="in", dpi = 300,limitsize = T)

# Save all_merged_subset_labelled ----
save(all_merged_subset_labelled,file = "data/all_merged_subset_labelled.RData")
load("data/all_merged_subset_labelled.RData")
# Mapping seurat object against a reference object (SHALEK LAB DATA)----
# Load reference seurat object
nasal_reference <- readRDS("Olympia/nasal_reference_REAL.rds")

view(nasal_reference@meta.data) 
# Taking a look
DimPlot(nasal_reference, group.by = "Detailed_Cell_Annotations",
        reduction = "umap",label = T)+theme_bw()+
  NoLegend()


# Mapping
anchors <- FindTransferAnchors(reference = nasal_reference,
                               query = all_merged_subset_labelled,
                               normalization.method = "SCT",
                               reference.reduction = "pca",
                               dims = 1:30)


all_merged_subset_labelled <- MapQuery(anchorset = anchors,
                              reference = nasal_reference,
                              query = all_merged_subset_labelled,
                              refdata = list(celltype="Detailed_Cell_Annotations"),
                              reference.reduction = "pca",
                              reduction.model = "umap")

view(all_merged_subset@meta.data)

# ----Dimplot
DimPlot(all_merged_subset_labelled,reduction = "umap.harmony",
        group.by = "predicted.celltype",label = F, repel = T,
        #split.by = ("HIV_Status")
)+theme_bw()+scale_fill_viridis()#+
NoLegend()

# Remove stressed cells ----


all_merged_subset_labelled_final <- subset(all_merged_subset_labelled,
                                           idents = c("Stressed cells"), 
                                           invert = T)
DimPlot(all_merged_subset_labelled_final,
        reduction = "umap.harmony",
        label = T)+
  theme_base()+
  theme(legend.position = "none")

# Save all merged subset labelled final object----

save(all_merged_subset_labelled_final, file = "data/all_merged_subset_labelled_final.RData")
load("data/all_merged_subset_labelled_final.RData")
# Dimplot of all_merged_subset_labelled_final

all_merged_umap <- DimPlot(all_merged_subset_labelled_final,
                           reduction = "umap.harmony",
                           #group.by = "predicted.celltype",
                           label = T, 
                           repel = T,
                           #pt.size = 2,
                           label.size = 6,
                           #split.by = ("HIV_Status")
                           )+
  labs(title = "",
       x="UMAP 1",
       y="UMAP 2")+
  theme(#axis.title = element_text(size = 20),
    #panel.border = element_rect(linewidth = 2),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank())+
  geom_hline(yintercept = 0,color="grey")+
  geom_vline(xintercept = 0,color="grey")+
  NoLegend()

all_merged_umap

# comvert dimplot to ggplot
all_merged_umap <- as.ggplot(all_merged_umap)

# save the saved converted dimplot 
ggsave(here("results","all_merged_umap.png"),
       plot = ((all_merged_umap|plot_layout(ncol = 3, nrow = 2, width = c(1,1,1,1,1,1,1,1)))),
       width = 33, height = 12, unit="in", dpi = 700,limitsize = F)



# Vlnplot
expression_allcellsln <- VlnPlot(all_merged_subset_labelled_final,
                                #idents = c("Monocytes/Macrophages","Neutrophils"),
                                features = c("nGene"),
                                #split.by = "HIV_Status",
                                split.plot = F,
                                #group.by = "Carriage_Status",
                                slot = "data",
                                log = F,
                                pt.size = 0.001)+
  theme_classic()+
  labs(y="Expression")+
  #scale_y_continuous(limits = c(0,8000))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold",
                                    size=25),
        plot.title = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   face = "bold",size=25))
expression_allcellsln




# comvert dimplot to ggplot
expression_allcellsln <- as.ggplot(expression_allcellsln)

# save the saved converted dimplot 
ggsave(here("results", "expression_allcellsln.png"),
       plot = ((expression_allcellsln|plot_layout(ncol = 2, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 25, height = 10, unit="in", dpi = 700,limitsize = F)



# Differential abundance ----
df_comp <- as.data.frame.matrix(table(all_merged_subset_labelled_final$sample, 
                                      all_merged_subset_labelled_final@active.ident))
df_comp$sample <- rownames(df_comp)

# Save

write.csv(df_comp,"results/df_comp_HIV.csv")


df_comp %>%
  pivot_longer(cols = c("Secretory cells":"Dendritic cells"),
               names_to = "Cell Type",
               values_to = "Frequency") %>%
  ggplot(aes(x=`sample`,y=Frequency,fill=`Cell Type`))+
  geom_bar(stat = "identity",position = "fill")+
  scale_fill_manual(values=distinct_palette(n=16,pal = "brewerPlus",add = "lightgrey"))+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = "right",
        axis.text = element_text(size = 20,face = "bold"),
        axis.title = element_text(size = 20,face = "bold"),
        legend.text = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(angle = 45,hjust = 1))


# Create HIV status column
df_comp$HIV_Status <- NA
df_comp$HIV_Status[which(str_detect(df_comp$sample, "^CUH124"))] <- "HIV+ ART>1 Year"
df_comp$HIV_Status[which(str_detect(df_comp$sample, "^CUG11X"))] <- "HIV+ ART>1 Year"
df_comp$HIV_Status[which(str_detect(df_comp$sample, "^CUF131"))] <- "HIV-"
df_comp$HIV_Status[which(str_detect(df_comp$sample, "^CUF130"))] <- "HIV+ ART>1 Year"
df_comp$HIV_Status[which(str_detect(df_comp$sample, "^CUF134"))] <- "HIV+ ART<3 Months"
df_comp$HIV_Status[which(str_detect(df_comp$sample, "^CUF13I"))] <- "HIV+ ART>1 Year"
df_comp$HIV_Status[which(str_detect(df_comp$sample, "^CUF135"))] <- "HIV+ ART>1 Year"
df_comp$HIV_Status[which(str_detect(df_comp$sample, "^CUF136"))] <- "HIV-"
df_comp$HIV_Status[which(str_detect(df_comp$sample, "^CUF13J"))] <- "HIV+ ART<3 Months"
df_comp$HIV_Status[which(str_detect(df_comp$sample, "^CUF13K"))] <- "HIV+ ART<3 Months"
df_comp$HIV_Status[which(str_detect(df_comp$sample, "^CUF137"))] <- "HIV-"
#df_comp$HIV_Status[which(str_detect(df_comp$sample, "^CUF12Y"))] <- "HIV+ ART>1 Year"


# Create Carriage status column
df_comp$Carriage_Status <- NA
df_comp$Carriage_Status[which(str_detect(df_comp$sample, "^CUH124"))] <- "SPN+"
df_comp$Carriage_Status[which(str_detect(df_comp$sample, "^CUG11X"))] <- "SPN-"
df_comp$Carriage_Status[which(str_detect(df_comp$sample, "^CUF131"))] <- "SPN+"
df_comp$Carriage_Status[which(str_detect(df_comp$sample, "^CUF130"))] <- "SPN-"
df_comp$Carriage_Status[which(str_detect(df_comp$sample, "^CUF134"))] <- "SPN+"
df_comp$Carriage_Status[which(str_detect(df_comp$sample, "^CUF13I"))] <- "SPN+"
df_comp$Carriage_Status[which(str_detect(df_comp$sample, "^CUF135"))] <- "SPN+"
df_comp$Carriage_Status[which(str_detect(df_comp$sample, "^CUF136"))] <- "SPN+"
df_comp$Carriage_Status[which(str_detect(df_comp$sample, "^CUF13J"))] <- "SPN+"
df_comp$Carriage_Status[which(str_detect(df_comp$sample, "^CUF13K"))] <- "SPN-"
df_comp$Carriage_Status[which(str_detect(df_comp$sample, "^CUF137"))] <- "SPN+"
#df_comp$Carriage_Status[which(str_detect(df_comp$sample, "^CUF12Y"))] <- "SPN+"

# Create Carriage density column
df_comp$Carriage_density <- NA
df_comp$Carriage_density[which(str_detect(df_comp$sample, "^CUH124"))] <- 5862500
df_comp$Carriage_density[which(str_detect(df_comp$sample, "^CUG11X"))] <- NA
df_comp$Carriage_density[which(str_detect(df_comp$sample, "^CUF131"))] <- 26800000
df_comp$Carriage_density[which(str_detect(df_comp$sample, "^CUF130"))] <- NA
df_comp$Carriage_density[which(str_detect(df_comp$sample, "^CUF134"))] <- 686750
df_comp$Carriage_density[which(str_detect(df_comp$sample, "^CUF13I"))] <- 670
df_comp$Carriage_density[which(str_detect(df_comp$sample, "^CUF135"))] <- 33500
df_comp$Carriage_density[which(str_detect(df_comp$sample, "^CUF136"))] <- 67
df_comp$Carriage_density[which(str_detect(df_comp$sample, "^CUF13J"))] <- 26800
df_comp$Carriage_density[which(str_detect(df_comp$sample, "^CUF13K"))] <- NA
df_comp$Carriage_density[which(str_detect(df_comp$sample, "^CUF137"))] <- 10720
#df_comp$Carriage_density[which(str_detect(df_comp$sample, "^CUF12Y"))] <- 2345

# Save. the df_comp
write.csv(df_comp,"results/df_comp.csv")

df_comp %>%
  filter(Carriage_Status=="SPN+") %>%
  ggplot(aes(Carriage_density,`Secretory cells`))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_x_log10()+
  stat_cor(method = "spearman")

df_comp %>%
  ggplot(aes(HIV_Status,`Neutrophils`))+
  geom_boxplot()+
  geom_jitter()+
  geom_pwc()+
  theme_base(base_size = 25)

df_comp %>%
  pivot_longer(cols = c("Secretory cells":"Dendritic cells"),
               names_to = "Sample Type",
               values_to = "Percentage of sample") %>%
  ggplot(aes(x=`Carriage_Status`,y=`Percentage of sample`, fill=`Sample Type`))+
  geom_bar(stat = "identity",position = "fill")+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1))

# load cell cluster relative abundance  csv file

cell_cluster_relative_abundance %>%
  filter(Carriage_Status=="SPN+") %>%
  ggplot(aes(Carriage_density,`Neutrophils`))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_x_log10()+
  stat_cor(method = "spearman")



cell_cluster_relative_abundance %>%
  pivot_longer(cols = c("Secretory cells":"Dendritic cells"),
               names_to = "Sample Type",
               values_to = "Percentage of sample") %>%
  ggplot(aes(x=`Sample Type`,y=`Percentage of sample`, fill=`HIV_Status`))+
  geom_bar(stat = "identity",position = "fill")+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1))

# Differential abundance using milo----
# convert seurat to a single celll experiment
all_merged_sce <- as.SingleCellExperiment(all_merged_subset_labelled_final)
all_merged_sce <- logNormCounts(all_merged_sce)

# Create a milo object
all_merged_milo <- Milo (all_merged_sce)

# construct KNN graph
all_merged_milo <- buildGraph(all_merged_milo,k=10,d=30)

all_merged_milo <- makeNhoods(all_merged_milo,prop = 0.1,
                              k=10,d=30)
plotNhoodSizeHist(all_merged_milo)

# Counting cells in neighbourhoods
all_merged_milo <- countCells(all_merged_milo,
                              meta.data = data.frame(colData(all_merged_milo)),
                              sample="sample")
view(all_merged_milo@colData)

head(nhoodCounts(all_merged_milo))
# Differential abundance testing

all_merged_design <- data.frame(colData(all_merged_milo))[,c("sample","HIV_Status","Carriage_Status","Carriage_density")]
all_merged_design <- distinct(all_merged_design)
all_merged_design


# store the distances between nearest neighbors in the Milo object.
all_merged_milo <- calcNhoodDistance(all_merged_milo,d=30)

# Now we can do the test, explicitly defining our experimental design.
da_results <- testNhoods(all_merged_milo,
                         design = ~ rownames(colData(all_merged_milo)),
                         design.df = all_merged_design)

da_results %>%
  arrange(- SpatialFDR) %>%
  head()



# New abundance testing ----
pt_HIV <- table(Idents(all_merged_subset_labelled),all_merged_subset_labelled$sample)
pt_HIV <- as.data.frame(pt_HIV)




# -----save
write.csv(pt_HIV,"results/pt_HIV_sample.csv")


pt_HIV_sample %>%
  filter(Carriage_Status=='SPN+' & `Cell Type`=='Neurons') %>%
  ggplot(aes(Carriage_density,Propotion))+
  geom_point(size=5,aes(color=HIV_Status))+
  geom_smooth(method = "lm")+
  scale_x_log10()+
  stat_cor(method = 'spearman',label.y = 0.4)+
  theme(legend.position = "none")


cell_prop_hiv <- pt_HIV_sample %>%
  filter(`Cell Type`!=146)%>%
  group_by(HIV_Status) %>%
  ggplot(aes(`Cell Type`,Propotion,color=HIV_Status))+
  geom_boxplot(outlier.shape = NA,linewidth=4)+
  geom_jitter(size=10,width = 0.2)+
  theme_bw()+
  labs(x="",y="Proportion of cells")+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 50),
        panel.border = element_rect(linewidth = 8))#+
  #scale_y_log10()
cell_prop_hiv

ggsave(here("results", "cell_prop_hiv.pdf"),
       plot = ((cell_prop_hiv|plot_layout(ncol = 1, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 50, height = 18, unit="in", dpi = 700,limitsize = F)  



# adding metadata 

# Create HIV status column
pt_HIV$HIV_Status <- NA
pt_HIV$HIV_Status[which(str_detect(pt_HIV$Var2, "^CUH124"))] <- "HIV+ ART>1 Year"
pt_HIV$HIV_Status[which(str_detect(pt_HIV$Var2, "^CUG11X"))] <- "HIV+ ART>1 Year"
pt_HIV$HIV_Status[which(str_detect(pt_HIV$Var2, "^CUF131"))] <- "HIV-"
pt_HIV$HIV_Status[which(str_detect(pt_HIV$Var2, "^CUF130"))] <- "HIV+ ART>1 Year"
pt_HIV$HIV_Status[which(str_detect(pt_HIV$Var2, "^CUF134"))] <- "HIV+ ART<3 Months"
pt_HIV$HIV_Status[which(str_detect(pt_HIV$Var2, "^CUF13I"))] <- "HIV+ ART>1 Year"
pt_HIV$HIV_Status[which(str_detect(pt_HIV$Var2, "^CUF135"))] <- "HIV+ ART>1 Year"
pt_HIV$HIV_Status[which(str_detect(pt_HIV$Var2, "^CUF136"))] <- "HIV-"
pt_HIV$HIV_Status[which(str_detect(pt_HIV$Var2, "^CUF13J"))] <- "HIV+ ART<3 Months"
pt_HIV$HIV_Status[which(str_detect(pt_HIV$Var2, "^CUF13K"))] <- "HIV+ ART<3 Months"
pt_HIV$HIV_Status[which(str_detect(pt_HIV$Var2, "^CUF137"))] <- "HIV-"
#pt_HIV$HIV_Status[which(str_detect(pt_HIV$Var2, "^CUF12Y"))] <- "HIV+ ART>1 Year"


# Create Carriage status column
pt_HIV$Carriage_Status <- NA
pt_HIV$Carriage_Status[which(str_detect(pt_HIV$Var2, "^CUH124"))] <- "SPN+"
pt_HIV$Carriage_Status[which(str_detect(pt_HIV$Var2, "^CUG11X"))] <- "SPN-"
pt_HIV$Carriage_Status[which(str_detect(pt_HIV$Var2, "^CUF131"))] <- "SPN+"
pt_HIV$Carriage_Status[which(str_detect(pt_HIV$Var2, "^CUF130"))] <- "SPN-"
pt_HIV$Carriage_Status[which(str_detect(pt_HIV$Var2, "^CUF134"))] <- "SPN+"
pt_HIV$Carriage_Status[which(str_detect(pt_HIV$Var2, "^CUF13I"))] <- "SPN+"
pt_HIV$Carriage_Status[which(str_detect(pt_HIV$Var2, "^CUF135"))] <- "SPN+"
pt_HIV$Carriage_Status[which(str_detect(pt_HIV$Var2, "^CUF136"))] <- "SPN+"
pt_HIV$Carriage_Status[which(str_detect(pt_HIV$Var2, "^CUF13J"))] <- "SPN+"
pt_HIV$Carriage_Status[which(str_detect(pt_HIV$Var2, "^CUF13K"))] <- "SPN-"
pt_HIV$Carriage_Status[which(str_detect(pt_HIV$Var2, "^CUF137"))] <- "SPN+"
#pt_HIV$Carriage_Status[which(str_detect(pt_HIV$Var2, "^CUF12Y"))] <- "SPN+"

# Create Carriage density column
pt_HIV$Carriage_density <- NA
pt_HIV$Carriage_density[which(str_detect(pt_HIV$Var2, "^CUH124"))] <- 5862500
pt_HIV$Carriage_density[which(str_detect(pt_HIV$Var2, "^CUG11X"))] <- NA
pt_HIV$Carriage_density[which(str_detect(pt_HIV$Var2, "^CUF131"))] <- 26800000
pt_HIV$Carriage_density[which(str_detect(pt_HIV$Var2, "^CUF130"))] <- NA
pt_HIV$Carriage_density[which(str_detect(pt_HIV$Var2, "^CUF134"))] <- 686750
pt_HIV$Carriage_density[which(str_detect(pt_HIV$Var2, "^CUF13I"))] <- 670
pt_HIV$Carriage_density[which(str_detect(pt_HIV$Var2, "^CUF135"))] <- 33500
pt_HIV$Carriage_density[which(str_detect(pt_HIV$Var2, "^CUF136"))] <- 67
pt_HIV$Carriage_density[which(str_detect(pt_HIV$Var2, "^CUF13J"))] <- 26800
pt_HIV$Carriage_density[which(str_detect(pt_HIV$Var2, "^CUF13K"))] <- NA
pt_HIV$Carriage_density[which(str_detect(pt_HIV$Var2, "^CUF137"))] <- 10720
#pt_HIV$Carriage_density[which(str_detect(pt_HIV$Var2, "^CUF12Y"))] <- 2345


# All cells in the clusters (frequency by sample)

distinct_palette(n=20,pal="brewerPlus")

CellType.sample <- ggplot(pt_HIV, aes(x = Var2, y = Freq, fill = Var1)) +
  #theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("Sample") +
  ylab("% of sample") +
  theme_bw(base_size = 30)+
  #facet_wrap(~ HIV_Status, scales = "free_x")+
  #scale_fill_manual(values = brewer.pal(16,"kelly")) +
  scale_fill_manual(values=distinct_palette(n=20,pal = "brewerPlus",add = "lightgrey"))+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle=45,hjust = 1,
                                   size = 20,face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 2),
        axis.title = element_text(size = 30))
CellType.sample
# Save celltype frequency in each sample plot ----

ggsave(here("results", "CellType.sample.tiff"),
       plot = ((CellType.sample|plot_layout(ncol = 2, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 30, height = 10, unit="in", dpi = 700,limitsize = F)

# All cells in the clusters (frequency by sample)

distinct_palette(n=20,pal="brewerPlus")

CellType.HIV <- ggplot(pt_HIV, aes(x = HIV_Status, y = Freq, fill = Var1)) +
  #theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("Sample") +
  ylab("% of sample") +
  theme_bw(base_size = 30)+
  #facet_wrap(~ HIV_Status, scales = "free_x")+
  #scale_fill_manual(values = brewer.pal(16,"kelly")) +
  scale_fill_manual(values=distinct_palette(n=20,pal = "brewerPlus",add = "lightgrey"))+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle=45,hjust = 1,
                                   size = 20,face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 2),
        axis.title = element_text(size = 30))
CellType.HIV
# Save celltype frequency in each sample plot ----

ggsave(here("results", "CellType.HIV.tiff"),
              plot = ((CellType.HIV|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
              width = 15, height = 10, unit="in", dpi = 700,limitsize = F)

immune_cells <- pt_HIV %>%
  filter(Var1=="CD3+ T cells" | Var1=="Phagocytes" |
           Var1=="Neutrophils" |
           Var1=="B cells" |
           Var1=="Dendritic cells") %>%
  ggplot(aes(x=Var2,y=Freq,fill=Var1))+
  #theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Proportion") +
  facet_wrap(~HIV_Status,scales = "free_x")+
  scale_fill_manual(values = brewer.pal(12,"Paired")) +
  theme_bw(base_size = 25)+
  theme(legend.title = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 45,hjust = 1))
immune_cells



immune_cells <- pt_HIV %>%
  filter(Var1=="CD3+ T cells" | 
           Var1=="Phagocytes" |
           Var1=="Neutrophils" |
           Var1=="B cells" |
         Var1=="Dendritic cells") %>%
  ggplot(aes(x=HIV_Status,y=Freq,fill=Var1))+
  #theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12,"Paired")) +
  theme_bw(base_size = 25)+
  theme(legend.title = element_blank(),
        legend.position = "right",
        axis.text = element_text(size = 50,face = "bold"),
        axis.title = element_text(size = 50,face = "bold"),
        legend.text = element_text(size = 50,face = "bold"),
        axis.text.x = element_text(angle = 45,hjust = 1))
immune_cells


epithelial_cells <- pt_HIV %>%
  filter(Var1!="CD3+ T cells" & 
           Var1!="Phagocytes" &
           Var1!="Neutrophils" &
           Var1!="B cells" &
           Var1!="Dendritic cells") %>%
  ggplot(aes(x=HIV_Status,y=Freq,fill=Var1))+
  theme_bw(base_size = 25) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(11,"Paired")) +
  theme(legend.title = element_blank(),
        legend.position = "right",
        axis.text = element_text(size = 50,face = "bold"),
        axis.title = element_text(size = 50,face = "bold"),
        legend.text = element_text(size = 50,face = "bold"),
        axis.text.x = element_text(angle = 45,hjust = 1))
  
epithelial_cells

ggsave(here("results", "cluster_abundance.pdf"),
       plot = ((immune_cells|epithelial_cells|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 50, height = 18, unit="in", dpi = 700,limitsize = F)

phagocytes_density <- pt_HIV %>%
  filter(Carriage_Status=="SPN+" & Var1=="Phagocytes")%>%
  ggplot(aes(Carriage_density,Freq)) +
  geom_point(size=5, aes(color=HIV_Status))+
  geom_smooth(method = "lm")+
  scale_x_log10()+
  stat_cor(method = "spearman",label.y = 350,size=15)+
  theme_bw(base_size = 25)+
  theme(legend.position = "none")
phagocytes_density  


Bcells_density <- pt_HIV %>%
  filter(Carriage_Status=="SPN+" & Var1=="Dendritic cells")%>%
  ggplot(aes(Carriage_density,Freq)) +
  geom_point(size=5, aes(color=HIV_Status))+
  geom_smooth(method = "lm")+
  scale_x_log10()+
  stat_cor(method = "spearman",label.y = 40,size=15)+
  theme_bw(base_size = 25)+
  theme(legend.position = "none")
Bcells_density

ggsave(here("results", "phagocytes_density.pdf"),
       plot = ((phagocytes_density|epithelial_cells|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 50, height = 18, unit="in", dpi = 700,limitsize = F)  

# save the abundance ggplots



# ----- Subsetting immune cells-----
all_merged_immune <- subset(all_merged_subset_labelled_final,
                            idents = c("B cells","CD3+ T cells","Phagocytes",
                                       "Neutrophils","Dendritic cells"), invert = F)

all_merged_immune$Cell_CLusters <- paste0(all_merged_immune@active.ident)

view(all_merged_immune@meta.data)
# Dimplot all merged immune

DimPlot(all_merged_immune,
        reduction = "umap.harmony")


all_merged_immune<-IntegrateLayers(object=all_merged_immune,
                                   method=HarmonyIntegration,
                                   group.by = "sample",
                                   dims = 1:20,
                                   orig.reduction ="pca",
                                   new.reduction="harmony",
                                   verbose = FALSE)


all_merged_immune<-FindNeighbors(all_merged_immune, 
                                 reduction = "harmony",
                                 dims = 1:20)

all_merged_immune<-FindClusters(all_merged_immune, resolution = 0.3,
                                cluster.name = "harmony_clusters",group.singletons = TRUE,
                                method = "igraph")

all_merged_immune<-RunUMAP(all_merged_immune, reduction = "harmony",
                           dims = 1:30,reduction.name = "umap.harmony")

all_merged_immune<-RunTSNE(all_merged_immune, reduction = "harmony",
                           dims = 1:30,reduction.name = "tsne.harmony")


# Find all markers for immune cells ----
all_merged_immune <- JoinLayers(all_merged_immune)
all_merged_immune_markers <- FindAllMarkers(
  all_merged_immune,
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


# save all_merged_immune_markers----
write.csv(all_merged_immune_markers,
          "results/all_merged_immune_markers.csv")




# Dimplot all all merged immune

all_merged_immune_umap <- DimPlot(all_merged_immune,
                           reduction = "umap.harmony",
                           label = T, 
                           repel = T,
                           #pt.size = 2,
                           #label.size = 20,
                           #group.by = ("Cell_CLusters"),
                           #split.by = ('HIV_Status')
)+
  labs(title = "",
       x="UMAP 1",y="UMAP 2")+
  theme_bw(base_size = 15)
all_merged_immune_umap

# comvert dimplot to ggplot
all_merged_immune_umap <- as.ggplot(all_merged_immune_umap)

# save the saved converted dimplot 
ggsave(here("results", "all_merged_immune_umap.pdf"),
       plot = ((all_merged_immune_umap|plot_layout(ncol = 1, nrow = 2, width = c(1,1,1,1,1,1,1,1)))),
       width = 45, height = 17, unit="in", dpi = 700,limitsize = F)

VlnPlot(all_merged_immune,
        features = c("CD19","CSF3R","CD2", 
                     "CD3D","CD8A","CD4",
                     "VCAN","FCGR3A","CD14",
                     "KLRB1","PLXDC2",
                     "CD163","C1QC"),
        log = T)
# ---- Save ----
save(all_merged_immune, file = "data/all_merged_immune.RData")
load("data/all_merged_immune.RData")


all_merged_immune <- JoinLayers(all_merged_immune)
all_merged_immune_final <- subset(all_merged_immune,
                            idents = c("0","8"), 
                            invert = T)


all_merged_immune_final <-IntegrateLayers(object=all_merged_immune_final,
                                   method=HarmonyIntegration,
                                   group.by = "sample",
                                   dims = 1:20,
                                   orig.reduction ="pca",
                                   new.reduction="harmony",
                                   verbose = FALSE)


all_merged_immune_final<-FindNeighbors(all_merged_immune_final, 
                                 reduction = "harmony",
                                 dims = 1:20)

all_merged_immune_final<-FindClusters(all_merged_immune_final, 
                                      resolution = 0.2,
                                      cluster.name = "harmony_clusters",
                                      group.singletons = TRUE,
                                      method = "igraph")

all_merged_immune_final<-RunUMAP(all_merged_immune_final, reduction = "harmony",
                           dims = 1:30,reduction.name = "umap.harmony")

all_merged_immune_final<-RunTSNE(all_merged_immune_final, reduction = "harmony",
                           dims = 1:30,reduction.name = "tsne")

DimPlot(all_merged_immune_final,
        reduction = "umap.harmony",
        label = T)

# finding all merged immune final markers ----
all_merged_immune_final_markers <- FindAllMarkers(
  all_merged_immune_final,
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

# save all_merged_immune_markers----
write.csv(all_merged_immune_final_markers,
          "results/all_merged_immune_final_markers.csv")
# Rename idents for all_merged_immune_final----
all_merged_immune_final_labelled <- RenameIdents(object = all_merged_immune_final,
                                           "0" = "T cells",
                                           "1" = "Neutrophils",
                                           "2" = "Monocytes/Macrophages",
                                           "3" = "NKT cells",
                                           "4" = "B cells",
                                           "5" = "Dendritic cells",
                                           "6" = "T cells")


# Save all_merged_immune_final_labelled object-----

save(all_merged_immune_final_labelled, file = "data/all_merged_immune_final_labelled.RData")
load("data/all_merged_immune_final_labelled.RData")

immune_cell_labelled_final_umap <- DimPlot(all_merged_immune_final_labelled,
                                           reduction = "umap.harmony",
                                           label = T,
                                           repel = T, 
                                           label.size = 5,
                                           #pt.size = 1.2,
                                           #group.by = c("predicted.celltype"),
                                           #alpha = 0.9,
                                           #split.by=c("sample")
)+
  theme_bw()+NoLegend()+
  labs(x="UMAP-1",y="UMAP-2",title = "")+
  theme(axis.title = element_text(size = 20,
                                  face = "bold"),
        axis.text = element_text(size = 10,
                                 face = "bold"),
        #panel.border = element_rect(linewidth = 2),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())+
  geom_hline(yintercept = 0,color="grey")+
  geom_vline(xintercept = 0,color="grey")

immune_cell_labelled_final_umap


immune_cell_labelled_final_umap <- as.ggplot(immune_cell_labelled_final_umap)
# Save the ggplot object
ggsave(here("results", "immune_cell_labelled_final_umap.pdf"),
       plot = ((immune_cell_labelled_final_umap|plot_layout(ncol = 3, nrow = 2, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 7, unit="in", dpi = 700,limitsize = F)

immune_cluster_final_labelled_heatmap <- DoHeatmap(subset(all_merged_immune_final_labelled,
                                    downsample=3000),
                             features = c("GNLY","GZMA","CCL5","CD3D","CD8A","KLRB1","CD2","NKG7","IL32","CD7", 
                                          "LUCAT1","CSF3R","AQP9","PROK2","FFAR2","IL1R2","CXCL8","FCGR3B","NAMPT","AC099489.1",
                                          "C1QC","C1QA","C1QB","APOC1","APOE","FCER1A","LYZ","TMEM176A","LGALS2","TMEM176B",
                                          "SCGB3A1","NEK10","MTRNR2L12","MTRNR2L8","DUOX2","MT-ND6","SNTN","MUC20","MT-ND5","CD24",
                                          "MS4A1","CD79A","IGKV3-20","BANK1","CD79B","CD19","PAX5","IGHM","SPIB","CD22",
                                          "CPA3","TPSAB1","KIT","GATA2","SLC24A3","SLC18A2","RHEX","THBS1","CCDC88A","CTTNBP2",
                                          "AL133268.4","ITSPAN8","CCN2","LRRIQ1","SELENBP1","ARHGEF28","DNAH11","USP53","ADH7","RBPMS"),
                             slot = "scale.data",
                             group.by = "ident",
                             group.bar = T,
                             angle = 90,
                             size = 4)+
  scale_fill_viridis()+
  #scale_fill_gradientn(colors = c("green", "white", "red"))+
  NoLegend()
immune_cluster_final_labelled_heatmap

immune_cluster_final_labelled_heatmap <- as.ggplot(immune_cluster_final_labelled_heatmap)
# Save the ggplot object
ggsave(here("results", "immune_cluster_final_labelled_heatmap.png"),
       plot = ((immune_cluster_final_labelled_heatmap|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 15, height = 10, unit="in", dpi = 700,limitsize = F)

expression_immunevln <- VlnPlot(all_merged_immune_final_labelled,
        #idents = c("Monocytes/Macrophages","Neutrophils"),
        features = c("nGene"),
        #split.by = "HIV_Status",
        split.plot = F,
        #group.by = "Carriage_Status",
        slot = "data",
        log = T,
        pt.size = 0.001)+
  theme_classic()+
  labs(y="Expression")+
  #scale_y_continuous(limits = c(0,8000))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold",
                                    size=25),
        plot.title = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   face = "bold",size=25))
expression_immunevln


expression_immunevln_myeloid <- VlnPlot(all_merged_immune_final_labelled,
                                idents = c("Monocytes/Macrophages","Neutrophils","Dendritic cells"),
                                features = c("ITGAM","CD63","SELL"),
                                #split.by = "HIV_Status",
                                split.plot = F,
                                #group.by = "HIV_Status",
                                slot = "data",
                                log = T,
                                pt.size = 0.001)+
NoLegend()#+
  #labs(y="Expression")#+
  #scale_y_continuous(limits = c(0,8000))+
  #theme(legend.position = "none",
        #axis.title.x = element_blank(),
        #axis.title.y = element_text(face = "bold",
                                    #size=25),
        #plot.title = element_blank(),
        #axis.text.x = element_text(angle = 45,
                                   #hjust = 1,
                                   #face = "bold",size=25))
expression_immunevln_myeloid


expression_immunevln_myeloid <- as.ggplot(expression_immunevln_myeloid)
# Save the ggplot object
ggsave(here("results", "expression_immunevln_myeloid.pdf"),
       plot = ((expression_immunevln_myeloid|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 7, unit="in", dpi = 700,limitsize = F)




FeaturePlot(all_merged_immune_final_labelled,
            reduction = "umap.harmony",
            features = c("PRF1","GZMB","MRC1"),
            split.by = "HIV_Status")+
  theme(legend.position = "none")




df_comp_immune <- as.data.frame.matrix(table(all_merged_immune_final_labelled$sample, 
                                      all_merged_immune_final_labelled@active.ident))
df_comp_immune$sample <- rownames(df_comp)


# Create HIV status column
df_comp_immune$HIV_Status <- NA
df_comp_immune$HIV_Status[which(str_detect(df_comp_immune$sample, "^CUH124"))] <- "HIV+ ART>1 Year"
df_comp_immune$HIV_Status[which(str_detect(df_comp_immune$sample, "^CUG11X"))] <- "HIV+ ART>1 Year"
df_comp_immune$HIV_Status[which(str_detect(df_comp_immune$sample, "^CUF131"))] <- "HIV-"
df_comp_immune$HIV_Status[which(str_detect(df_comp_immune$sample, "^CUF130"))] <- "HIV+ ART>1 Year"
df_comp_immune$HIV_Status[which(str_detect(df_comp_immune$sample, "^CUF134"))] <- "HIV+ ART<3 Months"
df_comp_immune$HIV_Status[which(str_detect(df_comp_immune$sample, "^CUF13I"))] <- "HIV+ ART>1 Year"
df_comp_immune$HIV_Status[which(str_detect(df_comp_immune$sample, "^CUF135"))] <- "HIV+ ART>1 Year"
df_comp_immune$HIV_Status[which(str_detect(df_comp_immune$sample, "^CUF136"))] <- "HIV-"
df_comp_immune$HIV_Status[which(str_detect(df_comp_immune$sample, "^CUF13J"))] <- "HIV+ ART<3 Months"
df_comp_immune$HIV_Status[which(str_detect(df_comp_immune$sample, "^CUF13K"))] <- "HIV+ ART<3 Months"
df_comp_immune$HIV_Status[which(str_detect(df_comp_immune$sample, "^CUF137"))] <- "HIV-"
#df_comp_immune$HIV_Status[which(str_detect(df_comp_immune$sample, "^CUF12Y"))] <- "HIV+ ART>1 Year"


# Create Carriage status column
df_comp_immune$Carriage_Status <- NA
df_comp_immune$Carriage_Status[which(str_detect(df_comp_immune$sample, "^CUH124"))] <- "SPN+"
df_comp_immune$Carriage_Status[which(str_detect(df_comp_immune$sample, "^CUG11X"))] <- "SPN-"
df_comp_immune$Carriage_Status[which(str_detect(df_comp_immune$sample, "^CUF131"))] <- "SPN+"
df_comp_immune$Carriage_Status[which(str_detect(df_comp_immune$sample, "^CUF130"))] <- "SPN-"
df_comp_immune$Carriage_Status[which(str_detect(df_comp_immune$sample, "^CUF134"))] <- "SPN+"
df_comp_immune$Carriage_Status[which(str_detect(df_comp_immune$sample, "^CUF13I"))] <- "SPN+"
df_comp_immune$Carriage_Status[which(str_detect(df_comp_immune$sample, "^CUF135"))] <- "SPN+"
df_comp_immune$Carriage_Status[which(str_detect(df_comp_immune$sample, "^CUF136"))] <- "SPN+"
df_comp_immune$Carriage_Status[which(str_detect(df_comp_immune$sample, "^CUF13J"))] <- "SPN+"
df_comp_immune$Carriage_Status[which(str_detect(df_comp_immune$sample, "^CUF13K"))] <- "SPN-"
df_comp_immune$Carriage_Status[which(str_detect(df_comp_immune$sample, "^CUF137"))] <- "SPN+"
#df_comp$Carriage_Status[which(str_detect(df_comp$sample, "^CUF12Y"))] <- "SPN+"

# Create Carriage density column
df_comp_immune$Carriage_density <- NA
df_comp_immune$Carriage_density[which(str_detect(df_comp_immune$sample, "^CUH124"))] <- 5862500
df_comp_immune$Carriage_density[which(str_detect(df_comp_immune$sample, "^CUG11X"))] <- NA
df_comp_immune$Carriage_density[which(str_detect(df_comp_immune$sample, "^CUF131"))] <- 26800000
df_comp_immune$Carriage_density[which(str_detect(df_comp_immune$sample, "^CUF130"))] <- NA
df_comp_immune$Carriage_density[which(str_detect(df_comp_immune$sample, "^CUF134"))] <- 686750
df_comp_immune$Carriage_density[which(str_detect(df_comp_immune$sample, "^CUF13I"))] <- 670
df_comp_immune$Carriage_density[which(str_detect(df_comp_immune$sample, "^CUF135"))] <- 33500
df_comp_immune$Carriage_density[which(str_detect(df_comp_immune$sample, "^CUF136"))] <- 67
df_comp_immune$Carriage_density[which(str_detect(df_comp_immune$sample, "^CUF13J"))] <- 26800
df_comp_immune$Carriage_density[which(str_detect(df_comp_immune$sample, "^CUF13K"))] <- NA
df_comp_immune$Carriage_density[which(str_detect(df_comp_immune$sample, "^CUF137"))] <- 10720
#df_comp_immune$Carriage_density[which(str_detect(df_comp_immune$sample, "^CUF12Y"))] <- 2345

abundance_immune <- df_comp_immune %>%
  pivot_longer(cols = c("T cells":"Dendritic cells"),
               names_to = "Cell Type",
               values_to = "Frequency") %>%
  ggplot(aes(x=`HIV_Status`,y=Frequency,fill=`Cell Type`))+
  geom_bar(stat = "identity",position = "fill")+
  scale_fill_manual(values=distinct_palette(n=16,pal = "brewerPlus",add = "lightgrey"))+
  theme_bw()+
  #geom_pwc()+
  labs(y="% of sample")+
  theme(legend.title = element_blank(),
        legend.position = "right",
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.text = element_text(size = 20,face = "bold"),
        axis.title = element_text(size = 20,face = "bold"),
        legend.text = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(angle = 45,hjust = 1))
abundance_immune


# Save the ggplot object----
ggsave(here("results", "abundance_immune.tiff"),
       plot = ((abundance_immune|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 15, height = 10, unit="in", dpi = 700,limitsize = F)


# ---- Save ----
save(all_merged_immune_final_labelled, file = "data/all_merged_immune_final_labelled.RData")

# save all merged_immune_final labelled as an H5seurat object
# ----- Subsetting epithelial cells -----
all_merged_epithelial <- subset(all_merged_subset_labelled_final,
                            idents = c("B cells","CD3+ T cells","Phagocytes",
                                       "Neutrophils","Dendritic cells"), invert = T)

all_merged_epithelial$Cell_CLusters <- paste0(all_merged_epithelial@active.ident)

view(all_merged_epithelial@meta.data)
# Dimplot all merged immune

DimPlot(all_merged_epithelial,
        reduction = "umap.harmony",
        #group.by = "Cell_CLusters"
        )


all_merged_epithelial<-IntegrateLayers(object=all_merged_epithelial,
                                   method=HarmonyIntegration,
                                   group.by = "sample",
                                   dims = 1:20,
                                   orig.reduction ="pca",
                                   new.reduction="harmony",
                                   verbose = TRUE)


all_merged_epithelial<-FindNeighbors(all_merged_epithelial, 
                                 reduction = "harmony",
                                 dims = 1:20)

all_merged_epithelial<-FindClusters(all_merged_epithelial, resolution = 0.5,
                                cluster.name = "harmony_clusters",group.singletons = TRUE,
                                method = "igraph")

all_merged_epithelial<-RunUMAP(all_merged_epithelial, reduction = "harmony",
                           dims = 1:20,reduction.name = "umap.harmony")

all_merged_epithelial<-RunTSNE(all_merged_epithelial, reduction = "harmony",
                           dims = 1:20,reduction.name = "tsne.harmony")

# Dimplot all merged epithelial

all_merged_epithelial_umap <- DimPlot(all_merged_epithelial,
                                  reduction = "umap.harmony",
                                  label = T, 
                                  repel = T,
                                  #pt.size = 2,
                                  #label.size = 15,
                                  #group.by = ("Cell_CLusters"),
                                  #split.by = ('HIV_Status')
)+
  theme_bw()+NoLegend()+
  labs(x="UMAP-1",y="UMAP-2")+
  theme(#axis.title = element_text(size = 20),
    #panel.border = element_rect(linewidth = 2),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank())+
  geom_hline(yintercept = 0,color="grey")+
  geom_vline(xintercept = 0,color="grey")

all_merged_epithelial_umap

# comvert dimplot to ggplot
all_merged_epithelial_umap <- as.ggplot(all_merged_epithelial_umap)

# save the saved converted dimplot 
ggsave(here("results", "all_merged_epithelial_umap.png"),
       plot = ((all_merged_epithelial_umap|plot_layout(ncol = 3, nrow = 2, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 7, unit="in", dpi = 700,limitsize = F)

FeaturePlot(all_merged_epithelial,
            reduction = "umap.harmony",
            features = c("KRT5","TUBA1A","MUC5AC","HELLS","KRT7"))

# find all_merged_epithelial markers
all_merged_epithelial <- JoinLayers(all_merged_epithelial)
all_merged_epithelial_markers <- FindAllMarkers(
  all_merged_epithelial,
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

# save all_merged_epithelial_markers-----
write.csv(all_merged_epithelial_markers,file = "data/all_merged_epithelial_markers.csv")

# Do heatmap of top 5 highly expressed genes from each epithelial cell cluster----
epithelial_cluster_heatmap <- DoHeatmap(subset(all_merged_epithelial,
                                                          downsample=500),
                                                   features = c("TMED7","LSM2","SLC26A2","TMTC3","SORD",
                                                                "CST1","FOS","HSD11B2","TNC","TP63",
                                                                "HLA-DPB1","MT-ND6","MTRNR2L12","UBTF","DDIT4",
                                                                "TFF3","SERPINB4","SERPINB3","SERPINF1","MT1X"),
                                                   slot = "data",
                                                   group.by = "ident",
                                                   group.bar = T,
                                                   angle = 90,
                                                   size = 4)+
  scale_fill_viridis()+
  #scale_fill_gradientn(colors = c("green", "white", "red"))+
  NoLegend()
epithelial_cluster_heatmap

epithelial_cluster_heatmap <- as.ggplot(epithelial_cluster_heatmap)
# Save the ggplot object
ggsave(here("results", "epithelial_cluster_heatmap.png"),
       plot = ((epithelial_cluster_heatmap|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 15, height = 10, unit="in", dpi = 700,limitsize = F)

VlnPlot(all_merged_epithelial,
        #group.by = c("Cell_CLusters"),
        features = c("BEST4","FOXJ1","SCEL","DEUP1",
                     "FOXI1","KRT15","MUC5AC","EGFR","FGFR3"),
        log = T,ncol = 5)


# ------save

save(all_merged_epithelial, file = "data/all_merged_epithelial.RData")
load("data/all_merged_epithelial.RData")

DimHeatmap(all_merged_epithelial,downsa)


# 
# Object ----
load('Thesis Figures/object.RData')

DimPlot(object)

# Immune Cell abundance
pt_immunineHIV <- table(Idents(object),object$HIV_Status)
pt_immunineHIV <- as.data.frame(pt_immunineHIV)


immune_cells_freq_plot <- pt_immunineHIV %>%
  ggplot(aes(Var2,Freq,fill=Var1))+
  geom_bar(stat = 'identity',position = 'fill')+
  labs(x='',
       y='Frequency of cells')+
  theme_bw()+
  theme(legend.position = 'right',
        panel.border = element_rect(),
        panel.grid = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_blank())
immune_cells_freq_plot



ggsave("Thesis Figures/immune_cells_freq_plot.pdf",
       plot = ((immune_cells_freq_plot|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 15, height = 7, unit="in", dpi = 700,limitsize = F)

