# Regulatory Networks in Neutrophils using SCENIC-----
library(SCENIC)
library(GENIE3)
library(AUCell)
library(RcisTarget)
library(SCopeLoomR)

# Load seurat files
load("HIV PAPER/Data/Single_Cell_Data/Immune_cells.RData")
Immune_cells_sce <-as.SingleCellExperiment(Immune_cells)
# Get cluster information from sce object
cellInfo <- colData(sce)
# Getting the expression matrix from the sce object
exprMat <- counts(sce)
# Saving into a loom file
loom <- SCopeLoomR::build_loom("HIV PAPER/Data/Single_Cell_Data/Immune_cells.loom",
                               dgem = exprMat)
# Add cell annotations
loom <- SCENIC::add_cell_annotation(loom,cellInfo)
SCopeLoomR::close_loom(loom)
# Load .loom file
loomPath <- system.file(package = "SCENIC","HIV PAPER/Data/Single_Cell_Data/Immune_cells.loom")
loom <- SCopeLoomR::open_loom("HIV PAPER/Data/Single_Cell_Data/Immune_cells.loom")
exprMat <- SCopeLoomR::get_dgem(loom)
cellInfo <- SCopeLoomR::get_cell_annotation(loom)
SCopeLoomR::close_loom(loom)
# Initialize settings
library(SCENIC)
library(RcisTarget)
data("motifAnnotations_hgnc_v9",package = "RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
ScenicOptions <- SCENIC::initializeScenic(org = 'hgnc',dbDir = "HIV PAPER/Data/Single_Cell_Data/cisTarget_databases",nCores = 1)
saveRDS(ScenicOptions,file = "HIV PAPER/Data/Single_Cell_Data/ScenicOptions.rds")

# Coexpression network
genesKept <- SCENIC::geneFiltering(exprMat,ScenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
SCENIC::runCorrelation(exprMat_filtered, ScenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, ScenicOptions)
exprMat <- getNormalizedMatrix(ScenicOptions)


# Build and score the GRN
exprMat_log <- log2(exprMat+1)
ScenicOptions@settings$dbs <- ScenicOptions@settings$dbs["10kb"]
#ScenicOptions <- SCENIC::runSCE (ScenicOptions)
ScenicOptions <- SCENIC::runSCENIC_1_coexNetwork2modules(ScenicOptions)
ScenicOptions <- SCENIC::runSCENIC_2createRegulons(ScenicOptions, coexMethod=c("top5perTarget"))
ScenicOptions <- SCENIC::runSCENIC_3_scoreCells(ScenicOptions, exprMat_log)

# Optional: Binarize activity
aucellApp <- plotTsne_AUCellApp(ScenicOptions, exprMat_log)
savedSelection <- shinny::runAPP(aucellApp)
newThresholds <- savedSelection$thresholds
ScenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(ScenicOptions, "aucell_thresholds"))
ScenicOption <- SCENIC::runSCENIC_4_aucell_binarize(ScenicOptions)
tsneAUC(ScenicOptions,aucType="AUC")

# Export:
SaveRDS(cellInfo, file=getDatasetInfo(ScenicOptions, "cellInfo"))
export2loom(ScenicOptions, exprMat)

# To save the current status, or any changes in settings, save the object agains
saveRDS(ScenicOptions,file = "HIV PAPER/Data/Single_Cell_Data/ScenicOptions.Rds")

# exploring output
motifEnrichment_selfMotifs_wGenes <- loadInt(ScenicOptions,"motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedRFs=="Sox8"]
viewMotifs(tableSubset)




