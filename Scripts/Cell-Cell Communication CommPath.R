# Load required packages and dependencies-----
library(CommPath)
library(Matrix)
library(circlize)
library(ggplot2)
library(dplyr)
library(reshape2)
library(GSVA)

load('data/all_merged_subset_labelled_final.RData')

# Getting expression values from a seurat object
expression_values <- GetAssayData(all_merged_subset_labelled,
                                  layer = "data")

#expression_values <- as.data.frame(as.matrix(expression_values))

# Getting cell cluster information
cell_clusters <- all_merged_subset_labelled@meta.data$Clusters
# Creating a commPath object
commpath_object <- CommPath::createCommPath(expr.mat = expression_values,
                                            cell.info = cell_clusters,
                                            species = "hsapiens")
#Identification of marker ligands and receptors
commpath_object <- CommPath::findLRmarker(object = commpath_object,
                                          method = 'wilcox.test')
# Statistical identification of potential ligand-receptor (LR) associations
# To find significant LR pairs
commpath_object <- findLRpairs(object = commpath_object,
                               logFC.thre = 0,
                               p.thre = 0.01)
# Then we can visualize all interactions through a circos plot
# To show the counts of LR associations among all clusters
# Here we set the parameter "filter" as FALSE, which means that those LR interactions are identified only based on their expression profiles, not filtered by pathways in the receiver cells (as described in the later sections)
pdf('circosPlot.count.nonfiltered.pdf',height=6,width=6)
circosPlot(object = commpath_object, filter=FALSE)
dev.off()

# To show the overall interaction intensity of LR interactions among all clusters
pdf('circosPlot.intensity.nonfiltered.pdf',height=6,width=6)
circosPlot(object = commpath_object, plot="intensity", filter=FALSE)
dev.off()

# Pathway enrichment analysis
# To find pathways of which the genesets show overlap with the marker ligands and receptors
# CommPath provides pathway annotations from KEGG pathways, WikiPathways, reactome pathways, and GO terms
# Here we take the KEGG pathways as an example
commpath_object <- findLRpath(object = commpath_object, category = "kegg")

# To compute pathway activation score by the gsva algorithm or in an average manner
# For more information about the gsva algorithm, see the GSVA package (PMID23323831)
commpath_object <- scorePath(object = commpath_object, method = "gsva", min.size = 10, parallel.sz = 4)

# modified workaround ----
expr.mat <- commpath_object@data
expr.mat <- as.matrix(expr.mat)
# Ensure all values are numeric
expr.mat <- apply(expr.mat, 2, as.numeric)
rownames(expr.mat) <- rownames(original_matrix)
path.list <- commpath_object@pathway

library(GSVA)

# Create parameter object
params <- list(
  expr = expr.mat,
  gset.idx.list = path.list,
  method = "gsva",
  min.sz = 10,
  parallel.sz = 4
)
# Run GSVA
gsva_result <- do.call(gsva, params)


