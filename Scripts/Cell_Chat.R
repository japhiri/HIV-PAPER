load("data/all_merged_subset_labelled_final.RData")

all_merged_subset_labelled_final$cell_clusters <- paste0(all_merged_subset_labelled_final@active.ident)
view(all_merged_subset_labelled_final@meta.data)


DimPlot(all_merged_subset_labelled_final,
        reduction = "umap.harmony",
        label = T)+NoLegend()
# Cell chat----
# Extracting input data from seurat v5 object----         
data.input<-LayerData(all_merged_subset_labelled_final, assay = "RNA",
                         layer = "data")

#labels <- Idents(filtered_seurat_subset_labelled)
labels <- all_merged_subset_labelled_final$cell_clusters


meta <- data.frame(group = labels, 
                   row.names = names(labels))

# Create CellChat object----
cellchat <- createCellChat(object = data.input,
                           group.by = "group")

cellchat <- addMeta(cellchat,meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat,
                     ident.use = "labels")
levels(cellchat@idents)

view(cellchat@meta$labels)
# Set the ligand-receptor interaction database----
CellChatDB <- CellChatDB.human

showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact")
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis----
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)

future::plan("multisession", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI
cellchat <- projectData(cellchat, PPI.human)

# Inference of cell-cell communication network
computeAveExpr(cellchat, 
               features = c("CXCL12","CXCR4","PDCD1","CD274","PDCD1LG2"), 
               type =  "truncatedMean", trim = 0.05)

# Compute the communication probability and infer cellular communication network----
cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups----
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract the inferred cellular communication network as a data frame----
df.net <- subsetCommunication(cellchat)
slot.name ="netP"
df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(3,4))
df.net <- subsetCommunication(cellchat, signaling = c("MK"))

# Infer the cell-cell communication at a signaling pathway level----
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network----
cellchat <- aggregateNet(cellchat)

view(cellchat@idents)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,3), xpd=TRUE)
netVisual_circle(cellchat@net$pval, vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = "Number of interactions")+
  facet_grid(.~HIV_Status)

netVisual_circle(cellchat@net$prob, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")



mat <- cellchat@net$weight
par(mfrow = c(2,4), xpd=TRUE)

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}


cellchat@netP$pathways


pathways.show <- c("MIF","MHC-I","CLEC","ALCAM","CD22","MHC-I","LCK","CD99",
                   "GALECTIN", "CD6","ADGRE5","CXCL","ANNEXIN"
                   ) 

vertex.receiver = seq(1,4)

netVisual_aggregate(cellchat, 
                    signaling = pathways.show,  
                    vertex.receiver = vertex.receiver)

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, 
                    signaling = pathways.show, 
                    layout = "circle")
# Choord diagram 
par(mfrow=c(1,2))
pdf(file ="cellchat.pdf", width = 20, height =16)

netVisual_aggregate(cellchat, 
                    signaling = pathways.show, 
                    layout = "chord")
dev.off()

# Heatmap
pathways.show <- c("MHC-I") 
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, #signaling = pathways.show, 
                  color.heatmap = "Reds",
                    comparison = C(1,2))


group.cellType <- c(rep("NKT cells", 4), rep("T cells", 4), rep("B cells", 4),
                    rep("Dendritic cells", 4), rep("Neutrophils", 4), rep("Monocytes/Macrophages", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, 
                     title.name = paste0(pathways.show, " signaling network"))

netAnalysis_contribution(cellchat, signaling = pathways.show)

# Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways ----
# Bubble plot
netVisual_bubble(cellchat, sources.use = 6, 
                 targets.use = c(2:11), remove.isolate = FALSE)

netVisual_bubble(cellchat, sources.use = 6, 
                 targets.use = c(5:11), 
                 signaling = c("MIF","GALECTIN"), remove.isolate = FALSE)

netVisual_chord_gene(cellchat, sources.use = 2, 
                     targets.use = c(2:11), 
                     lab.cex = 0.5,legend.pos.y = 40)

netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4, 5,6), targets.use = 2, legend.pos.x = 15)


netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4,5,6), targets.use = c(2:11), slot.name = "netP", legend.pos.x = 10)

# Plot the signaling gene expression distribution using violin/dot plot----

plotGeneExpression(cellchat, signaling="MIF")

plotGeneExpression(cellchat, signaling = "MIF", enriched.only = FALSE)

# Compute and visualize the network centrality scores----
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Visualize the dominant senders (sources) and receivers (targets) in a 2D space----
gg1 <- netAnalysis_signalingRole_scatter(cellchat)

gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("ADGRE5", "MIP"))

gg1 + gg2

# Identify signals contributing most to outgoing or incoming signaling of certain cell groups----

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2


# Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together----
selectK(cellchat, pattern ="outgoing")
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, 
                                          pattern = "outgoing", 
                                          k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing")

netAnalysis_dot(cellchat, pattern = "outgoing")
# Identify and visualize incoming communication pattern of target cells----
selectK(cellchat, pattern = "incoming")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, 
                                          pattern = "incoming", 
                                          k = nPatterns)

netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")


# Identify signaling groups based on their functional similarity----
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional",
                         umap.method = c("umap-learn","uwot"))




