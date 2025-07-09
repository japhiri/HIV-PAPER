# load seurat object -----
load("data/all_merged_subset_labelled_final.RData")
load('Thesis Figures/object.RData')

all_merged_subset_labelled_final$Cell_clusters <- paste0(all_merged_subset_labelled_final@active.ident)
view(all_merged_subset_labelled_final@meta.data)

object$Cell_clusters <- paste0(object@active.ident)
#DimPlot
DimPlot(object,
        reduction = "umap.harmony",
        #group.by = "Cell_clusters",
        label = T)+NoLegend()
# Perform DE analysis within the same cell type across conditions----
all_merged_subset_labelled_final$celltype.HIV.status <- paste(all_merged_subset_labelled_final$Cell_clusters,
                                                              all_merged_subset_labelled_final$HIV_Status, sep = "_")

object$celltype.HIV.status <- paste(object$Cell_clusters,
                                    object$HIV_Status, sep = "_")

Idents(all_merged_subset_labelled_final) <- "celltype.HIV.status"
Idents(object) <- "celltype.HIV.status"
T_cells_de <- FindMarkers(all_merged_subset_labelled_final,
                          ident.1 = "CD3+ T cells_HIV+ ART<3 Months",
                          ident.2 = "CD3+ T cells_HIV-")


T_cells_de <- FindMarkers(object,
                          ident.1 = "T cells_HIV+ ART<3 Months",
                          ident.2 = "T cells_HIV-")
head(T_cells_de,n=10)

# Perform DE analysis after pseudobulking-----

pseudo_exp <- AggregateExpression(all_merged_subset_labelled_final,
                                  assays = "RNA",
                                  slot = "counts",
                                  return.seurat = T,
                                  group.by = c("sample",
                                               "HIV_Status",
                                               "Cell_clusters"))

pseudo_exp <- AggregateExpression(object,
                                  assays = "RNA",
                                  slot = "counts",
                                  return.seurat = T,
                                  group.by = c("sample",
                                               "HIV_Status",
                                               "Cell_clusters"))




tail(Cells(pseudo_exp))

pseudo_exp$Cell_clusters <- sapply(strsplit(Cells(pseudo_exp),
                                            split = "_"),"[",3)

pseudo_exp$sample <- sapply(strsplit(Cells(pseudo_exp),
                                     split = "_"),
                            "[",1)

pseudo_exp$HIV_Status <- sapply(strsplit(Cells(pseudo_exp),
                                         split = "_"),
                                "[",2)

pseudo_exp$Celltype.condition <- paste(pseudo_exp$Cell_clusters,
                                       pseudo_exp$HIV_Status,
                                       sep = "_")



Idents(pseudo_exp) <- "Celltype.condition"

# save pseudo_exp object-----
save(pseudo_exp,file = 'data/pseudo_exp.RData')
load('data/pseudo_exp.RData')
#T cell DEG----
#3M vs HIV-
bulk.T.cells.de_3MvsNeg <- FindMarkers(object = pseudo_exp,
                               ident.1 = "T cells_HIV+ ART<3 Months",
                               ident.2 = "T cells_HIV-",
                               min.pct = 0.25,
                               test.use = "DESeq2")
write.csv(bulk.T.cells.de_3MvsNeg,file = 'data/bulk.T.cells.de_3MvsNeg.csv')

#bulk.T.cells.de_3MvsNeg %>%
 #mutate(pct.diff=pct.1-pct.2)
head(bulk.T.cells.de_3MvsNeg,n=50)
view(bulk.T.cells.de_3MvsNeg)

#1Y vs HIV-
bulk.T.cells.de_1YvsNeg <- FindMarkers(object = pseudo_exp,
                                       ident.1 = "T cells_HIV+ ART>1 Year",
                                       ident.2 = "T cells_HIV-",
                                       min.pct = 0.25,
                                       test.use = "DESeq2")
write.csv(bulk.T.cells.de_1YvsNeg,file = 'data/bulk.T.cells.de_1YvsNeg.csv')
head(bulk.T.cells.de_1YvsNeg,n=50)


#3M vs 1Y
bulk.T.cells.de_3Mvs1Y <- FindMarkers(object = pseudo_exp,
                                      ident.1 = "T cells_HIV+ ART<3 Months",
                                      ident.2 = "T cells_HIV+ ART>1 Year",
                                      min.pct = 0.25,
                                      test.use = "DESeq2")
write.csv(bulk.T.cells.de_3Mvs1Y,file = 'data/bulk.T.cells.de_3Mvs1Y.csv')
head(bulk.T.cells.de_3Mvs1Y,n=50)


#Plot a Volcano plot 3M vs HIV- ----

bulk.T.cells.de_3MvsNeg.df <- as.data.frame(bulk.T.cells.de_3MvsNeg)
bulk.T.cells.de_3MvsNeg.df$genes <-rownames(bulk.T.cells.de_3MvsNeg)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.T.cells.de_3MvsNeg.df$DiffExpressed <- "NO"

bulk.T.cells.de_3MvsNeg.df$DiffExpressed[bulk.T.cells.de_3MvsNeg.df$avg_log2FC >=1.5 &
                                           bulk.T.cells.de_3MvsNeg.df$p_val<0.05] <- "UP"

bulk.T.cells.de_3MvsNeg.df$DiffExpressed[bulk.T.cells.de_3MvsNeg.df$avg_log2FC < -1.5 &
                                           bulk.T.cells.de_3MvsNeg.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.T.cells.de_3MvsNeg.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.T.cells.de_3MvsNeg.df$diffLabel[bulk.T.cells.de_3MvsNeg.df$DiffExpressed !="NO"] <-
  bulk.T.cells.de_3MvsNeg.df$genes[bulk.T.cells.de_3MvsNeg.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Tcells3MvsNeg <- bulk.T.cells.de_3MvsNeg.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "T cells",
       subtitle = "HIV+ ART<3 Months vs HIV-")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        panel.border = element_rect(),
        panel.grid = element_line(linewidth = 0.05),
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Tcells3MvsNeg <- as.ggplot(Tcells3MvsNeg)
Tcells3MvsNeg



# Diff genes bargraph 3MvsNeg----
bulk.T.cells.de_3MvsNeg.df_top20 <- bulk.T.cells.de_3MvsNeg.df %>%
  filter(abs(avg_log2FC)>=1.5,p_val<0.05,genes!="NA") 
bulk.T.cells.de_3MvsNeg.df_top20 <- bulk.T.cells.de_3MvsNeg.df_top20[order(abs(bulk.T.cells.de_3MvsNeg.df_top20$avg_log2FC)),]
bulk.T.cells.de_3MvsNeg.df_top20 <- head(bulk.T.cells.de_3MvsNeg.df_top20,50)


bulk.T.cells.de_3MvsNeg.df_top20$diffLabel <- reorder(bulk.T.cells.de_3MvsNeg.df_top20$diffLabel,
                                                      bulk.T.cells.de_3MvsNeg.df_top20$avg_log2FC)



T_cells_3MvsNegcol <- bulk.T.cells.de_3MvsNeg.df_top20 %>%
  ggplot(aes(diffLabel,avg_log2FC, fill=avg_log2FC
             ))+
  geom_col()+
  coord_flip()+
  theme_pubr()+
  #scale_y_continuous(limits = c(-8,8))+
  scale_fill_gradient2(low = "red",mid = "blue",high = "red")+
  labs(x="",y="Average log2FoldChange",
       title = "T cells",
       subtitle = "HIV+ ART<3 Months vs HIV-")+
  guides(fill=guide_legend(title="Average log2FC"))+
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 20),
        plot.subtitle = element_text(hjust = 0.5,face = 'bold',size = 15,color='grey'))
T_cells_3MvsNegcol


# save the saved converted dimplot 
ggsave("results/T_cells_3MvsNegcol.pdf",
       plot = ((T_cells_3MvsNegcol|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 25, height = 12, unit="in", dpi = 700,limitsize = F)



#Plot a Volcano plot 1Y vs HIV- ----

bulk.T.cells.de_1YvsNeg.df <- as.data.frame(bulk.T.cells.de_1YvsNeg)
bulk.T.cells.de_1YvsNeg.df$genes <-rownames(bulk.T.cells.de_1YvsNeg)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.T.cells.de_1YvsNeg.df$DiffExpressed <- "NO"

bulk.T.cells.de_1YvsNeg.df$DiffExpressed[bulk.T.cells.de_1YvsNeg.df$avg_log2FC >=1.5 &
                                           bulk.T.cells.de_1YvsNeg.df$p_val<0.05] <- "UP"

bulk.T.cells.de_1YvsNeg.df$DiffExpressed[bulk.T.cells.de_1YvsNeg.df$avg_log2FC < -1.5 &
                                           bulk.T.cells.de_1YvsNeg.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.T.cells.de_1YvsNeg.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.T.cells.de_1YvsNeg.df$diffLabel[bulk.T.cells.de_1YvsNeg.df$DiffExpressed !="NO"] <-
  bulk.T.cells.de_1YvsNeg.df$genes[bulk.T.cells.de_1YvsNeg.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Tcells1YvsNeg <- bulk.T.cells.de_1YvsNeg.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "T cells",
       subtitle = "HIV+ ART>1 Year vs HIV-")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Tcells1YvsNeg <- as.ggplot(Tcells1YvsNeg)
Tcells1YvsNeg


# Diff genes bargraph 1YvsNeg----
bulk.T.cells.de_1YvsNeg.df_top20 <- bulk.T.cells.de_1YvsNeg.df %>%
  filter(abs(avg_log2FC)>=1.5,p_val<0.05) 
bulk.T.cells.de_1YvsNeg.df_top20 <- bulk.T.cells.de_1YvsNeg.df_top20[order(bulk.T.cells.de_3MvsNeg.df_top20$avg_log2FC, decreasing = T),]

bulk.T.cells.de_1YvsNeg.df_top20$diffLabel <- reorder(bulk.T.cells.de_1YvsNeg.df_top20$diffLabel,
                                                      bulk.T.cells.de_1YvsNeg.df_top20$avg_log2FC)



T_cells_1YvsNegcol <- bulk.T.cells.de_1YvsNeg.df_top20 %>%
  ggplot(aes(diffLabel,avg_log2FC, fill=avg_log2FC
  ))+
  geom_col()+
  coord_flip()+
  theme_pubr()+
  #scale_y_continuous(limits = c(-8,8))+
  scale_fill_gradient2(low = "red",mid = "blue",high = "red")+
  labs(x="",y="Average log2FoldChange",
       title = "T cells",
       subtitle = "HIV+ ART>1 Year vs HIV-")+
  guides(fill=guide_legend(title="Average log2FC"))+
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 20),
        plot.subtitle = element_text(hjust = 0.5,face = 'bold',size = 15,color='grey'))
T_cells_1YvsNegcol

#Plot a Volcano plot 3M vs 1Y ----

bulk.T.cells.de_3Mvs1Y.df <- as.data.frame(bulk.T.cells.de_3Mvs1Y)
bulk.T.cells.de_3Mvs1Y.df$genes <-rownames(bulk.T.cells.de_3Mvs1Y)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.T.cells.de_3Mvs1Y.df$DiffExpressed <- "NO"

bulk.T.cells.de_3Mvs1Y.df$DiffExpressed[bulk.T.cells.de_3Mvs1Y.df$avg_log2FC >=1.5 &
                                          bulk.T.cells.de_3Mvs1Y.df$p_val<0.05] <- "UP"

bulk.T.cells.de_3Mvs1Y.df$DiffExpressed[bulk.T.cells.de_3Mvs1Y.df$avg_log2FC < -1.5 &
                                          bulk.T.cells.de_3Mvs1Y.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.T.cells.de_3Mvs1Y.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.T.cells.de_3Mvs1Y.df$diffLabel[bulk.T.cells.de_3Mvs1Y.df$DiffExpressed !="NO"] <-
  bulk.T.cells.de_3Mvs1Y.df$genes[bulk.T.cells.de_3Mvs1Y.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Tcells3Mvs1Y <- bulk.T.cells.de_3Mvs1Y.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "T cells",
       subtitle = "HIV+ ART<3 months vs HIV+ ART>1 Year")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Tcells3Mvs1Y <- as.ggplot(Tcells3Mvs1Y)
Tcells3Mvs1Y

Tcells3MvsNeg+Tcells1YvsNeg+Tcells3Mvs1Y

# save the saved converted dimplot 
ggsave("results/TcellsVolcanoPlots.pdf",
       plot = ((Tcells3MvsNeg|Tcells1YvsNeg|Tcells3Mvs1Y|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 8, unit="in", dpi = 700,limitsize = F)


# Diff genes bargraph 3Mvs1Y----
bulk.T.cells.de_3Mvs1Y.df_top20 <- bulk.T.cells.de_3Mvs1Y.df %>%
  filter(abs(avg_log2FC)>=1.5,p_val<0.05) 
bulk.T.cells.de_3Mvs1Y.df_top20 <- bulk.T.cells.de_3Mvs1Y.df_top20[order(bulk.T.cells.de_3MvsNeg.df_top20$avg_log2FC, decreasing = T),]

bulk.T.cells.de_3Mvs1Y.df_top20$diffLabel <- reorder(bulk.T.cells.de_3Mvs1Y.df_top20$diffLabel,
                                                      bulk.T.cells.de_3Mvs1Y.df_top20$avg_log2FC)



T_cells_3Mvs1Ycol <- bulk.T.cells.de_3Mvs1Y.df_top20 %>%
  ggplot(aes(diffLabel,avg_log2FC, fill=avg_log2FC
  ))+
  geom_col()+
  coord_flip()+
  theme_pubr()+
  #scale_y_continuous(limits = c(-8,8))+
  scale_fill_gradient2(low = "red",mid = "blue",high = "red")+
  labs(x="",y="Average log2FoldChange",
       title = "T cells",
       subtitle = "HIV+ ART<3 Months vs HIV+ ART>1 Year")+
  guides(fill=guide_legend(title="Average log2FC"))+
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 20),
        plot.subtitle = element_text(hjust = 0.5,face = 'bold',size = 15,color='grey'))
T_cells_3Mvs1Ycol

# save the saved converted dimplot ----
ggsave("results/TcelldiffbarPlots.pdf",
       plot = ((T_cells_3MvsNegcol|T_cells_1YvsNegcol|T_cells_3Mvs1Ycol|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 10, unit="in", dpi = 700,limitsize = F)

# Plot violins of significant genes----
VlnPlot(subset(all_merged_subset_labelled_final,downsample=200),
        features = c("CCL4","NKG7","KLRB1","LAG3","CD74","GZMK","GZMH"),
        idents = c("CD3+ T cells_HIV-","CD3+ T cells_HIV+ ART<3 Months",
                   group.by="HIV_Status"))


DoHeatmap(subset(all_merged_subset_labelled_final,downsample=100),
          features = c("CCL4","NKG7","KLRB1","LAG3","CD74","GZMK","GZMH",
                       "HLA-DRB1","HLA-A ","TRAV29DV5","CTTNBP2","HDC","GBP5",
                       "GBP4","BMP2K","KIT","STAT1","MX1","IRF1"),
          group.by = "HIV_Status",
          angle = 90)+
  scale_fill_viridis()



#B cell DEG----
#3M vs HIV-
bulk.B.cells.de_3MvsNeg <- FindMarkers(object = pseudo_exp,
                                       ident.1 = "B cells_HIV+ ART<3 Months",
                                       ident.2 = "B cells_HIV-",
                                       test.use = "DESeq2")
head(bulk.B.cells.de_3MvsNeg,n=50)

#1Y vs HIV-
bulk.B.cells.de_1YvsNeg <- FindMarkers(object = pseudo_exp,
                                       ident.1 = "B cells_HIV+ ART>1 Year",
                                       ident.2 = "B cells_HIV-",
                                       test.use = "DESeq2")
head(bulk.B.cells.de_1YvsNeg,n=50)



bulk.B.cells.de_3Mvs1Y <- FindMarkers(object = pseudo_exp,
                                      ident.1 = "B cells_HIV+ ART<3 Months",
                                      ident.2 = "B cells_HIV+ ART>1 Year",
                                      test.use = "DESeq2")
head(bulk.B.cells.de_3Mvs1Y,n=50)


#Plot a Volcano plot 3M vs HIV- ----

bulk.B.cells.de_3MvsNeg.df <- as.data.frame(bulk.B.cells.de_3MvsNeg)
bulk.B.cells.de_3MvsNeg.df$genes <-rownames(bulk.B.cells.de_3MvsNeg)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.B.cells.de_3MvsNeg.df$DiffExpressed <- "NO"

bulk.B.cells.de_3MvsNeg.df$DiffExpressed[bulk.B.cells.de_3MvsNeg.df$avg_log2FC >=1.5 &
                                           bulk.B.cells.de_3MvsNeg.df$p_val<0.05] <- "UP"

bulk.B.cells.de_3MvsNeg.df$DiffExpressed[bulk.B.cells.de_3MvsNeg.df$avg_log2FC < -1.5 &
                                           bulk.B.cells.de_3MvsNeg.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.B.cells.de_3MvsNeg.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.B.cells.de_3MvsNeg.df$diffLabel[bulk.B.cells.de_3MvsNeg.df$DiffExpressed !="NO"] <-
  bulk.B.cells.de_3MvsNeg.df$genes[bulk.B.cells.de_3MvsNeg.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Bcells3MvsNeg <- bulk.B.cells.de_3MvsNeg.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "B cells",
       subtitle = "HIV+ ART<3 Months vs HIV-")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Bcells3MvsNeg <- as.ggplot(Bcells3MvsNeg)
Bcells3MvsNeg




#Plot a Volcano plot 1Y vs HIV- ----

bulk.B.cells.de_1YvsNeg.df <- as.data.frame(bulk.B.cells.de_1YvsNeg)
bulk.B.cells.de_1YvsNeg.df$genes <-rownames(bulk.B.cells.de_1YvsNeg)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.B.cells.de_1YvsNeg.df$DiffExpressed <- "NO"

bulk.B.cells.de_1YvsNeg.df$DiffExpressed[bulk.B.cells.de_1YvsNeg.df$avg_log2FC >=1.5 &
                                           bulk.B.cells.de_1YvsNeg.df$p_val<0.05] <- "UP"

bulk.B.cells.de_1YvsNeg.df$DiffExpressed[bulk.B.cells.de_1YvsNeg.df$avg_log2FC < -1.5 &
                                           bulk.B.cells.de_1YvsNeg.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.B.cells.de_1YvsNeg.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.B.cells.de_1YvsNeg.df$diffLabel[bulk.B.cells.de_1YvsNeg.df$DiffExpressed !="NO"] <-
  bulk.B.cells.de_1YvsNeg.df$genes[bulk.B.cells.de_1YvsNeg.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Bcells1YvsNeg <- bulk.B.cells.de_1YvsNeg.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "B cells",
       subtitle = "HIV+ ART>1 Year vs HIV-")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Bcells1YvsNeg <- as.ggplot(Bcells1YvsNeg)
Bcells1YvsNeg

#Plot a Volcano plot 3M vs 1Y ----

bulk.B.cells.de_3Mvs1Y.df <- as.data.frame(bulk.B.cells.de_3Mvs1Y)
bulk.B.cells.de_3Mvs1Y.df$genes <-rownames(bulk.B.cells.de_3Mvs1Y)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.B.cells.de_3Mvs1Y.df$DiffExpressed <- "NO"

bulk.B.cells.de_3Mvs1Y.df$DiffExpressed[bulk.B.cells.de_3Mvs1Y.df$avg_log2FC >=1.5 &
                                          bulk.B.cells.de_3Mvs1Y.df$p_val<0.05] <- "UP"

bulk.B.cells.de_3Mvs1Y.df$DiffExpressed[bulk.B.cells.de_3Mvs1Y.df$avg_log2FC < -1.5 &
                                          bulk.B.cells.de_3Mvs1Y.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.B.cells.de_3Mvs1Y.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.B.cells.de_3Mvs1Y.df$diffLabel[bulk.B.cells.de_3Mvs1Y.df$DiffExpressed !="NO"] <-
  bulk.B.cells.de_3Mvs1Y.df$genes[bulk.B.cells.de_3Mvs1Y.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Bcells3Mvs1Y <- bulk.B.cells.de_3Mvs1Y.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "B cells",
       subtitle = "HIV+ ART<3 months vs HIV+ ART>1 Year")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Bcells3Mvs1Y <- as.ggplot(Bcells3Mvs1Y)
Bcells3Mvs1Y


Bcells3MvsNeg+Bcells1YvsNeg+Bcells3Mvs1Y

# save the saved converted dimplot 
ggsave("results/BcellsVolcanoPlots.pdf",
       plot = ((Bcells3MvsNeg|Bcells1YvsNeg|Bcells3Mvs1Y|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 8, unit="in", dpi = 700,limitsize = F)


#Neutrophils cell DEG----
#3M vs HIV-
bulk.Neutrophils.cells.de_3MvsNeg <- FindMarkers(object = pseudo_exp,
                                       ident.1 = "Neutrophils_HIV+ ART<3 Months",
                                       ident.2 = "Neutrophils_HIV-",
                                       test.use = "DESeq2")
head(bulk.Neutrophils.cells.de_3MvsNeg,n=50)

#1Y vs HIV-
bulk.Neutrophils.cells.de_1YvsNeg <- FindMarkers(object = pseudo_exp,
                                       ident.1 = "Neutrophils_HIV+ ART>1 Year",
                                       ident.2 = "Neutrophils_HIV-",
                                       test.use = "DESeq2")
head(bulk.Neutrophils.cells.de_1YvsNeg,n=50)



bulk.Neutrophils.cells.de_3Mvs1Y <- FindMarkers(object = pseudo_exp,
                                      ident.1 = "Neutrophils_HIV+ ART<3 Months",
                                      ident.2 = "Neutrophils_HIV+ ART>1 Year",
                                      test.use = "DESeq2")
head(bulk.Neutrophils.cells.de_3Mvs1Y,n=50)


#Plot a Volcano plot 3M vs HIV- ----

bulk.Neutrophils.cells.de_3MvsNeg.df <- as.data.frame(bulk.Neutrophils.cells.de_3MvsNeg)
bulk.Neutrophils.cells.de_3MvsNeg.df$genes <-rownames(bulk.Neutrophils.cells.de_3MvsNeg)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Neutrophils.cells.de_3MvsNeg.df$DiffExpressed <- "NO"

bulk.Neutrophils.cells.de_3MvsNeg.df$DiffExpressed[bulk.Neutrophils.cells.de_3MvsNeg.df$avg_log2FC >=1.5 &
                                                     bulk.Neutrophils.cells.de_3MvsNeg.df$p_val<0.05] <- "UP"

bulk.Neutrophils.cells.de_3MvsNeg.df$DiffExpressed[bulk.Neutrophils.cells.de_3MvsNeg.df$avg_log2FC < -1.5 &
                                                     bulk.Neutrophils.cells.de_3MvsNeg.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Neutrophils.cells.de_3MvsNeg.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Neutrophils.cells.de_3MvsNeg.df$diffLabel[bulk.Neutrophils.cells.de_3MvsNeg.df$DiffExpressed !="NO"] <-
  bulk.Neutrophils.cells.de_3MvsNeg.df$genes[bulk.Neutrophils.cells.de_3MvsNeg.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Neutrophilscells3MvsNeg <- bulk.Neutrophils.cells.de_3MvsNeg.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 20)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Neutrophils",
       subtitle = "HIV+ ART<3 Months vs HIV-")+
  #scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Neutrophilscells3MvsNeg
Neutrophilscells3MvsNeg <- as.ggplot(Neutrophilscells3MvsNeg)

# Diff genes bargraph 3MvsNeg----
bulk.Neutrophils.cells.de_3MvsNeg.df_top20 <- bulk.Neutrophils.cells.de_3MvsNeg.df %>%
  filter(abs(avg_log2FC)>=1.5, p_val<0.05,genes!="NA") 
bulk.Neutrophils.cells.de_3MvsNeg.df_top20 <- bulk.Neutrophils.cells.de_3MvsNeg.df_top20[order(abs(bulk.Neutrophils.cells.de_3MvsNeg.df_top20$avg_log2FC)),]
#bulk.Neutrophils.cells.de_3MvsNeg.df_top20 <- head(bulk.Neutrophils.cells.de_3MvsNeg.df_top20,100)


bulk.Neutrophils.cells.de_3MvsNeg.df_top20$diffLabel <- reorder(bulk.Neutrophils.cells.de_3MvsNeg.df_top20$diffLabel,
                                                                bulk.Neutrophils.cells.de_3MvsNeg.df_top20$avg_log2FC)


Neutrophils_3MvsNegcol <- bulk.Neutrophils.cells.de_3MvsNeg.df_top20 %>%
  ggplot(aes(diffLabel,avg_log2FC, fill=avg_log2FC
  ))+
  geom_col()+
  coord_flip()+
  theme_pubr()+
  #scale_y_continuous(limits = c(-8,8))+
  scale_fill_gradient2(low = "red",mid = "blue",high = "red")+
  labs(x="",y="Average log2FoldChange",
       title = "Neutrophils",
       subtitle = "HIV+ ART<3 Months vs HIV-")+
  guides(fill=guide_legend(title="Average log2FC"))+
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 20),
        plot.subtitle = element_text(hjust = 0.5,face = 'bold',size = 15,color='grey'))
Neutrophils_3MvsNegcol


#Plot a Volcano plot 1Y vs HIV- ----

bulk.Neutrophils.cells.de_1YvsNeg.df <- as.data.frame(bulk.Neutrophils.cells.de_1YvsNeg)
bulk.Neutrophils.cells.de_1YvsNeg.df$genes <-rownames(bulk.Neutrophils.cells.de_1YvsNeg)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Neutrophils.cells.de_1YvsNeg.df$DiffExpressed <- "NO"

bulk.Neutrophils.cells.de_1YvsNeg.df$DiffExpressed[bulk.Neutrophils.cells.de_1YvsNeg.df$avg_log2FC >=1.5 &
                                                     bulk.Neutrophils.cells.de_1YvsNeg.df$p_val<0.05] <- "UP"

bulk.Neutrophils.cells.de_1YvsNeg.df$DiffExpressed[bulk.Neutrophils.cells.de_1YvsNeg.df$avg_log2FC < -1.5 &
                                                     bulk.Neutrophils.cells.de_1YvsNeg.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Neutrophils.cells.de_1YvsNeg.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Neutrophils.cells.de_1YvsNeg.df$diffLabel[bulk.Neutrophils.cells.de_1YvsNeg.df$DiffExpressed !="NO"] <-
  bulk.Neutrophils.cells.de_1YvsNeg.df$genes[bulk.Neutrophils.cells.de_1YvsNeg.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Neutrophilscells1YvsNeg <- bulk.Neutrophils.cells.de_1YvsNeg.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Neutrophils",
       subtitle = "HIV+ ART>1 Year vs HIV-")+
  #scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Neutrophilscells1YvsNeg
Neutrophilscells1YvsNeg <- as.ggplot(Neutrophilscells1YvsNeg)

# Diff genes bargraph 1YvsNeg----
bulk.Neutrophils.cells.de_1YvsNeg.df_top20 <- bulk.Neutrophils.cells.de_1YvsNeg.df %>%
  filter(abs(avg_log2FC)>1.5,p_val<0.05,genes!="NA") 
bulk.Neutrophils.cells.de_1YvsNeg.df_top20 <- bulk.Neutrophils.cells.de_1YvsNeg.df_top20[order(abs(bulk.Neutrophils.cells.de_1YvsNeg.df_top20$avg_log2FC)),]
bulk.Neutrophils.cells.de_1YvsNeg.df_top20 <- head(bulk.Neutrophils.cells.de_1YvsNeg.df_top20,50)


bulk.Neutrophils.cells.de_1YvsNeg.df_top20$diffLabel <- reorder(bulk.Neutrophils.cells.de_1YvsNeg.df_top20$diffLabel,
                                                                bulk.Neutrophils.cells.de_1YvsNeg.df_top20$avg_log2FC)


Neutrophils_1YvsNegcol <- bulk.Neutrophils.cells.de_1YvsNeg.df_top20 %>%
  ggplot(aes(diffLabel,avg_log2FC, fill=avg_log2FC
  ))+
  geom_col()+
  coord_flip()+
  theme_pubr()+
  #scale_y_continuous(limits = c(-8,8))+
  scale_fill_gradient2(low = "red",mid = "blue",high = "red")+
  labs(x="",y="Average log2FoldChange",
       title = "Neutrophils",
       subtitle = "HIV+ ART>1 Year vs HIV-")+
  guides(fill=guide_legend(title="Average log2FC"))+
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 20),
        plot.subtitle = element_text(hjust = 0.5,face = 'bold',size = 15,color='grey'))
Neutrophils_1YvsNegcol


#Plot a Volcano plot 3M vs 1Y ----

bulk.Neutrophils.cells.de_3Mvs1Y.df <- as.data.frame(bulk.Neutrophils.cells.de_3Mvs1Y)
bulk.Neutrophils.cells.de_3Mvs1Y.df$genes <-rownames(bulk.Neutrophils.cells.de_3Mvs1Y)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Neutrophils.cells.de_3Mvs1Y.df$DiffExpressed <- "NO"

bulk.Neutrophils.cells.de_3Mvs1Y.df$DiffExpressed[bulk.Neutrophils.cells.de_3Mvs1Y.df$avg_log2FC >=1.5 &
                                                    bulk.Neutrophils.cells.de_3Mvs1Y.df$p_val<0.05] <- "UP"

bulk.Neutrophils.cells.de_3Mvs1Y.df$DiffExpressed[bulk.Neutrophils.cells.de_3Mvs1Y.df$avg_log2FC <= -1.5 &
                                                    bulk.Neutrophils.cells.de_3Mvs1Y.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Neutrophils.cells.de_3Mvs1Y.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Neutrophils.cells.de_3Mvs1Y.df$diffLabel[bulk.Neutrophils.cells.de_3Mvs1Y.df$DiffExpressed !="NO"] <-
  bulk.Neutrophils.cells.de_3Mvs1Y.df$genes[bulk.Neutrophils.cells.de_3Mvs1Y.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Neutrophilscells3Mvs1Y <- bulk.Neutrophils.cells.de_3Mvs1Y.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Neutrophils",
       subtitle = "HIV+ ART<3 months vs HIV+ ART>1 Year")+
  #scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Neutrophilscells3Mvs1Y
Neutrophilscells3Mvs1Y <- as.ggplot(Neutrophilscells3Mvs1Y)


# Diff genes bargraph 3Mvs1Y----
bulk.Neutrophils.cells.de_3Mvs1Y.df_top20 <- bulk.Neutrophils.cells.de_3Mvs1Y.df %>%
  filter(abs(avg_log2FC)>1.5,p_val<0.05,genes!="NA") 
bulk.Neutrophils.cells.de_3Mvs1Y.df_top20 <- bulk.Neutrophils.cells.de_3Mvs1Y.df_top20[order(abs(bulk.Neutrophils.cells.de_3Mvs1Y.df_top20$avg_log2FC)),]
#bulk.Neutrophils.cells.de_3Mvs1Y.df_top20 <- head(bulk.Neutrophils.cells.de_3Mvs1Y.df_top20,100)


bulk.Neutrophils.cells.de_3Mvs1Y.df_top20$diffLabel <- reorder(bulk.Neutrophils.cells.de_3Mvs1Y.df_top20$diffLabel,
                                                               bulk.Neutrophils.cells.de_3Mvs1Y.df_top20$avg_log2FC)


Neutrophils_3Mvs1Ycol <- bulk.Neutrophils.cells.de_3Mvs1Y.df_top20 %>%
  ggplot(aes(diffLabel,avg_log2FC, fill=avg_log2FC
  ))+
  geom_col()+
  coord_flip()+
  theme_pubr()+
  #scale_y_continuous(limits = c(-8,8))+
  scale_fill_gradient2(low = "red",mid = "blue",high = "red")+
  labs(x="",y="Average log2FoldChange",
       title = "Neutrophils",
       subtitle = "HIV+ ART<3 Months vs HIV+ ART>1 Year")+
  guides(fill=guide_legend(title="Average log2FC"))+
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 20),
        plot.subtitle = element_text(hjust = 0.5,face = 'bold',size = 15,color='grey'))
Neutrophils_3Mvs1Ycol

Neutrophilscells3MvsNeg+Neutrophilscells1YvsNeg+Neutrophilscells3Mvs1Y

# save the saved converted dimplot 
ggsave("results/NeutrophilsVolcanoPlots.pdf",
       plot = ((Neutrophilscells3MvsNeg|Neutrophilscells1YvsNeg|Neutrophilscells3Mvs1Y|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 8, unit="in", dpi = 700,limitsize = F)

# save the saved converted dimplot 
ggsave("results/NeutrophilsbarPlots.pdf",
       plot = ((Neutrophils_3MvsNegcol|Neutrophils_1YvsNegcol|Neutrophils_3Mvs1Ycol|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 8, unit="in", dpi = 700,limitsize = F)








#Goblet cell DEG----
#3M vs HIV-
bulk.Goblet.cells.de_3MvsNeg <- FindMarkers(object = pseudo_exp,
                                                 ident.1 = "Goblet cells_HIV+ ART<3 Months",
                                                 ident.2 = "Goblet cells_HIV-",
                                                 test.use = "DESeq2")
head(bulk.Goblet.cells.de_3MvsNeg,n=50)

#1Y vs HIV-
bulk.Goblet.cells.de_1YvsNeg <- FindMarkers(object = pseudo_exp,
                                                 ident.1 = "Goblet cells_HIV+ ART>1 Year",
                                                 ident.2 = "Goblet cells_HIV-",
                                                 test.use = "DESeq2")
head(bulk.Goblet.cells.de_1YvsNeg,n=50)


#3M vs 1Y
bulk.Goblet.cells.de_3Mvs1Y <- FindMarkers(object = pseudo_exp,
                                                ident.1 = "Goblet cells_HIV+ ART<3 Months",
                                                ident.2 = "Goblet cells_HIV+ ART>1 Year",
                                                test.use = "DESeq2")
head(bulk.Goblet.cells.de_3Mvs1Y,n=50)


#Plot a Volcano plot 3M vs HIV- ----

bulk.Goblet.cells.de_3MvsNeg.df <- as.data.frame(bulk.Goblet.cells.de_3MvsNeg)
bulk.Goblet.cells.de_3MvsNeg.df$genes <-rownames(bulk.Goblet.cells.de_3MvsNeg)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Goblet.cells.de_3MvsNeg.df$DiffExpressed <- "NO"

bulk.Goblet.cells.de_3MvsNeg.df$DiffExpressed[bulk.Goblet.cells.de_3MvsNeg.df$avg_log2FC >=1.5 &
                                                bulk.Goblet.cells.de_3MvsNeg.df$p_val<0.05] <- "UP"

bulk.Goblet.cells.de_3MvsNeg.df$DiffExpressed[bulk.Goblet.cells.de_3MvsNeg.df$avg_log2FC <= -1.5 &
                                                bulk.Goblet.cells.de_3MvsNeg.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Goblet.cells.de_3MvsNeg.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Goblet.cells.de_3MvsNeg.df$diffLabel[bulk.Goblet.cells.de_3MvsNeg.df$DiffExpressed !="NO"] <-
  bulk.Goblet.cells.de_3MvsNeg.df$genes[bulk.Goblet.cells.de_3MvsNeg.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Gobletcells3MvsNeg <- bulk.Goblet.cells.de_3MvsNeg.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  geom_vline(xintercept = c(-1.5,1.5), col="red")+
  geom_hline(yintercept = -log10(0.2), col="blue")+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Goblet cells",
       subtitle = "HIV+ ART<3 Months vs HIV-")+
  scale_y_continuous(limits = c(0,8))+
  scale_x_continuous(limits = c(-15,15))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Gobletcells3MvsNeg <- as.ggplot(Gobletcells3MvsNeg)
Gobletcells3MvsNeg

# Diff genes bargraph 3MvsNeg----
bulk.Goblet.cells.de_3MvsNeg.df_top20 <- bulk.Goblet.cells.de_3MvsNeg.df %>%
  filter(abs(avg_log2FC)>1.5,p_val<0.05,genes!="NA") 
bulk.Goblet.cells.de_3MvsNeg.df_top20 <- bulk.Goblet.cells.de_3MvsNeg.df_top20[order(abs(bulk.Goblet.cells.de_3MvsNeg.df_top20$avg_log2FC)),]
#bulk.Goblet.cells.de_3MvsNeg.df_top20 <- head(bulk.Goblet.cells.de_3MvsNeg.df_top20,50)


bulk.Goblet.cells.de_3MvsNeg.df_top20$diffLabel <- reorder(bulk.Goblet.cells.de_3MvsNeg.df_top20$diffLabel,
                                                           bulk.Goblet.cells.de_3MvsNeg.df_top20$avg_log2FC)



Goblet_3MvsNegcol <- bulk.Goblet.cells.de_3MvsNeg.df_top20 %>%
  ggplot(aes(diffLabel,avg_log2FC, fill=avg_log2FC
  ))+
  geom_col()+
  coord_flip()+
  theme_pubr()+
  #scale_y_continuous(limits = c(-8,8))+
  scale_fill_gradient2(low = "red",mid = "blue",high = "red")+
  labs(x="",y="Average log2FoldChange",
       title = "Goblet cells",
       subtitle = "HIV+ ART<3 Months vs HIV-")+
  guides(fill=guide_legend(title="Average log2FC"))+
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 20),
        plot.subtitle = element_text(hjust = 0.5,face = 'bold',size = 15,color='grey'))
Goblet_3MvsNegcol


#Plot a Volcano plot 1Y vs HIV- ----

bulk.Goblet.cells.de_1YvsNeg.df <- as.data.frame(bulk.Goblet.cells.de_1YvsNeg)
bulk.Goblet.cells.de_1YvsNeg.df$genes <-rownames(bulk.Goblet.cells.de_1YvsNeg)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Goblet.cells.de_1YvsNeg.df$DiffExpressed <- "NO"

bulk.Goblet.cells.de_1YvsNeg.df$DiffExpressed[bulk.Goblet.cells.de_1YvsNeg.df$avg_log2FC >=1.5 &
                                                bulk.Goblet.cells.de_1YvsNeg.df$p_val<0.05] <- "UP"

bulk.Goblet.cells.de_1YvsNeg.df$DiffExpressed[bulk.Goblet.cells.de_1YvsNeg.df$avg_log2FC <= -1.5 &
                                                bulk.Goblet.cells.de_1YvsNeg.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Goblet.cells.de_1YvsNeg.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Goblet.cells.de_1YvsNeg.df$diffLabel[bulk.Goblet.cells.de_1YvsNeg.df$DiffExpressed !="NO"] <-
  bulk.Goblet.cells.de_1YvsNeg.df$genes[bulk.Goblet.cells.de_1YvsNeg.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Gobletcells1YvsNeg <- bulk.Goblet.cells.de_1YvsNeg.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Goblet cells",
       subtitle = "HIV+ ART>1 Year vs HIV-")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Gobletcells1YvsNeg <- as.ggplot(Gobletcells1YvsNeg)
Gobletcells1YvsNeg

# Diff genes bargraph 1YvsNeg----
bulk.Goblet.cells.de_1YvsNeg.df_top20 <- bulk.Goblet.cells.de_1YvsNeg.df %>%
  filter(abs(avg_log2FC)>1.5,p_val<0.05,genes!="NA") 
bulk.Goblet.cells.de_1YvsNeg.df_top20 <- bulk.Goblet.cells.de_1YvsNeg.df_top20[order(abs(bulk.Goblet.cells.de_1YvsNeg.df_top20$avg_log2FC)),]
#bulk.Goblet.cells.de_1YvsNeg.df_top20 <- head(bulk.Goblet.cells.de_1YvsNeg.df_top20,50)


bulk.Goblet.cells.de_1YvsNeg.df_top20$diffLabel <- reorder(bulk.Goblet.cells.de_1YvsNeg.df_top20$diffLabel,
                                                           bulk.Goblet.cells.de_1YvsNeg.df_top20$avg_log2FC)



Goblet_1YvsNegcol <- bulk.Goblet.cells.de_1YvsNeg.df_top20 %>%
  ggplot(aes(diffLabel,avg_log2FC, fill=avg_log2FC
  ))+
  geom_col()+
  coord_flip()+
  theme_pubr()+
  #scale_y_continuous(limits = c(-8,8))+
  scale_fill_gradient2(low = "red",mid = "blue",high = "red")+
  labs(x="",y="Average log2FoldChange",
       title = "Goblet cells",
       subtitle = "HIV+ ART>1 Year vs HIV-")+
  guides(fill=guide_legend(title="Average log2FC"))+
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 20),
        plot.subtitle = element_text(hjust = 0.5,face = 'bold',size = 15,color='grey'))
Goblet_1YvsNegcol


#Plot a Volcano plot 3M vs 1Y ----

bulk.Goblet.cells.de_3Mvs1Y.df <- as.data.frame(bulk.Goblet.cells.de_3Mvs1Y)
bulk.Goblet.cells.de_3Mvs1Y.df$genes <-rownames(bulk.Goblet.cells.de_3Mvs1Y)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Goblet.cells.de_3Mvs1Y.df$DiffExpressed <- "NO"

bulk.Goblet.cells.de_3Mvs1Y.df$DiffExpressed[bulk.Goblet.cells.de_3Mvs1Y.df$avg_log2FC >=1.5 &
                                               bulk.Goblet.cells.de_3Mvs1Y.df$p_val<0.05] <- "UP"

bulk.Goblet.cells.de_3Mvs1Y.df$DiffExpressed[bulk.Goblet.cells.de_3Mvs1Y.df$avg_log2FC <= -1.5 &
                                               bulk.Goblet.cells.de_3Mvs1Y.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Goblet.cells.de_3Mvs1Y.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Goblet.cells.de_3Mvs1Y.df$diffLabel[bulk.Goblet.cells.de_3Mvs1Y.df$DiffExpressed !="NO"] <-
  bulk.Goblet.cells.de_3Mvs1Y.df$genes[bulk.Goblet.cells.de_3Mvs1Y.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Gobletcells3Mvs1Y <- bulk.Goblet.cells.de_3Mvs1Y.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Goblet cells",
       subtitle = "HIV+ ART<3 months vs HIV+ ART>1 Year")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Gobletcells3Mvs1Y <- as.ggplot(Gobletcells3Mvs1Y)
Gobletcells3Mvs1Y
Gobletcells3MvsNeg+Gobletcells1YvsNeg+Gobletcells3Mvs1Y


# Diff genes bargraph 3Mvs1Y----
bulk.Goblet.cells.de_3Mvs1Y.df_top20 <- bulk.Goblet.cells.de_3Mvs1Y.df %>%
  filter(abs(avg_log2FC)>1.5,p_val<0.05,genes!="NA") 
bulk.Goblet.cells.de_3Mvs1Y.df_top20 <- bulk.Goblet.cells.de_3Mvs1Y.df_top20[order(abs(bulk.Goblet.cells.de_3Mvs1Y.df_top20$avg_log2FC)),]
#bulk.Goblet.cells.de_3Mvs1Y.df_top20 <- head(bulk.Goblet.cells.de_3Mvs1Y.df_top20,50)


bulk.Goblet.cells.de_3Mvs1Y.df_top20$diffLabel <- reorder(bulk.Goblet.cells.de_3Mvs1Y.df_top20$diffLabel,
                                                          bulk.Goblet.cells.de_3Mvs1Y.df_top20$avg_log2FC)



Goblet_3Mvs1Ycol <- bulk.Goblet.cells.de_3Mvs1Y.df_top20 %>%
  ggplot(aes(diffLabel,avg_log2FC, fill=avg_log2FC
  ))+
  geom_col()+
  coord_flip()+
  theme_pubr()+
  #scale_y_continuous(limits = c(-8,8))+
  scale_fill_gradient2(low = "red",mid = "blue",high = "red")+
  labs(x="",y="Average log2FoldChange",
       title = "Goblet cells",
       subtitle = "HIV+ ART<3 Months vs HIV+ ART>1 Year")+
  guides(fill=guide_legend(title="Average log2FC"))+
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 20),
        plot.subtitle = element_text(hjust = 0.5,face = 'bold',size = 15,color='grey'))
Goblet_3Mvs1Ycol


# save the saved converted dimplot 
ggsave("results/GobletVolcanoPlots.pdf",
       plot = ((Gobletcells3MvsNeg|Gobletcells1YvsNeg|Gobletcells3Mvs1Y|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 8, unit="in", dpi = 700,limitsize = F)

# save the saved converted dimplot ----
ggsave("results/GobletbarPlots.pdf",
       plot = ((Goblet_3MvsNegcol|Goblet_1YvsNegcol|Goblet_3Mvs1Ycol|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 10, unit="in", dpi = 700,limitsize = F)






#Secretory cell DEG----
#3M vs HIV-
bulk.Secretory.cells.de_3MvsNeg <- FindMarkers(object = pseudo_exp,
                                            ident.1 = "Secretory cells_HIV+ ART<3 Months",
                                            ident.2 = "Secretory cells_HIV-",
                                            test.use = "DESeq2")
bulk.Secretory.cells.de_3MvsNeg$genes <- rownames(bulk.Secretory.cells.de_3MvsNeg)

head(bulk.Secretory.cells.de_3MvsNeg,n=50)

#1Y vs HIV-
bulk.Secretory.cells.de_1YvsNeg <- FindMarkers(object = pseudo_exp,
                                            ident.1 = "Secretory cells_HIV+ ART>1 Year",
                                            ident.2 = "Secretory cells_HIV-",
                                            test.use = "DESeq2")
bulk.Secretory.cells.de_1YvsNeg$genes <- rownames(bulk.Secretory.cells.de_1YvsNeg)

head(bulk.Secretory.cells.de_1YvsNeg,n=50)


#3M vs 1Y
bulk.Secretory.cells.de_3Mvs1Y <- FindMarkers(object = pseudo_exp,
                                           ident.1 = "Secretory cells_HIV+ ART<3 Months",
                                           ident.2 = "Secretory cells_HIV+ ART>1 Year",
                                           test.use = "DESeq2")
bulk.Secretory.cells.de_3Mvs1Y$genes <- rownames(bulk.Secretory.cells.de_3Mvs1Y)

head(bulk.Secretory.cells.de_3Mvs1Y,n=50)


#Plot a Volcano plot 3M vs HIV- ----

bulk.Secretory.cells.de_3MvsNeg.df <- as.data.frame(bulk.Secretory.cells.de_3MvsNeg)
bulk.Secretory.cells.de_3MvsNeg.df$genes <-rownames(bulk.Secretory.cells.de_3MvsNeg)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Secretory.cells.de_3MvsNeg.df$DiffExpressed <- "NO"

bulk.Secretory.cells.de_3MvsNeg.df$DiffExpressed[bulk.Secretory.cells.de_3MvsNeg.df$avg_log2FC >=1.5 &
                                                   bulk.Secretory.cells.de_3MvsNeg.df$p_val<0.05] <- "UP"

bulk.Secretory.cells.de_3MvsNeg.df$DiffExpressed[bulk.Secretory.cells.de_3MvsNeg.df$avg_log2FC <= -1.5 &
                                                   bulk.Secretory.cells.de_3MvsNeg.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Secretory.cells.de_3MvsNeg.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Secretory.cells.de_3MvsNeg.df$diffLabel[bulk.Secretory.cells.de_3MvsNeg.df$DiffExpressed !="NO"] <-
  bulk.Secretory.cells.de_3MvsNeg.df$genes[bulk.Secretory.cells.de_3MvsNeg.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Secretorycells3MvsNeg <- bulk.Secretory.cells.de_3MvsNeg.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.2), col="blue")+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Secretory cells",
       subtitle = "HIV+ ART<3 Months vs HIV-")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Secretorycells3MvsNeg <- as.ggplot(Secretorycells3MvsNeg)
Secretorycells3MvsNeg

# Diff genes bargraph 1YvsNeg----
bulk.Secretory.cells.de_3MvsNeg.df_top20 <- bulk.Secretory.cells.de_3MvsNeg.df %>%
  filter(abs(avg_log2FC)>1.5,p_val<0.05) 
bulk.Secretory.cells.de_3MvsNeg.df_top20 <- bulk.Secretory.cells.de_3MvsNeg.df_top20[order(bulk.Secretory.cells.de_3MvsNeg.df_top20$avg_log2FC, decreasing = T),]

bulk.Secretory.cells.de_3MvsNeg.df_top20$diffLabel <- reorder(bulk.Secretory.cells.de_3MvsNeg.df_top20$diffLabel,
                                                              bulk.Secretory.cells.de_3MvsNeg.df_top20$avg_log2FC)



Secretory_cells_3MvsNegcol <- bulk.Secretory.cells.de_3MvsNeg.df_top20 %>%
  ggplot(aes(diffLabel,avg_log2FC, fill=avg_log2FC
  ))+
  geom_col()+
  coord_flip()+
  theme_pubr()+
  #scale_y_continuous(limits = c(-8,8))+
  scale_fill_gradient2(low = "red",mid = "blue",high = "red")+
  labs(x="",y="Average log2FoldChange",
       title = "Secretory cells",
       subtitle = "HIV+ ART<3 months vs HIV-")+
  guides(fill=guide_legend(title="Average log2FC"))+
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 20),
        plot.subtitle = element_text(hjust = 0.5,face = 'bold',size = 15,color='grey'))
Secretory_cells_3MvsNegcol



#Plot a Volcano plot 1Y vs HIV- ----

bulk.Secretory.cells.de_1YvsNeg.df <- as.data.frame(bulk.Secretory.cells.de_1YvsNeg)
bulk.Secretory.cells.de_1YvsNeg.df$genes <-rownames(bulk.Secretory.cells.de_1YvsNeg)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Secretory.cells.de_1YvsNeg.df$DiffExpressed <- "NO"

bulk.Secretory.cells.de_1YvsNeg.df$DiffExpressed[bulk.Secretory.cells.de_1YvsNeg.df$avg_log2FC >=1.5 &
                                                   bulk.Secretory.cells.de_1YvsNeg.df$p_val<0.05] <- "UP"

bulk.Secretory.cells.de_1YvsNeg.df$DiffExpressed[bulk.Secretory.cells.de_1YvsNeg.df$avg_log2FC <= -1.5 &
                                                   bulk.Secretory.cells.de_1YvsNeg.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Secretory.cells.de_1YvsNeg.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Secretory.cells.de_1YvsNeg.df$diffLabel[bulk.Secretory.cells.de_1YvsNeg.df$DiffExpressed !="NO"] <-
  bulk.Secretory.cells.de_1YvsNeg.df$genes[bulk.Secretory.cells.de_1YvsNeg.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Secretorycells1YvsNeg <- bulk.Secretory.cells.de_1YvsNeg.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Secretory cells",
       subtitle = "HIV+ ART>1 Year vs HIV-")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Secretorycells1YvsNeg <- as.ggplot(Secretorycells1YvsNeg)
Secretorycells1YvsNeg

# Diff genes bargraph 1YvsNeg----
bulk.Secretory.cells.de_1YvsNeg.df_top20 <- bulk.Secretory.cells.de_1YvsNeg.df %>%
  filter(abs(avg_log2FC)>1.5,p_val<0.05) 
bulk.Secretory.cells.de_1YvsNeg.df_top20 <- bulk.Secretory.cells.de_1YvsNeg.df_top20[order(bulk.Secretory.cells.de_1YvsNeg.df_top20$avg_log2FC, decreasing = T),]

bulk.Secretory.cells.de_1YvsNeg.df_top20$diffLabel <- reorder(bulk.Secretory.cells.de_1YvsNeg.df_top20$diffLabel,
                                                              bulk.Secretory.cells.de_1YvsNeg.df_top20$avg_log2FC)



Secretory_cells_1YvsNegcol <- bulk.Secretory.cells.de_1YvsNeg.df_top20 %>%
  ggplot(aes(diffLabel,avg_log2FC, fill=avg_log2FC
  ))+
  geom_col()+
  coord_flip()+
  theme_pubr()+
  #scale_y_continuous(limits = c(-8,8))+
  scale_fill_gradient2(low = "red",mid = "blue",high = "red")+
  labs(x="",y="Average log2FoldChange",
       title = "Secretory cells",
       subtitle = "HIV+ ART>1 Year vs HIV-")+
  guides(fill=guide_legend(title="Average log2FC"))+
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 20),
        plot.subtitle = element_text(hjust = 0.5,face = 'bold',size = 15,color='grey'))
Secretory_cells_1YvsNegcol

#Plot a Volcano plot 3M vs 1Y ----

bulk.Secretory.cells.de_3Mvs1Y.df <- as.data.frame(bulk.Secretory.cells.de_3Mvs1Y)
bulk.Secretory.cells.de_3Mvs1Y.df$genes <-rownames(bulk.Secretory.cells.de_3Mvs1Y)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Secretory.cells.de_3Mvs1Y.df$DiffExpressed <- "NO"

bulk.Secretory.cells.de_3Mvs1Y.df$DiffExpressed[bulk.Secretory.cells.de_3Mvs1Y.df$avg_log2FC >=1.5 &
                                                  bulk.Secretory.cells.de_3Mvs1Y.df$p_val<0.05] <- "UP"

bulk.Secretory.cells.de_3Mvs1Y.df$DiffExpressed[bulk.Secretory.cells.de_3Mvs1Y.df$avg_log2FC <= -1.5 &
                                                  bulk.Secretory.cells.de_3Mvs1Y.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Secretory.cells.de_3Mvs1Y.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Secretory.cells.de_3Mvs1Y.df$diffLabel[bulk.Secretory.cells.de_3Mvs1Y.df$DiffExpressed !="NO"] <-
  bulk.Secretory.cells.de_3Mvs1Y.df$genes[bulk.Secretory.cells.de_3Mvs1Y.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Secretorycells3Mvs1Y <- bulk.Secretory.cells.de_3Mvs1Y.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Secretory cells",
       subtitle = "HIV+ ART<3 months vs HIV+ ART>1 Year")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Secretorycells3Mvs1Y <- as.ggplot(Secretorycells3Mvs1Y)
Secretorycells3Mvs1Y
Secretorycells3MvsNeg+Secretorycells1YvsNeg+Secretorycells3Mvs1Y


# Diff genes bargraph 3Mvs1Y----
bulk.Secretory.cells.de_3Mvs1Y.df_top20 <- bulk.Secretory.cells.de_3Mvs1Y.df %>%
  filter(abs(avg_log2FC)>1.5,p_val<0.05) 
bulk.Secretory.cells.de_3Mvs1Y.df_top20 <- bulk.Secretory.cells.de_3Mvs1Y.df_top20[order(bulk.Secretory.cells.de_3Mvs1Y.df_top20$avg_log2FC, decreasing = T),]

bulk.Secretory.cells.de_3Mvs1Y.df_top20$diffLabel <- reorder(bulk.Secretory.cells.de_3Mvs1Y.df_top20$diffLabel,
                                                             bulk.Secretory.cells.de_3Mvs1Y.df_top20$avg_log2FC)



Secretory_cells_3Mvs1Ycol <- bulk.Secretory.cells.de_3Mvs1Y.df_top20 %>%
  ggplot(aes(diffLabel,avg_log2FC, fill=avg_log2FC
  ))+
  geom_col()+
  coord_flip()+
  theme_pubr()+
  #scale_y_continuous(limits = c(-8,8))+
  scale_fill_gradient2(low = "red",mid = "blue",high = "red")+
  labs(x="",y="Average log2FoldChange",
       title = "Secretory cells",
       subtitle = "HIV+ ART<3 months vs HIV+ ART>1 Year")+
  guides(fill=guide_legend(title="Average log2FC"))+
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 20),
        plot.subtitle = element_text(hjust = 0.5,face = 'bold',size = 15,color='grey'))
Secretory_cells_3Mvs1Ycol


# save the saved converted dimplot 
ggsave("results/SecretoryVolcanoPlots.pdf",
       plot = ((Secretorycells3MvsNeg|Secretorycells1YvsNeg|Secretorycells3Mvs1Y|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 8, unit="in", dpi = 700,limitsize = F)

# save the saved converted dimplot ----
ggsave("results/SecretorybarPlots.pdf",
       plot = ((Secretory_cells_3MvsNegcol|Secretory_cells_1YvsNegcol|Secretory_cells_3Mvs1Ycol|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 10, unit="in", dpi = 700,limitsize = F)








#Basal cell DEG----
#3M vs HIV-
bulk.Basal.cells.de_3MvsNeg <- FindMarkers(object = pseudo_exp,
                                               ident.1 = "Basal cells_HIV+ ART<3 Months",
                                               ident.2 = "Basal cells_HIV-",
                                               test.use = "DESeq2")
head(bulk.Basal.cells.de_3MvsNeg,n=50)

#1Y vs HIV-
bulk.Basal.cells.de_1YvsNeg <- FindMarkers(object = pseudo_exp,
                                               ident.1 = "Basal cells_HIV+ ART>1 Year",
                                               ident.2 = "Basal cells_HIV-",
                                               test.use = "DESeq2")
head(bulk.Basal.cells.de_1YvsNeg,n=50)


#3M vs 1Y
bulk.Basal.cells.de_3Mvs1Y <- FindMarkers(object = pseudo_exp,
                                              ident.1 = "Basal cells_HIV+ ART<3 Months",
                                              ident.2 = "Basal cells_HIV+ ART>1 Year",
                                              test.use = "DESeq2")
head(bulk.Basal.cells.de_3Mvs1Y,n=50)


#Plot a Volcano plot 3M vs HIV- ----

bulk.Basal.cells.de_3MvsNeg.df <- as.data.frame(bulk.Basal.cells.de_3MvsNeg)
bulk.Basal.cells.de_3MvsNeg.df$genes <-rownames(bulk.Basal.cells.de_3MvsNeg)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Basal.cells.de_3MvsNeg.df$DiffExpressed <- "NO"

bulk.Basal.cells.de_3MvsNeg.df$DiffExpressed[bulk.Basal.cells.de_3MvsNeg.df$avg_log2FC >=1.5 &
                                               bulk.Basal.cells.de_3MvsNeg.df$p_val<0.05] <- "UP"

bulk.Basal.cells.de_3MvsNeg.df$DiffExpressed[bulk.Basal.cells.de_3MvsNeg.df$avg_log2FC <= -1.5 &
                                               bulk.Basal.cells.de_3MvsNeg.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Basal.cells.de_3MvsNeg.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Basal.cells.de_3MvsNeg.df$diffLabel[bulk.Basal.cells.de_3MvsNeg.df$DiffExpressed !="NO"] <-
  bulk.Basal.cells.de_3MvsNeg.df$genes[bulk.Basal.cells.de_3MvsNeg.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Basalcells3MvsNeg <- bulk.Basal.cells.de_3MvsNeg.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.2), col="blue")+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Basal cells",
       subtitle = "HIV+ ART<3 Months vs HIV-")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Basalcells3MvsNeg <- as.ggplot(Basalcells3MvsNeg)
Basalcells3MvsNeg
# Diff genes bargraph 3MvsNeg----
bulk.Basal.cells.de_3MvsNeg.df_top20 <- bulk.Basal.cells.de_3MvsNeg.df %>%
  filter(abs(avg_log2FC)>1.5,p_val<0.05,genes!="NA") 
bulk.Basal.cells.de_3MvsNeg.df_top20 <- bulk.Basal.cells.de_3MvsNeg.df_top20[order(abs(bulk.Basal.cells.de_3MvsNeg.df_top20$avg_log2FC)),]
bulk.Basal.cells.de_3MvsNeg.df_top20 <- head(bulk.Basal.cells.de_3MvsNeg.df_top20,50)


bulk.Basal.cells.de_3MvsNeg.df_top20$diffLabel <- reorder(bulk.Basal.cells.de_3MvsNeg.df_top20$diffLabel,
                                                          bulk.Basal.cells.de_3MvsNeg.df_top20$avg_log2FC)



Basal_cells_3MvsNegcol <- bulk.Basal.cells.de_3MvsNeg.df_top20 %>%
  ggplot(aes(diffLabel,avg_log2FC, fill=avg_log2FC
  ))+
  geom_col()+
  coord_flip()+
  theme_pubr()+
  #scale_y_continuous(limits = c(-8,8))+
  scale_fill_gradient2(low = "red",mid = "blue",high = "red")+
  labs(x="",y="Average log2FoldChange",
       title = "Basal cells",
       subtitle = "HIV+ ART<3 Months vs HIV-")+
  guides(fill=guide_legend(title="Average log2FC"))+
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 20),
        plot.subtitle = element_text(hjust = 0.5,face = 'bold',size = 15,color='grey'))
Basal_cells_3MvsNegcol


#Plot a Volcano plot 1Y vs HIV- ----

bulk.Basal.cells.de_1YvsNeg.df <- as.data.frame(bulk.Basal.cells.de_1YvsNeg)
bulk.Basal.cells.de_1YvsNeg.df$genes <-rownames(bulk.Basal.cells.de_1YvsNeg)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Basal.cells.de_1YvsNeg.df$DiffExpressed <- "NO"

bulk.Basal.cells.de_1YvsNeg.df$DiffExpressed[bulk.Basal.cells.de_1YvsNeg.df$avg_log2FC >=1.5 &
                                               bulk.Basal.cells.de_1YvsNeg.df$p_val<0.05] <- "UP"

bulk.Basal.cells.de_1YvsNeg.df$DiffExpressed[bulk.Basal.cells.de_1YvsNeg.df$avg_log2FC <= -1.5 &
                                               bulk.Basal.cells.de_1YvsNeg.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Basal.cells.de_1YvsNeg.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Basal.cells.de_1YvsNeg.df$diffLabel[bulk.Basal.cells.de_1YvsNeg.df$DiffExpressed !="NO"] <-
  bulk.Basal.cells.de_1YvsNeg.df$genes[bulk.Basal.cells.de_1YvsNeg.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Basalcells1YvsNeg <- bulk.Basal.cells.de_1YvsNeg.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Basal cells",
       subtitle = "HIV+ ART>1 Year vs HIV-")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Basalcells1YvsNeg <- as.ggplot(Basalcells1YvsNeg)
Basalcells1YvsNeg

# Diff genes bargraph 1YvsNeg----
bulk.Basal.cells.de_1YvsNeg.df_top20 <- bulk.Basal.cells.de_1YvsNeg.df %>%
  filter(abs(avg_log2FC)>1.5,p_val<0.05) 
bulk.Basal.cells.de_1YvsNeg.df_top20 <- bulk.Basal.cells.de_1YvsNeg.df_top20[order(bulk.Basal.cells.de_1YvsNeg.df_top20$avg_log2FC, decreasing = T),]

bulk.Basal.cells.de_1YvsNeg.df_top20$diffLabel <- reorder(bulk.Basal.cells.de_1YvsNeg.df_top20$diffLabel,
                                                          bulk.Basal.cells.de_1YvsNeg.df_top20$avg_log2FC)



Basal_cells_1YvsNegcol <- bulk.Basal.cells.de_1YvsNeg.df_top20 %>%
  ggplot(aes(diffLabel,avg_log2FC, fill=avg_log2FC
  ))+
  geom_col()+
  coord_flip()+
  theme_pubr()+
  #scale_y_continuous(limits = c(-8,8))+
  scale_fill_gradient2(low = "red",mid = "blue",high = "red")+
  labs(x="",y="Average log2FoldChange",
       title = "Basal cells",
       subtitle = "HIV+ ART>1 Year vs HIV-")+
  guides(fill=guide_legend(title="Average log2FC"))+
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 20),
        plot.subtitle = element_text(hjust = 0.5,face = 'bold',size = 15,color='grey'))
Basal_cells_1YvsNegcol


#Plot a Volcano plot 3M vs 1Y ----

bulk.Basal.cells.de_3Mvs1Y.df <- as.data.frame(bulk.Basal.cells.de_3Mvs1Y)
bulk.Basal.cells.de_3Mvs1Y.df$genes <-rownames(bulk.Basal.cells.de_3Mvs1Y)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Basal.cells.de_3Mvs1Y.df$DiffExpressed <- "NO"

bulk.Basal.cells.de_3Mvs1Y.df$DiffExpressed[bulk.Basal.cells.de_3Mvs1Y.df$avg_log2FC >=1.5 &
                                              bulk.Basal.cells.de_3Mvs1Y.df$p_val<0.05] <- "UP"

bulk.Basal.cells.de_3Mvs1Y.df$DiffExpressed[bulk.Basal.cells.de_3Mvs1Y.df$avg_log2FC <= -1.5 &
                                              bulk.Basal.cells.de_3Mvs1Y.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Basal.cells.de_3Mvs1Y.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Basal.cells.de_3Mvs1Y.df$diffLabel[bulk.Basal.cells.de_3Mvs1Y.df$DiffExpressed !="NO"] <-
  bulk.Basal.cells.de_3Mvs1Y.df$genes[bulk.Basal.cells.de_3Mvs1Y.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Basalcells3Mvs1Y <- bulk.Basal.cells.de_3Mvs1Y.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Basal cells",
       subtitle = "HIV+ ART<3 months vs HIV+ ART>1 Year")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Basalcells3Mvs1Y <- as.ggplot(Basalcells3Mvs1Y)
Basalcells3Mvs1Y
Basalcells3MvsNeg+Basalcells1YvsNeg+Basalcells3Mvs1Y


# Diff genes bargraph 3Mvs1Y----
bulk.Basal.cells.de_3Mvs1Y.df_top20 <- bulk.Basal.cells.de_3Mvs1Y.df %>%
  filter(abs(avg_log2FC)>1.5,p_val<0.05) 
bulk.Basal.cells.de_3Mvs1Y.df_top20 <- bulk.Basal.cells.de_3Mvs1Y.df_top20[order(bulk.Basal.cells.de_3Mvs1Y.df_top20$avg_log2FC, decreasing = T),]

bulk.Basal.cells.de_3Mvs1Y.df_top20$diffLabel <- reorder(bulk.Basal.cells.de_3Mvs1Y.df_top20$diffLabel,
                                                         bulk.Basal.cells.de_3Mvs1Y.df_top20$avg_log2FC)



Basal_cells_3Mvs1Ycol <- bulk.Basal.cells.de_3Mvs1Y.df_top20 %>%
  ggplot(aes(diffLabel,avg_log2FC, fill=avg_log2FC
  ))+
  geom_col()+
  coord_flip()+
  theme_pubr()+
  #scale_y_continuous(limits = c(-8,8))+
  scale_fill_gradient2(low = "red",mid = "blue",high = "red")+
  labs(x="",y="Average log2FoldChange",
       title = "Basal cells",
       subtitle = "HIV+ ART<3 Months vs HIV+ ART>1 Year")+
  guides(fill=guide_legend(title="Average log2FC"))+
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 20),
        plot.subtitle = element_text(hjust = 0.5,face = 'bold',size = 15,color='grey'))
Basal_cells_3Mvs1Ycol



# save the saved converted dimplot 
ggsave("results/BasalVolcanoPlots.pdf",
       plot = ((Basalcells3MvsNeg|Basalcells1YvsNeg|Basalcells3Mvs1Y|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 8, unit="in", dpi = 700,limitsize = F)


#Dendritic cell DEG----
#3M vs HIV-
bulk.Dendritic.cells.de_3MvsNeg <- FindMarkers(object = pseudo_exp,
                                           ident.1 = "Dendritic cells_HIV+ ART<3 Months",
                                           ident.2 = "Dendritic cells_HIV-",
                                           test.use = "DESeq2")
head(bulk.Dendritic.cells.de_3MvsNeg,n=50)

#1Y vs HIV-
bulk.Dendritic.cells.de_1YvsNeg <- FindMarkers(object = pseudo_exp,
                                           ident.1 = "Dendritic cells_HIV+ ART>1 Year",
                                           ident.2 = "Dendritic cells_HIV-",
                                           test.use = "DESeq2")
head(bulk.Dendritic.cells.de_1YvsNeg,n=50)


#3M vs 1Y
bulk.Dendritic.cells.de_3Mvs1Y <- FindMarkers(object = pseudo_exp,
                                          ident.1 = "Dendritic cells_HIV+ ART<3 Months",
                                          ident.2 = "Dendritic cells_HIV+ ART>1 Year",
                                          test.use = "DESeq2")
head(bulk.Dendritic.cells.de_3Mvs1Y,n=50)


#Plot a Volcano plot 3M vs HIV- ----

bulk.Dendritic.cells.de_3MvsNeg.df <- as.data.frame(bulk.Dendritic.cells.de_3MvsNeg)
bulk.Dendritic.cells.de_3MvsNeg.df$genes <-rownames(bulk.Dendritic.cells.de_3MvsNeg)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Dendritic.cells.de_3MvsNeg.df$DiffExpressed <- "NO"

bulk.Dendritic.cells.de_3MvsNeg.df$DiffExpressed[bulk.Dendritic.cells.de_3MvsNeg.df$avg_log2FC >=1.5 &
                                                   bulk.Dendritic.cells.de_3MvsNeg.df$p_val<0.05] <- "UP"

bulk.Dendritic.cells.de_3MvsNeg.df$DiffExpressed[bulk.Dendritic.cells.de_3MvsNeg.df$avg_log2FC <= -1.5 &
                                                   bulk.Dendritic.cells.de_3MvsNeg.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Dendritic.cells.de_3MvsNeg.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Dendritic.cells.de_3MvsNeg.df$diffLabel[bulk.Dendritic.cells.de_3MvsNeg.df$DiffExpressed !="NO"] <-
  bulk.Dendritic.cells.de_3MvsNeg.df$genes[bulk.Dendritic.cells.de_3MvsNeg.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Dendriticcells3MvsNeg <- bulk.Dendritic.cells.de_3MvsNeg.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Dendritic cells",
       subtitle = "HIV+ ART<3 Months vs HIV-")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Dendriticcells3MvsNeg <- as.ggplot(Dendriticcells3MvsNeg)
Dendriticcells3MvsNeg



#Plot a Volcano plot 1Y vs HIV- ----

bulk.Dendritic.cells.de_1YvsNeg.df <- as.data.frame(bulk.Dendritic.cells.de_1YvsNeg)
bulk.Dendritic.cells.de_1YvsNeg.df$genes <-rownames(bulk.Dendritic.cells.de_1YvsNeg)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Dendritic.cells.de_1YvsNeg.df$DiffExpressed <- "NO"

bulk.Dendritic.cells.de_1YvsNeg.df$DiffExpressed[bulk.Dendritic.cells.de_1YvsNeg.df$avg_log2FC >=1.5 &
                                                   bulk.Dendritic.cells.de_1YvsNeg.df$p_val<0.05] <- "UP"

bulk.Dendritic.cells.de_1YvsNeg.df$DiffExpressed[bulk.Dendritic.cells.de_1YvsNeg.df$avg_log2FC <= -1.5 &
                                                   bulk.Dendritic.cells.de_1YvsNeg.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Dendritic.cells.de_1YvsNeg.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Dendritic.cells.de_1YvsNeg.df$diffLabel[bulk.Dendritic.cells.de_1YvsNeg.df$DiffExpressed !="NO"] <-
  bulk.Dendritic.cells.de_1YvsNeg.df$genes[bulk.Dendritic.cells.de_1YvsNeg.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Dendriticcells1YvsNeg <- bulk.Dendritic.cells.de_1YvsNeg.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Dendritic cells",
       subtitle = "HIV+ ART>1 Year vs HIV-")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Dendriticcells1YvsNeg <- as.ggplot(Dendriticcells1YvsNeg)
Dendriticcells1YvsNeg

#Plot a Volcano plot 3M vs 1Y ----

bulk.Dendritic.cells.de_3Mvs1Y.df <- as.data.frame(bulk.Dendritic.cells.de_3Mvs1Y)
bulk.Dendritic.cells.de_3Mvs1Y.df$genes <-rownames(bulk.Dendritic.cells.de_3Mvs1Y)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Dendritic.cells.de_3Mvs1Y.df$DiffExpressed <- "NO"

bulk.Dendritic.cells.de_3Mvs1Y.df$DiffExpressed[bulk.Dendritic.cells.de_3Mvs1Y.df$avg_log2FC >=1.5 &
                                                  bulk.Dendritic.cells.de_3Mvs1Y.df$p_val<0.05] <- "UP"

bulk.Dendritic.cells.de_3Mvs1Y.df$DiffExpressed[bulk.Dendritic.cells.de_3Mvs1Y.df$avg_log2FC <= -1.5 &
                                                  bulk.Dendritic.cells.de_3Mvs1Y.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Dendritic.cells.de_3Mvs1Y.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Dendritic.cells.de_3Mvs1Y.df$diffLabel[bulk.Dendritic.cells.de_3Mvs1Y.df$DiffExpressed !="NO"] <-
  bulk.Dendritic.cells.de_3Mvs1Y.df$genes[bulk.Dendritic.cells.de_3Mvs1Y.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Dendriticcells3Mvs1Y <- bulk.Dendritic.cells.de_3Mvs1Y.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Dendritic cells",
       subtitle = "HIV+ ART<3 months vs HIV+ ART>1 Year")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Dendriticcells3Mvs1Y <- as.ggplot(Dendriticcells3Mvs1Y)
Dendriticcells3Mvs1Y
Dendriticcells3MvsNeg+Dendriticcells1YvsNeg+Dendriticcells3Mvs1Y

# save the saved converted dimplot 
ggsave("results/DendriticVolcanoPlots.pdf",
       plot = ((Dendriticcells3MvsNeg|Dendriticcells1YvsNeg|Dendriticcells3Mvs1Y|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 8, unit="in", dpi = 700,limitsize = F)























#Phagocytes cell DEG----
#3M vs HIV-
bulk.Phagocytes.cells.de_3MvsNeg <- FindMarkers(object = pseudo_exp,
                                               ident.1 = "Phagocytes_HIV+ ART<3 Months",
                                               ident.2 = "Phagocytes_HIV-",
                                               test.use = "DESeq2")
head(bulk.Phagocytes.cells.de_3MvsNeg,n=50)

#1Y vs HIV-
bulk.Phagocytes.cells.de_1YvsNeg <- FindMarkers(object = pseudo_exp,
                                               ident.1 = "Phagocytes_HIV+ ART>1 Year",
                                               ident.2 = "Phagocytes_HIV-",
                                               test.use = "DESeq2")
head(bulk.Phagocytes.cells.de_1YvsNeg,n=50)


#3M vs 1Y
bulk.Phagocytes.cells.de_3Mvs1Y <- FindMarkers(object = pseudo_exp,
                                              ident.1 = "Phagocytes_HIV+ ART<3 Months",
                                              ident.2 = "Phagocytes_HIV+ ART>1 Year",
                                              test.use = "DESeq2")
head(bulk.Phagocytes.cells.de_3Mvs1Y,n=50)


#Plot a Volcano plot 3M vs HIV- ----

bulk.Phagocytes.cells.de_3MvsNeg.df <- as.data.frame(bulk.Phagocytes.cells.de_3MvsNeg)
bulk.Phagocytes.cells.de_3MvsNeg.df$genes <-rownames(bulk.Phagocytes.cells.de_3MvsNeg)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Phagocytes.cells.de_3MvsNeg.df$DiffExpressed <- "NO"

bulk.Phagocytes.cells.de_3MvsNeg.df$DiffExpressed[bulk.Phagocytes.cells.de_3MvsNeg.df$avg_log2FC >=1.5 &
                                                    bulk.Phagocytes.cells.de_3MvsNeg.df$p_val<0.05] <- "UP"

bulk.Phagocytes.cells.de_3MvsNeg.df$DiffExpressed[bulk.Phagocytes.cells.de_3MvsNeg.df$avg_log2FC <= -1.5 &
                                                   bulk.Phagocytes.cells.de_3MvsNeg.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Phagocytes.cells.de_3MvsNeg.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Phagocytes.cells.de_3MvsNeg.df$diffLabel[bulk.Phagocytes.cells.de_3MvsNeg.df$DiffExpressed !="NO"] <-
  bulk.Phagocytes.cells.de_3MvsNeg.df$genes[bulk.Phagocytes.cells.de_3MvsNeg.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Phagocytescells3MvsNeg <- bulk.Phagocytes.cells.de_3MvsNeg.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  geom_vline(xintercept = c(-1.5,1.5), col="red")+
  geom_hline(yintercept = -log10(0.2), col="blue")+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Phagocytes cells",
       subtitle = "HIV+ ART<3 Months vs HIV-")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Phagocytescells3MvsNeg <- as.ggplot(Phagocytescells3MvsNeg)
Phagocytescells3MvsNeg



#Plot a Volcano plot 1Y vs HIV- ----

bulk.Phagocytes.cells.de_1YvsNeg.df <- as.data.frame(bulk.Phagocytes.cells.de_1YvsNeg)
bulk.Phagocytes.cells.de_1YvsNeg.df$genes <-rownames(bulk.Phagocytes.cells.de_1YvsNeg)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Phagocytes.cells.de_1YvsNeg.df$DiffExpressed <- "NO"

bulk.Phagocytes.cells.de_1YvsNeg.df$DiffExpressed[bulk.Phagocytes.cells.de_1YvsNeg.df$avg_log2FC >=1.5 &
                                                    bulk.Phagocytes.cells.de_1YvsNeg.df$p_val<0.05] <- "UP"

bulk.Phagocytes.cells.de_1YvsNeg.df$DiffExpressed[bulk.Phagocytes.cells.de_1YvsNeg.df$avg_log2FC <= -1.5 &
                                                   bulk.Phagocytes.cells.de_1YvsNeg.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Phagocytes.cells.de_1YvsNeg.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Phagocytes.cells.de_1YvsNeg.df$diffLabel[bulk.Phagocytes.cells.de_1YvsNeg.df$DiffExpressed !="NO"] <-
  bulk.Phagocytes.cells.de_1YvsNeg.df$genes[bulk.Phagocytes.cells.de_1YvsNeg.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Phagocytescells1YvsNeg <- bulk.Phagocytes.cells.de_1YvsNeg.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Phagocytes cells",
       subtitle = "HIV+ ART>1 Year vs HIV-")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Phagocytescells1YvsNeg <- as.ggplot(Phagocytescells1YvsNeg)
Phagocytescells1YvsNeg

#Plot a Volcano plot 3M vs 1Y ----

bulk.Phagocytes.cells.de_3Mvs1Y.df <- as.data.frame(bulk.Phagocytes.cells.de_3Mvs1Y)
bulk.Phagocytes.cells.de_3Mvs1Y.df$genes <-rownames(bulk.Phagocytes.cells.de_3Mvs1Y)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Phagocytes.cells.de_3Mvs1Y.df$DiffExpressed <- "NO"

bulk.Phagocytes.cells.de_3Mvs1Y.df$DiffExpressed[bulk.Phagocytes.cells.de_3Mvs1Y.df$avg_log2FC >=1.5 &
                                                  bulk.Phagocytes.cells.de_3Mvs1Y.df$p_val<0.05] <- "UP"

bulk.Phagocytes.cells.de_3Mvs1Y.df$DiffExpressed[bulk.Phagocytes.cells.de_3Mvs1Y.df$avg_log2FC <= -1.5 &
                                                  bulk.Phagocytes.cells.de_3Mvs1Y.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Phagocytes.cells.de_3Mvs1Y.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Phagocytes.cells.de_3Mvs1Y.df$diffLabel[bulk.Phagocytes.cells.de_3Mvs1Y.df$DiffExpressed !="NO"] <-
  bulk.Phagocytes.cells.de_3Mvs1Y.df$genes[bulk.Phagocytes.cells.de_3Mvs1Y.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Phagocytescells3Mvs1Y <- bulk.Phagocytes.cells.de_3Mvs1Y.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Phagocytes cells",
       subtitle = "HIV+ ART<3 months vs HIV+ ART>1 Year")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Phagocytescells3Mvs1Y <- as.ggplot(Phagocytescells3Mvs1Y)
Phagocytescells3Mvs1Y
Dendriticcells3MvsNeg+Dendriticcells1YvsNeg+Phagocytescells3Mvs1Y

# save the saved converted dimplot 
ggsave("results/PhagocytesVolcanoPlots.pdf",
       plot = ((Phagocytescells3MvsNeg|Phagocytescells1YvsNeg|Phagocytescells3Mvs1Y|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 8, unit="in", dpi = 700,limitsize = F)


#Squamous cell DEG----
#3M vs HIV-
bulk.Squamous.cells.de_3MvsNeg <- FindMarkers(object = pseudo_exp,
                                                ident.1 = "Squamous cells_HIV+ ART<3 Months",
                                                ident.2 = "Squamous cells_HIV-",
                                                test.use = "DESeq2")
head(bulk.Squamous.cells.de_3MvsNeg,n=50)

#1Y vs HIV-
bulk.Squamous.cells.de_1YvsNeg <- FindMarkers(object = pseudo_exp,
                                                ident.1 = "Squamous cells_HIV+ ART>1 Year",
                                                ident.2 = "Squamous cells_HIV-",
                                                test.use = "DESeq2")
head(bulk.Squamous.cells.de_1YvsNeg,n=50)


#3M vs 1Y
bulk.Squamous.cells.de_3Mvs1Y <- FindMarkers(object = pseudo_exp,
                                               ident.1 = "Squamous cells_HIV+ ART<3 Months",
                                               ident.2 = "Squamous cells_HIV+ ART>1 Year",
                                               test.use = "DESeq2")
head(bulk.Squamous.cells.de_3Mvs1Y,n=50)


#Plot a Volcano plot 3M vs HIV- ----

bulk.Squamous.cells.de_3MvsNeg.df <- as.data.frame(bulk.Squamous.cells.de_3MvsNeg)
bulk.Squamous.cells.de_3MvsNeg.df$genes <-rownames(bulk.Squamous.cells.de_3MvsNeg)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Squamous.cells.de_3MvsNeg.df$DiffExpressed <- "NO"

bulk.Squamous.cells.de_3MvsNeg.df$DiffExpressed[bulk.Squamous.cells.de_3MvsNeg.df$avg_log2FC >=1.5 &
                                                  bulk.Squamous.cells.de_3MvsNeg.df$p_val<0.05] <- "UP"

bulk.Squamous.cells.de_3MvsNeg.df$DiffExpressed[bulk.Squamous.cells.de_3MvsNeg.df$avg_log2FC <= -1.5 &
                                                  bulk.Squamous.cells.de_3MvsNeg.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Squamous.cells.de_3MvsNeg.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Squamous.cells.de_3MvsNeg.df$diffLabel[bulk.Squamous.cells.de_3MvsNeg.df$DiffExpressed !="NO"] <-
  bulk.Squamous.cells.de_3MvsNeg.df$genes[bulk.Squamous.cells.de_3MvsNeg.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Squamouscells3MvsNeg <- bulk.Squamous.cells.de_3MvsNeg.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  geom_vline(xintercept = c(-1.5,1.5), col="red")+
  geom_hline(yintercept = -log10(0.2), col="blue")+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Squamous cells",
       subtitle = "HIV+ ART<3 Months vs HIV-")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Squamouscells3MvsNeg <- as.ggplot(Squamouscells3MvsNeg)
Squamouscells3MvsNeg



#Plot a Volcano plot 1Y vs HIV- ----

bulk.Squamous.cells.de_1YvsNeg.df <- as.data.frame(bulk.Squamous.cells.de_1YvsNeg)
bulk.Squamous.cells.de_1YvsNeg.df$genes <-rownames(bulk.Squamous.cells.de_1YvsNeg)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Squamous.cells.de_1YvsNeg.df$DiffExpressed <- "NO"

bulk.Squamous.cells.de_1YvsNeg.df$DiffExpressed[bulk.Squamous.cells.de_1YvsNeg.df$avg_log2FC >=1.5 &
                                                  bulk.Squamous.cells.de_1YvsNeg.df$p_val<0.05] <- "UP"

bulk.Squamous.cells.de_1YvsNeg.df$DiffExpressed[bulk.Squamous.cells.de_1YvsNeg.df$avg_log2FC <= -1.5 &
                                                  bulk.Squamous.cells.de_1YvsNeg.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Squamous.cells.de_1YvsNeg.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Squamous.cells.de_1YvsNeg.df$diffLabel[bulk.Squamous.cells.de_1YvsNeg.df$DiffExpressed !="NO"] <-
  bulk.Squamous.cells.de_1YvsNeg.df$genes[bulk.Squamous.cells.de_1YvsNeg.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Squamouscells1YvsNeg <- bulk.Squamous.cells.de_1YvsNeg.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Squamous cells",
       subtitle = "HIV+ ART>1 Year vs HIV-")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Squamouscells1YvsNeg <- as.ggplot(Squamouscells1YvsNeg)
Squamouscells1YvsNeg

#Plot a Volcano plot 3M vs 1Y ----

bulk.Squamous.cells.de_3Mvs1Y.df <- as.data.frame(bulk.Squamous.cells.de_3Mvs1Y)
bulk.Squamous.cells.de_3Mvs1Y.df$genes <-rownames(bulk.Squamous.cells.de_3Mvs1Y)
#Create new columns based on UP or Down regulated genes based on P value and log2FoldChange 

bulk.Squamous.cells.de_3Mvs1Y.df$DiffExpressed <- "NO"

bulk.Squamous.cells.de_3Mvs1Y.df$DiffExpressed[bulk.Squamous.cells.de_3Mvs1Y.df$avg_log2FC >=1.5 &
                                                 bulk.Squamous.cells.de_3Mvs1Y.df$p_val<0.05] <- "UP"

bulk.Squamous.cells.de_3Mvs1Y.df$DiffExpressed[bulk.Squamous.cells.de_3Mvs1Y.df$avg_log2FC <= -1.5 &
                                                 bulk.Squamous.cells.de_3Mvs1Y.df$p_val<0.05] <- "DOWN"

#defining that gene labels should only be added  differentially expressed genes. if not, dont add label
#Creating an empty column
bulk.Squamous.cells.de_3Mvs1Y.df$diffLabel <- NA
# Adding new labels to differentially expressed genes in the new column. if not, dont add label
bulk.Squamous.cells.de_3Mvs1Y.df$diffLabel[bulk.Squamous.cells.de_3Mvs1Y.df$DiffExpressed !="NO"] <-
  bulk.Squamous.cells.de_3Mvs1Y.df$genes[bulk.Squamous.cells.de_3Mvs1Y.df$DiffExpressed !="NO"]
# Define the coloring of the graph
mycolors <- c("red","blue","grey")
names(mycolors) <- c("UP","DOWN","NO")



# Volcano plot
Squamouscells3Mvs1Y <- bulk.Squamous.cells.de_3Mvs1Y.df%>%
  #filter(p_val<0.2) %>%
  #is.na() %>%
  ggplot(aes(avg_log2FC,-log10(p_val),
             col=DiffExpressed,
             label=diffLabel))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  labs(x="Log2 Fold change",
       y="-log10(P Value)",
       title = "Squamous cells",
       subtitle = "HIV+ ART<3 months vs HIV+ ART>1 Year")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        plot.subtitle = element_text(hjust = 0.5, face = "bold",color = "grey"),
        axis.title = element_text(hjust = 0.5,face = "bold",size = 15))
Squamouscells3Mvs1Y <- as.ggplot(Squamouscells3Mvs1Y)
Squamouscells3Mvs1Y
Squamouscells3MvsNeg+Squamouscells1YvsNeg+Squamouscells3Mvs1Y

# save the saved converted dimplot 
ggsave("results/SquamousVolcanoPlots.pdf",
       plot = ((Squamouscells3MvsNeg|Squamouscells1YvsNeg|Squamouscells3Mvs1Y|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 20, height = 8, unit="in", dpi = 700,limitsize = F)


# Plot violins of significant genes----
VlnPlot(subset(all_merged_subset_labelled_final,downsample=200),
        features = c("NKG7","LAG3","CD74","HLA-DRB1","GZMH","CCL4","HLA-A","IFITM1","GZMB","STAT1",
                     "CD79A","MS4A1","IGKV3-20","CDV3","HLA-DRB1","HLA-DQA1","HLA-A","BANK1","ANXA1",
                     "RAMP1","CXCL2","IDO1","DNAH5","VNN3","MUC13","TRIM31","SYNE1","MAP2",
                     "IL32","PTGS2","TMCC3","IFI6","VNN2","PLAT","SLC7A2","GNGT1","RND1","IP6K3"),
        idents = c("CD3+ T cells_HIV-",
                   "CD3+ T cells_HIV+ ART<3 Months",
                   "CD3+ T cells_HIV+ ART>1 Year",
                   "B cells_HIV-",
                   "B cells_HIV+ ART<3 Months",
                   "B cells_HIV+ ART>1 Year",
                   "Goblet cells_HIV-",
                   "Goblet cells_HIV+ ART<3 Months",
                   "Goblet cells_HIV+ ART>1 Year",
                   "Secretory cells_HIV-",
                   "Secretory cells_HIV+ ART<3 Months",
                   "Secretory cells_HIV+ ART>1 Year"),
                   group.by="HIV_Status",
        stack = T,
        log = F,
        raster = T,
        same.y.lims = T)+
  NoLegend()


DotPlot(subset(all_merged_subset_labelled_final,downsample=200),
        features = c("NKG7","LAG3","CD74","HLA-DRB1","GZMH","CCL4","HLA-A","IFITM1","GZMB","STAT1",
                     "CD79A","MS4A1","IGKV3-20","CDV3","HLA-DRB1","HLA-DQA1","HLA-A","BANK1","ANXA1",
                     "RAMP1","CXCL2","IDO1","DNAH5","VNN3","MUC13","TRIM31","SYNE1","MAP2",
                     "IL32","PTGS2","TMCC3","IFI6","VNN2","PLAT","SLC7A2","GNGT1","RND1","IP6K3"),
        idents = c("CD3+ T cells_HIV-",
                   "CD3+ T cells_HIV+ ART<3 Months",
                   "CD3+ T cells_HIV+ ART>1 Year",
                   "B cells_HIV-",
                   "B cells_HIV+ ART<3 Months",
                   "B cells_HIV+ ART>1 Year",
                   "Goblet cells_HIV-",
                   "Goblet cells_HIV+ ART<3 Months",
                   "Goblet cells_HIV+ ART>1 Year",
                   "Secretory cells_HIV-",
                   "Secretory cells_HIV+ ART<3 Months",
                   "Secretory cells_HIV+ ART>1 Year"),
        group.by="HIV_Status")+
  NoLegend()



bulk.Neutrophils.cells.de_3MvsNegsig <- bulk.Neutrophils.cells.de_3MvsNeg%>%
  filter(p_val>0.05,avg_log2FC>1.5)
bulk.Neutrophils.cells.de_1YvsNegsig <- bulk.Neutrophils.cells.de_1YvsNeg%>%
  filter(p_val<0.05,avg_log2FC>1.5)

DoHeatmap(all_merged_subset_labelled_final,
          features = c("NLRP3","FFAR2","PPP1R15A","CD14","SPAG9","MARCKS","STX11","CD69","NFKBIA","NFKBIZ",
                       "FCGBP","CCDC146","CRIP1","RARRES1","MT-ATP8","ERCC1","TUBA1A","MUC5AC",
                       "AGR2","VAMP8","VMO1","ADIRF","SLPI","EPAS1","CD177","CST1","LGALS7B","CPA4"), 
          assay = "RNA", #slot = "counts", 
          angle = 45,size = 4, group.by = "HIV_Status") + 
  scale_fill_viridis_c() #+ NoLegend()




# GENE MODULE IDENTIFICATION
# Identification of Gene modules in HIV infection in T cells cells----

Tcells_res <- unique(rbind(bulk.T.cells.de_3MvsNeg.df,
                           bulk.T.cells.de_1YvsNeg.df,
                           bulk.T.cells.de_3Mvs1Y.df)) %>%
  filter(p_val<0.05)

T_HIV_rlog_out <- rlog(Tcellsdds, blind=F) #get normalized count data from dds object
T_HIV_mat<-assay(T_HIV_rlog_out) %>%
  as.data.frame()
T_HIV_mat <- rownames(T_HIV_mat %in% Tcells_res$genes)
T_HIV_mat<-assay(T_HIV_rlog_out)[rownames(Tcells_res), rownames(ColData)] #sig genes x samples
T_HIV_mat_cem <- as.data.frame(T_HIV_mat)

sample_annot <- read.csv("data/cemitool_HIV_sample_annot.csv",
                         header = TRUE#,sep = ","
)

cem <-  (T_HIV_mat_cem,
                sample_annot,
                force_beta = T)
#sample_annotation(cem, 
#sample_name_column="SampleName", 
#class_column="Class") <- sample_annot

#sample_annot <- colnames(Goblet_HIV_mat_cem)

nmodules(cem)
cem <- mod_gsea(cem)
cem <- plot_gsea(cem)
show_plot(cem, "gsea")

cem <- plot_profile(cem)
plots <- show_plot(cem,"profile")
plots[2]

# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)

cem <- mod_ora(cem, gmt_in)
cem <- plot_ora(cem)
plots <- show_plot(cem, "ora")
plots[2]

# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)

# plot interactions
library(ggplot2)
interactions_data(cem) <- int_df # add interactions
cem <- plot_interactions(cem) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module
plots[4]
plots[1]
plots[2]
plots[3]

# run cemitool
library(ggplot2)
cem <- cemitool(T_HIV_mat_cem, sample_annot, gmt_in, interactions=int_df,
                force_beta = T,
                filter=TRUE, plot=TRUE, verbose=TRUE)
# create report as html document
generate_report(cem, directory="./results/immuneReport")

# write analysis results into files
write_files(cem, directory="./results/immuneTables")

# save all plots
save_plots(cem, "all", directory="./results/immunePlots")




# VENN DIAGRAMS----
# T cells ----
# T cells genes differentially expressed in HIV+ ART<3 Months
bulk.T.cells.de_3M1.df_diff <- bulk.T.cells.de_3MvsNeg %>%
  filter(avg_log2FC>1.5,p_val<0.05)
bulk.T.cells.de_3M3.df_diff <- bulk.T.cells.de_3Mvs1Y %>%
  filter(avg_log2FC>1.5,p_val<0.05)
Diffin3M.Tcells <- dplyr::combine(bulk.T.cells.de_3M1.df_diff,
                           bulk.T.cells.de_3M3.df_diff)


# T cells genes differentially expressed in HIV+ ART>1 year
bulk.T.cells.de_1Y2.df_diff <- bulk.T.cells.de_1YvsNeg %>%
  filter(avg_log2FC>1.5,p_val<0.05)
bulk.T.cells.de.1Y3.df.diff <- bulk.T.cells.de_3Mvs1Y %>%
  filter(avg_log2FC<1.5,p_val<0.05)
diffin1Y.Tcells <- dplyr::combine(bulk.T.cells.de_1Y2.df_diff,
                                  bulk.T.cells.de.1Y3.df.diff)


# T cells genes differentially expressed in HIV-
bulk.T.cells.de_Neg1.df_diff <- bulk.T.cells.de_3MvsNeg %>%
  filter(avg_log2FC<1.5,p_val<0.05)
bulk.T.cells.de_Neg2.df_diff <- bulk.T.cells.de_1YvsNeg %>%
  filter(avg_log2FC<1.5,p_val<0.05)
diffinNeg.Tcells <- dplyr::combine(bulk.T.cells.de_Neg1.df_diff,
                                   bulk.T.cells.de_Neg2.df_diff)

# T cells venn diagram
venn.diagram(x=list(rownames(Diffin3M.Tcells),
                    rownames(diffin1Y.Tcells),
                    rownames(diffinNeg.Tcells)
                    ),
             category.names = c('PLHIV ART<3 Months',
                                'PLHIV ART>1 Year',
                                'HIV-'),
             filename = 'Thesis Figures/Tcells_venn_diagram.png',
             output=TRUE,
             imagetype='png',
             height = 480,
             width = 480,
             resolution = 300,
             compression = 'lwz',
             lwd=1,
             col=c('red','blue','green'),
             fill=c(alpha('red',0.3),alpha('blue',0.3),alpha('blue',0.3)),
             cex=0.5,
             fontfamily='sans',
             cat.cex=0.3,
             cat.default.pos='outer',
             cat.pos=c(-15,15,180),
             cat.dist=c(0.055,0.055,0.045),
             cat.fontfamily='sans',
             cat.col=c('red','blue','green'),
             rotation=1)


# selecting genes unique to each element in a venn diagram


Tcelloverlap <- calculate.overlap(x=list(rownames(Diffin3M.Tcells),
                                   rownames(diffin1Y.Tcells),
                                   rownames(diffinNeg.Tcells)))

unique.Tcell.genes.in.3M <- Tcelloverlap$a1 %>%
  as.data.frame()

Diffin3M.Tcells$genes <- rownames(Diffin3M.Tcells)
unique.Tcell.genes.in.3M.df <- Diffin3M.Tcells %>%
  filter(genes %in% unique.Tcell.genes.in.3M$.)

unique.Tcell.genes.in.1Y <- Tcelloverlap$a3 %>%
  as.data.frame()

diffin1Y.Tcells$genes <- rownames(diffin1Y.Tcells)
unique.Tcell.genes.in.1Y.df <- diffin1Y.Tcells %>%
  filter(genes %in% unique.Tcell.genes.in.1Y$.)

unique.Tcell.genes.in.Neg <- Tcelloverlap$a7 %>%
  as.data.frame()

diffinNeg.Tcells$genes <- rownames(diffinNeg.Tcells)
unique.Tcell.genes.in.Neg.df <- diffinNeg.Tcells %>%
  filter(genes %in% unique.Tcell.genes.in.Neg$.)

combined.unique.genes.in.Tcells <- combine(unique.Tcell.genes.in.3M.df,
                                           unique.Tcell.genes.in.1Y.df,
                                           unique.Tcell.genes.in.Neg.df)

# B cells ----
# B cells genes differentially expressed in HIV+ ART<3 Months
bulk.B.cells.de_3M1.df_diff <- bulk.B.cells.de_3MvsNeg %>%
  filter(avg_log2FC>1.5,p_val<0.05)
bulk.B.cells.de_3M3.df_diff <- bulk.B.cells.de_3Mvs1Y %>%
  filter(avg_log2FC>1.5,p_val<0.05)
Diffin3M.Bcells <- combine(bulk.B.cells.de_3M1.df_diff,
                           bulk.B.cells.de_3M3.df_diff)


# B cells genes differentially expressed in HIV+ ART>1 year
bulk.B.cells.de_1Y2.df_diff <- bulk.B.cells.de_1YvsNeg %>%
  filter(avg_log2FC>1.5,p_val<0.05)
bulk.B.cells.de.1Y3.df.diff <- bulk.B.cells.de_3Mvs1Y %>%
  filter(avg_log2FC<1.5,p_val<0.05)
diffin1Y.Bcells <- combine(bulk.B.cells.de_1Y2.df_diff,
                           bulk.B.cells.de.1Y3.df.diff)


# B cells genes differentially expressed in HIV-
bulk.B.cells.de_Neg1.df_diff <- bulk.B.cells.de_3MvsNeg %>%
  filter(avg_log2FC<1.5,p_val<0.05)
bulk.B.cells.de_Neg2.df_diff <- bulk.B.cells.de_1YvsNeg %>%
  filter(avg_log2FC<1.5,p_val<0.05)
diffinNeg.Bcells <- combine(bulk.B.cells.de_Neg1.df_diff,
                            bulk.B.cells.de_Neg2.df_diff)

# B cells venn diagram
venn.diagram(x=list(rownames(Diffin3M.Bcells),
                    rownames(diffin1Y.Bcells),
                    rownames(diffinNeg.Bcells)
),
category.names = c('HIV+ ART<3 Months',
                   'HIV+ ART>1 Year',
                   'HIV-'),
filename = 'results/Bcells_venn_diagram.png',
output=TRUE,
imagetype='png',
height = 480,
width = 480,
resolution = 300,
compression = 'lwz',
lwd=1,
col=c('red','blue','green'),
fill=c(alpha('red',0.3),alpha('blue',0.3),alpha('blue',0.3)),
cex=0.5,
fontfamily='sans',
cat.cex=0.3,
cat.default.pos='outer',
cat.pos=c(-15,15,180),
cat.dist=c(0.055,0.055,0.045),
cat.fontfamily='sans',
cat.col=c('red','blue','green'),
rotation=1)


# selecting genes unique to each element in a venn diagram


Bcelloverlap <- calculate.overlap(x=list(rownames(Diffin3M.Bcells),
                                         rownames(diffin1Y.Bcells),
                                         rownames(diffinNeg.Bcells)))

unique.Bcell.genes.in.3M <- Bcelloverlap$a1 %>%
  as.data.frame()

Diffin3M.Bcells$genes <- rownames(Diffin3M.Bcells)
unique.Bcell.genes.in.3M.df <- Diffin3M.Bcells %>%
  filter(genes %in% unique.Bcell.genes.in.3M$.)

unique.Bcell.genes.in.1Y <- Bcelloverlap$a3 %>%
  as.data.frame()

diffin1Y.Bcells$genes <- rownames(diffin1Y.Bcells)
unique.Bcell.genes.in.1Y.df <- diffin1Y.Bcells %>%
  filter(genes %in% unique.Bcell.genes.in.1Y$.)

unique.Bcell.genes.in.Neg <- Bcelloverlap$a7 %>%
  as.data.frame()

diffinNeg.Bcells$genes <- rownames(diffinNeg.Bcells)
unique.Bcell.genes.in.Neg.df <- diffinNeg.Bcells %>%
  filter(genes %in% unique.Bcell.genes.in.Neg$.)

combined.unique.genes.in.Bcells <- combine(unique.Bcell.genes.in.3M.df,
                                           unique.Bcell.genes.in.1Y.df,
                                           unique.Bcell.genes.in.Neg.df)


# Neutrophils ----
  # Neutrophil genes differentially expressed in HIV+ ART<3 Months
bulk.Neutrophils.cells.de_3M1.df_diff <- bulk.Neutrophils.cells.de_3MvsNeg %>%
  filter(avg_log2FC>=1.5,p_val<0.05)
bulk.Neutrophils.cells.de_3M3.df_diff <- bulk.Neutrophils.cells.de_3Mvs1Y %>%
  filter(avg_log2FC>=1.5,p_val<0.05)
Diffin3M.Neutrophils <- combine(bulk.Neutrophils.cells.de_3M1.df_diff,
                           bulk.Neutrophils.cells.de_3M3.df_diff)

# Neutrophil genes differentially expressed in HIV+ ART>1 year
bulk.Neutrophils.cells.de_1Y2.df_diff <- bulk.Neutrophils.cells.de_1YvsNeg %>%
  filter(avg_log2FC>=1.5,p_val<0.05)
bulk.Neutrophils.cells.de.1Y3.df.diff <- bulk.Neutrophils.cells.de_3Mvs1Y %>%
  filter(avg_log2FC<1.5,p_val<0.05)
diffin1Y.Neutrophils <- combine(bulk.Neutrophils.cells.de_1Y2.df_diff,
                           bulk.Neutrophils.cells.de.1Y3.df.diff)

# Neurtophil genes differentially expressed in HIV-
bulk.Neutrophils.cells.de_Neg1.df_diff <- bulk.Neutrophils.cells.de_3MvsNeg %>%
  filter(avg_log2FC<1.5,p_val<0.05)
bulk.Neutrophils.cells.de_Neg2.df_diff <- bulk.Neutrophils.cells.de_1YvsNeg %>%
  filter(avg_log2FC<1.5,p_val<0.05)
diffinNeg.Neutrophils <- combine(bulk.Neutrophils.cells.de_Neg1.df_diff,
                            bulk.Neutrophils.cells.de_Neg2.df_diff)

# Neutrophil venn diagram
venn.diagram(x=list(rownames(Diffin3M.Neutrophils),
                    rownames(diffin1Y.Neutrophils),
                    rownames(diffinNeg.Neutrophils)),
             category.names = c('HIV+ ART<3 Months',
                                'HIV+ ART>1 Year',
                                'HIV-'),
             filename = 'results/Neutrophils_venn_diagram.png',
             output=TRUE,
             imagetype='png',
             height = 480,
             width = 480,
             resolution = 300,
             compression = 'lwz',
             lwd=1,
             col=c('red','blue','green'),
             fill=c(alpha('red',0.3),alpha('blue',0.3),alpha('blue',0.3)),
             cex=0.5,
             fontfamily='sans',
             cat.cex=0.3,
             cat.default.pos='outer',
             cat.pos=c(-15,15,180),
             cat.dist=c(0.055,0.055,0.045),
             cat.fontfamily='sans',
             cat.col=c('red','blue','green'),
             rotation=1)


# selecting genes unique to each element in a venn diagram


neutrophilsoverlap <- calculate.overlap(x=list(rownames(Diffin3M.Neutrophils),
                                    rownames(diffin1Y.Neutrophils),
                                    rownames(diffinNeg.Neutrophils)))

unique.neutrophil.genes.in.3M <- neutrophilsoverlap$a1 %>%
  as.data.frame()

Diffin3M.Neutrophils$genes <- rownames(Diffin3M.Neutrophils)
unique.neutrophil.genes.in.3M.df <- Diffin3M.Neutrophils %>%
  filter(genes %in% unique.neutrophil.genes.in.3M$.)


unique.neutrophil.genes.in.1Y <- neutrophilsoverlap$a3 %>%
  as.data.frame()

diffin1Y.Neutrophils$genes <- rownames(diffin1Y.Neutrophils)
unique.neutrophil.genes.in.1Y.df <- diffin1Y.Neutrophils %>%
  filter(genes %in% unique.neutrophil.genes.in.1Y$.)


unique.neutrophil.genes.in.Neg <- neutrophilsoverlap$a7 %>%
  as.data.frame()

diffinNeg.Neutrophils$genes <- rownames(diffinNeg.Neutrophils)
unique.neutrophil.genes.in.Neg.df <- diffinNeg.Neutrophils %>%
  filter(genes %in% unique.neutrophil.genes.in.Neg$.)

combined.unique.genes.in.neutrophils <- combine(unique.neutrophil.genes.in.3M.df,
                                                unique.neutrophil.genes.in.1Y.df,
                                                unique.neutrophil.genes.in.Neg.df)


# Goblet cells ----
# Goblet genes differentially expressed in HIV+ ART<3 Months
bulk.Goblet.cells.de_3M1.df_diff <- bulk.Goblet.cells.de_3MvsNeg %>%
  filter(avg_log2FC>=1.5,p_val<0.05)
bulk.Goblet.cells.de_3M3.df_diff <- bulk.Goblet.cells.de_3Mvs1Y %>%
  filter(avg_log2FC>=1.5,p_val<0.05)
Diffin3M.Goblet <- combine(bulk.Goblet.cells.de_3M1.df_diff,
                           bulk.Goblet.cells.de_3M3.df_diff)

# Goblet genes differentially expressed in HIV+ ART>1 year
bulk.Goblet.cells.de_1Y2.df_diff <- bulk.Goblet.cells.de_1YvsNeg %>%
  filter(avg_log2FC>=1.5,p_val<0.05)
bulk.Goblet.cells.de.1Y3.df.diff <- bulk.Goblet.cells.de_3Mvs1Y %>%
  filter(avg_log2FC<1.5,p_val<0.05)
diffin1Y.Goblet <- combine(bulk.Goblet.cells.de_1Y2.df_diff,
                           bulk.Goblet.cells.de.1Y3.df.diff)

# Goblet genes differentially expressed in HIV-
bulk.Goblet.cells.de_Neg1.df_diff <- bulk.Goblet.cells.de_3MvsNeg %>%
  filter(avg_log2FC<1.5,p_val<0.05)
bulk.Goblet.cells.de_Neg2.df_diff <- bulk.Goblet.cells.de_1YvsNeg %>%
  filter(avg_log2FC<1.5,p_val<0.05)
diffinNeg.Goblet <- combine(bulk.Goblet.cells.de_Neg1.df_diff,
                            bulk.Goblet.cells.de_Neg2.df_diff)

diffinNeg.Goblet$genes <- rownames(diffinNeg.Goblet)
diffinNeg.Goblet<- diffinNeg.Goblet%>%
  filter(genes!="WARS",genes!="PARP14")

view(Gobletoverlap$a5)


# Goblet venn diagram
venn.diagram(x=list(rownames(Diffin3M.Goblet),
                    rownames(diffin1Y.Goblet),
                    rownames(diffinNeg.Goblet)),
             category.names = c('HIV+ ART<3 Months',
                                'HIV+ ART>1 Year',
                                'HIV-'),
             filename = 'results/Goblet_venn_diagram.png',
             output=TRUE,
             imagetype='png',
             height = 480,
             width = 480,
             resolution = 300,
             compression = 'lwz',
             lwd=1,
             col=c('red','blue','green'),
             fill=c(alpha('red',0.3),alpha('blue',0.3),alpha('blue',0.3)),
             cex=0.5,
             fontfamily='sans',
             cat.cex=0.3,
             cat.default.pos='outer',
             cat.pos=c(-15,15,180),
             cat.dist=c(0.055,0.055,0.045),
             cat.fontfamily='sans',
             cat.col=c('red','blue','green'),
             rotation=1)


# selecting genes unique to each element in a venn diagram

Gobletoverlap <- calculate.overlap(x=list(rownames(Diffin3M.Goblet),
                                               rownames(diffin1Y.Goblet),
                                               rownames(diffinNeg.Goblet)))

unique.Goblet.genes.in.3M <- Gobletoverlap$a1 %>%
  as.data.frame()

Diffin3M.Goblet$genes <- rownames(Diffin3M.Goblet)
unique.Goblet.genes.in.3M.df <- Diffin3M.Goblet %>%
  filter(genes %in% unique.Goblet.genes.in.3M$.)


unique.Goblet.genes.in.1Y <- Gobletoverlap$a3 %>%
  as.data.frame()

diffin1Y.Goblet$genes <- rownames(diffin1Y.Goblet)
unique.Goblet.genes.in.1Y.df <- diffin1Y.Goblet %>%
  filter(genes %in% unique.Goblet.genes.in.1Y$.)


unique.Goblet.genes.in.Neg <- Gobletoverlap$a7 %>%
  as.data.frame()

diffinNeg.Goblet$genes <- rownames(diffinNeg.Goblet)
unique.Goblet.genes.in.Neg.df <- diffinNeg.Goblet %>%
  filter(genes %in% unique.Goblet.genes.in.Neg$.)

combined.unique.genes.in.Goblet <- combine(unique.Goblet.genes.in.3M.df,
                                           unique.Goblet.genes.in.1Y.df,
                                           unique.Goblet.genes.in.Neg.df)

# Secretory cells ----
# Secretory genes differentially expressed in HIV+ ART<3 Months
bulk.Secretory.cells.de_3M1.df_diff <- bulk.Secretory.cells.de_3MvsNeg %>%
  filter(avg_log2FC>=1.5,p_val<0.05)
bulk.Secretory.cells.de_3M3.df_diff <- bulk.Secretory.cells.de_3Mvs1Y %>%
  filter(avg_log2FC>=1.5,p_val<0.05)
Diffin3M.Secretory <- combine(bulk.Secretory.cells.de_3M1.df_diff,
                              bulk.Secretory.cells.de_3M3.df_diff)


# Secretory genes differentially expressed in HIV+ ART>1 year
bulk.Secretory.cells.de_1Y2.df_diff <- bulk.Secretory.cells.de_1YvsNeg %>%
  filter(avg_log2FC>=1.5,p_val<0.05)
bulk.Secretory.cells.de.1Y3.df.diff <- bulk.Secretory.cells.de_3Mvs1Y %>%
  filter(avg_log2FC<1.5,p_val<0.05)
diffin1Y.Secretory <- combine(bulk.Secretory.cells.de_1Y2.df_diff,
                           bulk.Secretory.cells.de.1Y3.df.diff
                           )

# Secretory genes differentially expressed in HIV-
bulk.Secretory.cells.de_Neg1.df_diff <- bulk.Secretory.cells.de_3MvsNeg %>%
  filter(avg_log2FC<1.5,p_val<0.05)
bulk.Secretory.cells.de_Neg2.df_diff <- bulk.Secretory.cells.de_1YvsNeg %>%
  filter(avg_log2FC<1.5,p_val<0.05)
diffinNeg.Secretory <- combine(bulk.Secretory.cells.de_Neg1.df_diff,
                            bulk.Secretory.cells.de_Neg2.df_diff) 

diffinNeg.Secretory$genes <- rownames(diffinNeg.Secretory)
diffinNeg.Secretory <- diffinNeg.Secretory%>%
  filter(genes!="DUSP5",genes!="IFI44L")



# Secretory venn diagram
venn.diagram(x=list(rownames(Diffin3M.Secretory),
                    rownames(diffin1Y.Secretory),
                    rownames(diffinNeg.Secretory)),
             category.names = c('HIV+ ART<3 Months',
                                'HIV+ ART>1 Year',
                                'HIV-'),
             filename = 'results/Secretory_venn_diagram.png',
             output=TRUE,
             imagetype='png',
             height = 480,
             width = 480,
             resolution = 300,
             compression = 'lwz',
             lwd=1,
             col=c('red','blue','green'),
             fill=c(alpha('red',0.3),alpha('blue',0.3),alpha('blue',0.3)),
             cex=0.5,
             fontfamily='sans',
             cat.cex=0.3,
             cat.default.pos='outer',
             cat.pos=c(-15,15,180),
             cat.dist=c(0.055,0.055,0.045),
             cat.fontfamily='sans',
             cat.col=c('red','blue','green'),
             rotation=1)


# selecting genes unique to each element in a venn diagram

Secretoryoverlap <- calculate.overlap(x=list(rownames(Diffin3M.Secretory),
                                          rownames(diffin1Y.Secretory),
                                          rownames(diffinNeg.Secretory)))

unique.Secretory.genes.in.3M <- Secretoryoverlap$a5 %>%
  as.data.frame()

Diffin3M.Secretory$genes <- rownames(Diffin3M.Secretory)
unique.Secretory.genes.in.3M.df <- Diffin3M.Secretory %>%
  filter(genes %in% unique.Secretory.genes.in.3M$.)


unique.Secretory.genes.in.1Y <- Secretoryoverlap$a3 %>%
  as.data.frame()

diffin1Y.Secretory$genes <- rownames(diffin1Y.Secretory)
unique.Secretory.genes.in.1Y.df <- diffin1Y.Secretory %>%
  filter(genes %in% unique.Secretory.genes.in.1Y$.)


unique.Secretory.genes.in.Neg <- Secretoryoverlap$a7 %>%
  as.data.frame()

diffinNeg.Secretory$genes <- rownames(diffinNeg.Secretory)
unique.Secretory.genes.in.Neg.df <- diffinNeg.Secretory %>%
  filter(genes %in% unique.Secretory.genes.in.Neg$.)

combined.unique.genes.in.Secretory <- combine(unique.Secretory.genes.in.3M.df,
                                           unique.Secretory.genes.in.1Y.df,
                                           unique.Secretory.genes.in.Neg.df)

# Dendritic cells ----
# Dendritic cells genes differentially expressed in HIV+ ART<3 Months
bulk.Dendritic.cells.de_3M1.df_diff <- bulk.Dendritic.cells.de_3MvsNeg %>%
  filter(avg_log2FC>1.5,p_val<0.05)
bulk.Dendritic.cells.de_3M3.df_diff <- bulk.Dendritic.cells.de_3Mvs1Y %>%
  filter(avg_log2FC>1.5,p_val<0.05)
Diffin3M.Dendriticcells <- dplyr::combine(bulk.Dendritic.cells.de_3M1.df_diff,
                                          bulk.Dendritic.cells.de_3M3.df_diff)


# Dendritic cells genes differentially expressed in HIV+ ART>1 year
bulk.Dendritic.cells.de_1Y2.df_diff <- bulk.Dendritic.cells.de_1YvsNeg %>%
  filter(avg_log2FC>1.5,p_val<0.05)
bulk.Dendritic.cells.de.1Y3.df.diff <- bulk.Dendritic.cells.de_3Mvs1Y %>%
  filter(avg_log2FC<1.5,p_val<0.05)
diffin1Y.Dendriticcells <- dplyr::combine(bulk.Dendritic.cells.de_1Y2.df_diff,
                                          bulk.Dendritic.cells.de.1Y3.df.diff)


# Dendritic cells genes differentially expressed in HIV-
bulk.Dendritic.cells.de_Neg1.df_diff <- bulk.Dendritic.cells.de_3MvsNeg %>%
  filter(avg_log2FC<1.5,p_val<0.05)
bulk.Dendritic.cells.de_Neg2.df_diff <- bulk.Dendritic.cells.de_1YvsNeg %>%
  filter(avg_log2FC<1.5,p_val<0.05)
diffinNeg.Dendriticcells <- dplyr::combine(bulk.Dendritic.cells.de_Neg1.df_diff,
                                           bulk.Dendritic.cells.de_Neg2.df_diff)

# T cells venn diagram
venn.diagram(x=list(rownames(Diffin3M.Dendriticcells),
                    rownames(diffin1Y.Dendriticcells),
                    rownames(diffinNeg.Dendriticcells)
),
category.names = c('HIV+ ART<3 Months',
                   'HIV+ ART>1 Year',
                   'HIV-'),
filename = 'results/Dendriticcells_venn_diagram.png',
output=TRUE,
imagetype='png',
height = 480,
width = 480,
resolution = 300,
compression = 'lwz',
lwd=1,
col=c('red','blue','green'),
fill=c(alpha('red',0.3),alpha('blue',0.3),alpha('blue',0.3)),
cex=0.5,
fontfamily='sans',
cat.cex=0.3,
cat.default.pos='outer',
cat.pos=c(-15,15,180),
cat.dist=c(0.055,0.055,0.045),
cat.fontfamily='sans',
cat.col=c('red','blue','green'),
rotation=1)


# selecting genes unique to each element in a venn diagram


Dendriticcelloverlap <- calculate.overlap(x=list(rownames(Diffin3M.Dendriticcells),
                                                 rownames(diffin1Y.Dendriticcells),
                                                 rownames(diffinNeg.Dendriticcells)))

unique.Dendriticcell.genes.in.3M <- Dendriticcelloverlap$a1 %>%
  as.data.frame()

Diffin3M.Dendriticcells$genes <- rownames(Diffin3M.Dendriticcells)
unique.Dendriticcell.genes.in.3M.df <- Diffin3M.Dendriticcells %>%
  filter(genes %in% unique.Dendriticcell.genes.in.3M$.)

unique.Dendriticcell.genes.in.1Y <- Dendriticcelloverlap$a3 %>%
  as.data.frame()

diffin1Y.Dendriticcells$genes <- rownames(diffin1Y.Dendriticcells)
unique.Dendriticcell.genes.in.1Y.df <- diffin1Y.Dendriticcells %>%
  filter(genes %in% unique.Dendriticcell.genes.in.1Y$.)

unique.Dendriticcell.genes.in.Neg <- Dendriticcelloverlap$a7 %>%
  as.data.frame()

diffinNeg.Dendriticcells$genes <- rownames(diffinNeg.Dendriticcells)
unique.Dendriticcell.genes.in.Neg.df <- diffinNeg.Dendriticcells %>%
  filter(genes %in% unique.Dendriticcell.genes.in.Neg$.)

combined.unique.genes.in.Dendriticcells <- combine(unique.Dendriticcell.genes.in.3M.df,
                                                   unique.Dendriticcell.genes.in.1Y.df,
                                                   unique.Dendriticcell.genes.in.Neg.df)

# GSEA----
# Prepare input data----

genelistGSEA <- combined.unique.genes.in.Dendriticcells
# we want the log2 fold change---- 
original_gene_list <- genelistGSEA$avg_log2FC

# name the vector----
names(original_gene_list) <- combined.unique.genes.in.Dendriticcells$genes

# omit any NA values----
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for cluster Profiler)----
gene_list = sort(gene_list, decreasing = TRUE)

head(gene_list)


gse <- gseGO(geneList = gene_list,
             ont = "BP",
             keyType = "SYMBOL",
             minGSSize = 3,
             maxGSSize =800,
             pvalueCutoff = 0.1,
             verbose = T,
             OrgDb = "org.Hs.eg.db",
             pAdjustMethod = "none")

view(gse@result)
print()

# Dotplot----
require(DOSE)
dotplot(gse, showCategory=20, split=".sign")+
  facet_grid(.~.sign)

enrichplot::gseaplot2(gse,geneSetID = c(1))

heatplot(gse, showCategory = 10)
heatplot(gse, foldChange=gene_list, showCategory=10)

gse2<-pairwise_termsim(gse, showCategory = 20)
treeplot(gse2)
treeplot(gse2,hclust_method="average")
upsetplot(gse)
upsetplot(kk2)

# Encrichment Map----
emapplot(gse2, showCategory=20)


tailgse <- gse@result %>%
  filter(setSize<5)

gse@result <- tailgse
# Category Netplot----
cnetplot(gse,categ3orySize="pvalue", foldChange=gene_list,
         showCategory=30, circular=T, colorEdge=T)

# GSEA Plot----
gseaplot(gse, by="runningScore", title=gse$Description[72], geneSetID = 10)


# Enrich GO----
df <- combined.unique.genes.in.Dendriticcells
original_gene_list <- df$avg_log2FC
names(original_gene_list) <- rownames(df)
gene_list <- na.omit(original_gene_list)
gene_list <- sort(gene_list,decreasing = TRUE)
sig_genes_df <- subset(df,p_val<0.05)
genes <- sig_genes_df$avg_log2FC
names(genes) <- rownames(sig_genes_df)
genes <- na.omit(genes)
genes <- names(genes)[abs(genes)>1.5]

organism <- "org.Hs.eg.db"
go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "BP",
                      minGSSize = 5,
                      pAdjustMethod = "none",
                      maxGSSize = 500#,
                      #pvalueCutoff = 0.1, 
                      #qvalueCutoff = 0.10
                      )

print(go_enrich)

# Upset plot : emphasized the genes overlapping among different sets
upsetplot(go_enrich)
# Barplot
GO3MvsNeg <- barplot(go_enrich,
                     drop=TRUE,
                     showCategory=20,
                     title="Go Biological pathaways",
                     font.size=5,
                     p.adjust.methods="none")+
  theme(plot.title = element_text(hjust = 0.5,face = 'bold',size = 12),
        axis.text.x = element_text(hjust = 0.5,face = 'bold',size = 8),
        axis.text.y = element_text(hjust = 0.5,face = 'bold'))
GO3MvsNeg


go_enrich_data <- go_enrich@result
  ggplot(go_enrich,aes(x=go_enrich@result$Description,y=go_enrich@result$Count))+
    geom_col()

cnetplot(go_enrich,
         categorySize="pvalue",
         foldChange=gene_list,
         font.size=1)


enrichplot::gseaplot2(gse2,geneSetID = c(6:8))
enrichplot::upsetplot(gse2)

# save the saved converted dimplot 
ggsave("results/GO3MvsNeg.pdf",
       plot = ((GO3MvsNeg|plot_layout(ncol = 3, nrow = 1, width = c(1,1,1,1,1,1,1,1)))),
       width = 10, height = 6, unit="in", dpi = 700,limitsize = F)


# CEMITOOL FOR MODULE IDENTIFICATION----

# Genes differentially expressed in T cells----
sample_annot <- read.csv("data/cemitool_HIV_all_sample_annot.csv",
                         header = TRUE#,sep = ","
)

# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)

# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)

# Preparing expression matrix
load('data/Tcellsdds.RData')
Tcells_HIV_rlog_out <- rlog(Tcellsdds, blind=FALSE) #get normalized count data from dds object
#Neutrophils_HIV_rlog_out <- Neutrophil_cells_dds_HIV
Tcells_HIV_mat<-assay(Tcells_HIV_rlog_out)#sig genes x samples
Tcells_HIV_mat <- na.omit(Tcells_HIV_mat)
Tcells_HIV_mat <- as.data.frame(Tcells_HIV_mat)
Tcells_HIV_mat <- Tcells_HIV_mat #%>%
#filter(rownames(Tcells_HIV_mat) %in% combined.unique.genes.in.Tcells$genes)

# Create a CEMITOOL object
Tcellscem <- cemitool(Tcells_HIV_mat,
                annot = sample_annot,
                gmt = gmt_in,
                interactions = int_df,
                filter = TRUE,
                #filter_pval = 0.05,
                #apply_vst = T,
                cor_method = "pearson",
                cor_function = "cor",
                network_type = "signed",
                force_beta = T,
                ora_pval = 0.05,
                min_ngen = 5,
                gsea_min_size = 3,
                gsea_max_size = 800,
                plot = T,
                plot_diagnostics = T)

# Plot gene set enrichment analysis modules
Tcellsmodule_plots <- show_plot(Tcellscem,"gsea")
Tcellsmodule_plots

# Plot over representation analysis
Tcellsora_plots <- show_plot(Tcellscem, "ora")
TcellsOraM1<-Tcellsora_plots[1]
TcellsOraM1
TcellsOraM2<-Tcellsora_plots[2]
TcellsOraM2
TcellsOraM3<-Tcellsora_plots[3]
TcellsOraM3
TcellsOraM4<-Tcellsora_plots[4]
TcellsOraM4

# plot interactions
Tcellsinteractions_plots <- show_plot(Tcellscem, "interaction")#+
 #geom_text_repel(max.overlaps = 10) # view the plot for the first module
TcellsInteractionsM1<-Tcellsinteractions_plots[1]
TcellsInteractionsM1
TcellsInteractionsM4<-Tcellsinteractions_plots[4]
TcellsInteractionsM4

# create report as html document
generate_report(Tcellscem, directory="./results/TcellReport",force=T)

# write analysis results into files
write_files(Tcellscem, directory="./results/TcellReport",force=T)

# save all plots
save_plots(Tcellscem, "all", directory="./results/TcellReport",force=T)

# Genes differentially expressed in Dendritic cells----
sample_annot <- read.csv("data/cemitool_HIV_all_sample_annot.csv",
                         header = TRUE#,sep = ","
)

# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)

# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)

# Preparing expression matrix
Dendriticcells_HIV_rlog_out <- rlog(Dendriticdds, blind=FALSE) #get normalized count data from dds object
#Neutrophils_HIV_rlog_out <- Neutrophil_cells_dds_HIV
Dendriticcells_HIV_mat<-assay(Dendriticcells_HIV_rlog_out)#sig genes x samples
Dendriticcells_HIV_mat <- na.omit(Dendriticcells_HIV_mat)
Dendriticcells_HIV_mat <- as.data.frame(Dendriticcells_HIV_mat)
Dendriticcells_HIV_mat <- Dendriticcells_HIV_mat #%>%
#filter(rownames(Tcells_HIV_mat) %in% combined.unique.genes.in.Tcells$genes)

# Create a CEMITOOL object
Dendriticcellscem <- cemitool(Dendriticcells_HIV_mat,
                      annot = sample_annot,
                      gmt = gmt_in,
                      interactions = int_df,
                      filter = TRUE,
                      #filter_pval = 0.05,
                      #apply_vst = T,
                      cor_method = "pearson",
                      cor_function = "cor",
                      network_type = "signed",
                      force_beta = T,
                      ora_pval = 0.05,
                      min_ngen = 20,
                      gsea_min_size = 3,
                      gsea_max_size = 800,
                      plot = T,
                      plot_diagnostics = T)

# Plot gene set enrichment analysis modules
Dendriticcellsmodule_plots <- show_plot(Dendriticcellscem,"gsea")
Dendriticcellsmodule_plots

# Plot over representation analysis
Dendriticcellsora_plots <- show_plot(Dendriticcellscem, "ora")
DendriticcellsOraM1<-Dendriticcellsora_plots[1]
DendriticcellsOraM1
DendriticcellsOraM3<-Dendriticcellsora_plots[3]
DendriticcellsOraM3
DendriticcellsOraM5<-Dendriticcellsora_plots[5]
DendriticcellsOraM5

# plot interactions
Dendriticcellsinteractions_plots <- show_plot(Dendriticcellscem, "interaction") # view the plot for the first module
DendriticcellsInteractionsM1<-Dendriticcellsinteractions_plots[1]
DendriticcellsInteractionsM1
DendriticcellsInteractionsM3<-Dendriticcellsinteractions_plots[3]
DendriticcellsInteractionsM3
DendriticcellsInteractionsM5<-Dendriticcellsinteractions_plots[5]
DendriticcellsInteractionsM5
# create report as html document
generate_report(Dendriticcellscem, directory="./results/DendriticcellReport",force=T)

# write analysis results into files
write_files(Dendriticcellscem, directory="./results/DendriticcellReport",force=T)

# save all plots
save_plots(Dendriticcellscem, "all", directory="./results/DendriticcellReport",force=T)
# Genes differentially expressed in B cells----
sample_annot <- read.csv("data/cemitool_HIV_all_sample_annot.csv",
                         header = TRUE#,sep = ","
)

# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)

# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)

# Preparing expression matrix
Bcells_HIV_rlog_out <- rlog(Bcellsdds, blind=FALSE) #get normalized count data from dds object
#Neutrophils_HIV_rlog_out <- Neutrophil_cells_dds_HIV
Bcells_HIV_mat<-assay(Bcells_HIV_rlog_out)#sig genes x samples
Bcells_HIV_mat <- na.omit(Bcells_HIV_mat)
Bcells_HIV_mat <- as.data.frame(Bcells_HIV_mat)
Bcells_HIV_mat <- Bcells_HIV_mat #%>%
#filter(rownames(Tcells_HIV_mat) %in% combined.unique.genes.in.Tcells$genes)

# Create a CEMITOOL object
Bcellscem <- cemitool(Bcells_HIV_mat,
                      annot = sample_annot,
                      gmt = gmt_in,
                      interactions = int_df,
                      filter = TRUE,
                      filter_pval = 0.05,
                      #apply_vst = T,
                      cor_method = "pearson",
                      cor_function = "cor",
                      network_type = "signed",
                      force_beta = T,
                      ora_pval = 0.05,
                      min_ngen = 5,
                      gsea_min_size = 3,
                      gsea_max_size = 800,
                      plot = T,
                      plot_diagnostics = T)

# Plot gene set enrichment analysis modules
Bcellsmodule_plots <- show_plot(Bcellscem,"gsea")
Bcellsmodule_plots

# Plot over representation analysis
Bcellsora_plots <- show_plot(Bcellscem, "ora")
BcellsOraM1<-Bcellsora_plots[1]
BcellsOraM1
BcellsOraM2<-Bcellsora_plots[2]
BcellsOraM2
BcellsOraM3<-Bcellsora_plots[3]
BcellsOraM3

# plot interactions
Bcellsinteractions_plots <- show_plot(Bcellscem, "interaction") # view the plot for the first module
BcellsInteractionsM1<-Bcellsinteractions_plots[1]
BcellsInteractionsM3<-Bcellsinteractions_plots[3]

# create report as html document
generate_report(Bcellscem, directory="./results/BcellReport",force=T)

# write analysis results into files
write_files(Tcellscem, directory="./results/BcellReport",force=T)

# save all plots
save_plots(Tcellscem, "all", directory="./results/TcellReport",force=T)

# Genes differentially expressed in Neutrophils ----

sample_annot <- read.csv("data/cemitool_HIV_all_sample_annot.csv",
                         header = TRUE#,sep = ","
)

# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)

# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)

# Create a CEMITOOL object
Neutrophils_HIV_rlog_out <- rlog(Neutrophil_dds, blind=FALSE) #get normalized count data from dds object
#Neutrophils_HIV_rlog_out <- Neutrophil_cells_dds_HIV
Neutrophils_HIV_mat<-assay(Neutrophil_dds)#sig genes x samples
Neutrophils_HIV_mat <- na.omit(Neutrophils_HIV_mat)
Neutrophils_HIV_mat <- as.data.frame(Neutrophils_HIV_mat)
#Neutrophils_HIV_mat <- Neutrophils_HIV_mat %>%
  #filter(rownames(Neutrophils_HIV_mat) %in% combined.unique.genes.in.neutrophils$genes)


Neutrophilscem <- cemitool(Neutrophils_HIV_mat,
                      annot = sample_annot,
                      gmt = gmt_in,
                      interactions = int_df,
                      filter = TRUE,
                      filter_pval = 0.05,
                      #apply_vst = T,
                      cor_method = "pearson",
                      cor_function = "cor",
                      network_type = "signed",
                      gsea_scale = T,
                      force_beta = T,
                      ora_pval = 0.05,
                      min_ngen = 5,
                      gsea_min_size = 3,
                      gsea_max_size = 800,
                      plot = T,
                      center_func = 'median',
                      plot_diagnostics = T)

# Plot gene set enrichment analysis modules
Neutrophilsmodule_plots <- show_plot(Neutrophilscem,"gsea")
Neutrophilsmodule_plots

# Plot over representation analysis
Neutrophilsora_plots <- show_plot(Neutrophilscem, "ora")
NeutrophilsOraM1<-Neutrophilsora_plots[1]
NeutrophilsOraM1
NeutrophilsOraM2<-Neutrophilsora_plots[2]
NeutrophilsOraM2
M2genes <- as.data.frame(Neutrophilscem@ora) %>%
  filter(Module=="M2",Description=="Activated TLR4 signalling")

NeutrophilsOraM3<-Neutrophilsora_plots[3]
NeutrophilsOraM3

# plot interactions
Neutrophilsinteractions_plots <- show_plot(Neutrophilscem, "interaction") # view the plot for the first module
NeutrophilsInteractionsM1<-Neutrophilsinteractions_plots[1]
NeutrophilsInteractionsM1
NeutrophilsInteractionsM2<-Neutrophilsinteractions_plots[2]
NeutrophilsInteractionsM2
NeutrophilsInteractionsM3<-Neutrophilsinteractions_plots[3]
NeutrophilsInteractionsM3

# create report as html document
generate_report(Neutrophilscem, directory="./results/NeutrophilsReport",force=T)

# write analysis results into files
write_files(Neutrophilscem, directory="./results/NeutrophilsReport",force=T)

# save all plots
save_plots(Neutrophilscem, "all", directory="./results/NeutrophilsReport",force=T)


# Genes differentially expressed in Goblet ----

sample_annot <- read.csv("data/cemitool_HIV_all_sample_annot.csv",
                         header = TRUE#,sep = ","
)

# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)

# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)

# Create a CEMITOOL object
Goblet_rlog_out <- rlog(Gobletdds, blind=FALSE) #get normalized count data from dds object
Goblet_mat<-assay(Goblet_rlog_out)#sig genes x samples
Goblet_mat <- na.omit(Goblet_mat)
Goblet_mat <- as.data.frame(Goblet_mat)
#Neutrophils_HIV_mat <- Neutrophils_HIV_mat #%>%
  #filter(rownames(Neutrophils_HIV_mat) %in% combined.unique.genes.in.neutrophils$genes)


Gobletcem <- cemitool(Goblet_mat,
                      annot = sample_annot,
                      gmt = gmt_in,
                      interactions = int_df,
                      filter = TRUE,
                      filter_pval = 0.05,
                      #apply_vst = T,
                      cor_method = "pearson",
                      cor_function = "cor",
                      network_type = "signed",
                      gsea_scale = T,
                      force_beta = T,
                      ora_pval = 0.05,
                      min_ngen = 30,
                      gsea_min_size = 5,
                      gsea_max_size = 800,
                      plot = T,
                      center_func = 'mean',
                      plot_diagnostics = T)

# Plot gene set enrichment analysis modules
Gobletmodule_plots <- show_plot(Gobletcem,"gsea")
Gobletmodule_plots

# Plot over representation analysis
Gobletora_plots <- show_plot(Gobletcem, "ora")
GobletOraM1<-Gobletora_plots[1]
GobletOraM1
GobletOraM2<-Gobletora_plots[2]
GobletOraM2
GobletOraM3<-Gobletora_plots[3]
GobletOraM3
GobletOraM4<-Gobletora_plots[4]
GobletOraM4

# plot interactions
Gobletinteractions_plots <- show_plot(Gobletcem, "interaction") # view the plot for the first module
GobletInteractionsM1<-Gobletinteractions_plots[1]
GobletInteractionsM1
GobletInteractionsM3<-Gobletinteractions_plots[3]
GobletInteractionsM3

# create report as html document
generate_report(Gobletcem, directory="./results/GobletReport",force=T)

# write analysis results into files
write_files(Gobletcem, directory="./results/GobletReport",force=T)

# save all plots
save_plots(Gobletcem, "all", directory="./results/GobletReport",force=T)

# Genes differentially expressed in Secretory cells ----

sample_annot <- read.csv("data/cemitool_HIV_all_sample_annot.csv",
                         header = TRUE#,sep = ","
)

# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)

# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)

# Create a CEMITOOL object
Secretory_rlog_out <- rlog(Secretorydds, blind=FALSE) #get normalized count data from dds object
Secretory_mat<-assay(Secretory_rlog_out)#sig genes x samples
Secretory_mat <- na.omit(Secretory_mat)
Secretory_mat <- as.data.frame(Secretory_mat)
#Neutrophils_HIV_mat <- Neutrophils_HIV_mat #%>%
#filter(rownames(Neutrophils_HIV_mat) %in% combined.unique.genes.in.neutrophils$genes)


Secretorycem <- cemitool(Secretory_mat,
                      annot = sample_annot,
                      gmt = gmt_in,
                      interactions = int_df,
                      filter = TRUE,
                      filter_pval = 0.05,
                      #apply_vst = T,
                      cor_method = "spearman",
                      cor_function = "cor",
                      network_type = "signed",
                      gsea_scale = T,
                      force_beta = T,
                      ora_pval = 0.05,
                      min_ngen = 20,
                      gsea_min_size = 5,
                      gsea_max_size = 800,
                      plot = T,
                      center_func = 'mean',
                      plot_diagnostics = T)

# Plot gene set enrichment analysis modules
Secretorymodule_plots <- show_plot(Secretorycem,"gsea")
Secretorymodule_plots

# Plot over representation analysis
Secretoryora_plots <- show_plot(Secretorycem, "ora")
SecretoryoraM1<-Secretoryora_plots[1]
SecretoryoraM1
SecretoryoraM4<-Secretoryora_plots[4]
SecretoryoraM4

# plot interactions
Secretoryinteractions_plots <- show_plot(Secretorycem, "interaction") # view the plot for the first module
SecretoryinteractionsM1<-Secretoryinteractions_plots[1]
SecretoryinteractionsM1
SecretoryinteractionsM3<-Secretoryinteractions_plots[3]
SecretoryinteractionsM3

# create report as html document
generate_report(Secretorycem, directory="./results/SecretoryReport",force=T)

# write analysis results into files
write_files(Secretorycem, directory="./results/SecretoryReport",force=T)

# save all plots
save_plots(Secretorycem, "all", directory="./results/SecretoryReport",force=T)


# Genes differentially expressed in Basal cells ----

sample_annot <- read.csv("data/cemitool_HIV_all_sample_annot.csv",
                         header = TRUE#,sep = ","
)

# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)

# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)

# Create a CEMITOOL object
Basal_rlog_out <- rlog(Basaldds, blind=FALSE) #get normalized count data from dds object
Basal_mat <- assay(Basal_rlog_out)#sig genes x samples
#Basal_mat <- Basal_mat[row.names(Basal_mat) %in% rownames(Basalres),]
Basal_mat <- na.omit(Basal_mat)
Basal_mat <- as.data.frame(Basal_mat)
#Neutrophils_HIV_mat <- Neutrophils_HIV_mat #%>%
#filter(rownames(Neutrophils_HIV_mat) %in% combined.unique.genes.in.neutrophils$genes)


Basalcem <- cemitool(Basal_mat,
                         annot = sample_annot,
                         gmt = gmt_in,
                         interactions = int_df,
                         filter = TRUE,
                         filter_pval = 0.05,
                         #apply_vst = T,
                         cor_method = "spearman",
                         cor_function = "cor",
                         network_type = "signed",
                         gsea_scale = T,
                         force_beta = T,
                         ora_pval = 0.05,
                         min_ngen = 30,
                         gsea_min_size = 10,
                         gsea_max_size = 800,
                         plot = T,
                         center_func = 'mean',
                         plot_diagnostics = T)

# Plot gene set enrichment analysis modules
Basalmodule_plots <- show_plot(Basalcem,"gsea")
Basalmodule_plots

# Plot over representation analysis
Basalora_plots <- show_plot(Basalcem, "ora")
BasaloraM1<-Basalora_plots[1]
BasaloraM1
BasaloraM5<-Basalora_plots[5]
BasaloraM5

# plot interactions
Basalinteractions_plots <- show_plot(Basalcem, "interaction") # view the plot for the first module
BasalinteractionsM1<-Basalinteractions_plots[1]
BasalinteractionsM1
BasalinteractionsM5<-Basalinteractions_plots[5]
BasalinteractionsM5

# create report as html document
generate_report(Basalcem, directory="./results/BasalReport",force=T)

# write analysis results into files
write_files(Basalcem, directory="./results/BasalReport",force=T)

# save all plots
save_plots(Basalcem, "all", directory="./results/BasalReport",force=T)


