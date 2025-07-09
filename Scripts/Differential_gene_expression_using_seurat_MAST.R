# Required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(readr)
library(tibble)

# Load data
load("Data/Single_Cell_Data/Immune_cells.RData")

#' Run Seurat::FindMarkers using MAST for two groups within a single cluster
#' @param seurat_obj A Seurat object
#' @param cluster_name The name of the cluster to subset
#' @param ident_1 First identity group for comparison
#' @param ident_2 Second identity group for comparison
#' @param label A short label used in output filenames
run_cluster_DE <- function(seurat_obj, cluster_name, ident_1, ident_2, label) {
  cluster_obj <- subset(seurat_obj, subset = Clusters == cluster_name)

  if (ncol(cluster_obj) < 20) {
    warning(paste("Cluster", cluster_name, "has fewer than 20 cells. Skipping."))
    return(NULL)
  }

  cell_counts <- table(cluster_obj@meta.data$HIV_Status)
  print(paste("Cluster:", cluster_name))
  print(cell_counts)

  if (length(cell_counts) < 3 || any(cell_counts == 0)) {
    warning(paste("Cluster", cluster_name, "missing some HIV_Status levels. Skipping."))
    return(NULL)
  }

  de_results <- FindMarkers(cluster_obj,
                            ident.1 = ident_1,
                            ident.2 = ident_2,
                            group.by = "HIV_Status",
                            test.use = "MAST",
                            logfc.threshold = 0.1,
                            min.pct = 0.05)

  de_results <- de_results %>%
    mutate(cluster = cluster_name,
           adj_pval = p.adjust(p_val, method = "BH")) %>%
    rownames_to_column("gene") %>%
    dplyr::select(cluster, gene, avg_log2FC, p_val, adj_pval, pct.1, pct.2)

  write.csv(de_results, paste0("scRNAseq_Results/Seurat_MAST_DE_Analysis/seurat_mast_", label, "_de_", gsub(" ", "_", cluster_name), ".csv"), row.names = FALSE)
  return(de_results)
}

# Define comparisons
comparisons <- list(
  list(ident_1 = "HIV+ ART<3 Months", ident_2 = "HIV-", label = "3mvsneg"),
  list(ident_1 = "HIV+ ART>1 Year", ident_2 = "HIV-", label = "1yvsneg"),
  list(ident_1 = "HIV+ ART<3 Months", ident_2 = "HIV+ ART>1 Year", label = "3mvs1y")
)

clusters <- unique(Immune_cells@meta.data$Clusters)
print(paste("Found", length(clusters), "clusters:", paste(clusters, collapse = ", ")))

for (comp in comparisons) {
  de_results_list <- lapply(clusters, function(cluster) {
    tryCatch({
      run_cluster_DE(Immune_cells, cluster, comp$ident_1, comp$ident_2, comp$label)
    }, error = function(e) {
      warning(paste("Error in cluster", cluster, ":", e$message))
      return(NULL)
    })
  })
  names(de_results_list) <- clusters
}

# Volcano plot helper function
plot_volcano <- function(csv_path, title, subtitle, filename_prefix) {
  df <- read_csv(csv_path) %>%
    mutate(Significance = ifelse(p_val < .05 & avg_log2FC > .25, "Upregulated",
                                 ifelse(p_val < .05 & avg_log2FC < -.25, "Downregulated", "Not significant")),
           label = ifelse(p_val < .05 & abs(avg_log2FC) > .25, gene, NA))

  fig <- ggplot(df, aes(avg_log2FC, -log10(p_val), color = Significance)) +
    geom_point(size = 2.5, alpha = 0.5) +
    geom_text_repel(aes(label = label), size = 8) +
    geom_hline(yintercept = 1.3, linetype = 'dashed') +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = 'dashed') +
    scale_color_manual(values = c("blue", "grey", "red")) +
    labs(x = expression("Average" ~ "log"[2] ~ "FoldChange"),
         y = expression("-" ~ "log"[10] ~ "adjusted P-value"),
         title = title, subtitle = subtitle) +
    scale_x_continuous(limits = c(-5, 3)) +
    theme_classic() +
    theme(plot.title = element_text(size = 60, face = "bold", hjust = -0.1),
          legend.position = 'none',
          axis.ticks.length = unit(0.5, 'cm'),
          plot.subtitle = element_text(size = 30, hjust = .5),
          axis.text = element_text(size = 30),
          axis.title = element_text(size = 35))

  ggsave(fig, filename = paste0("New_Figures/Figure 6/", filename_prefix, ".png"), width = 14, height = 18, dpi = 300)
  ggsave(fig, filename = paste0("New_Figures/Figure 6/", filename_prefix, ".pdf"), width = 14, height = 18, dpi = 300)
  return(fig)
}

# Generate plots
Fig6b_new <- plot_volcano("scRNAseq_Results/Seurat_MAST_DE_Analysis/seurat_mast_3mvsneg_de_Neutrophils.csv",
                          title = 'b.', subtitle = "PLHIV-ART<3m vs HIV-", filename_prefix = "Fig6b_new")
Fig6a_new <- plot_volcano("scRNAseq_Results/Seurat_MAST_DE_Analysis/seurat_mast_1yvsneg_de_Neutrophils.csv",
                          title = 'a.', subtitle = "PLHIV-ART>1yr vs HIV-", filename_prefix = "Fig6a_new")
Fig6c_new <- plot_volcano("scRNAseq_Results/Seurat_MAST_DE_Analysis/seurat_mast_3mvs1y_de_Neutrophils.csv",
                          title = 'c.', subtitle = "PLHIV-ART<3m vs PLHIV-ART>1yr", filename_prefix = "Fig6c_new")
