
# T cell Module Analysis using CEMiTool
# Author: [Joseph Aston Phiri]
# Date: [22nd December 2023]
# Description: This script performs module enrichment analysis on CD3+ T cells using CEMiTool.

# Required Packages
library(Seurat)
library(dplyr)
library(stringr)
library(tibble)
library(CEMiTool)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggpubr)
library(readr)
library(scales)

# Load the annotated Seurat object
all_merged_subset_labelled_new <- readRDS('Data/Single_Cell_Data/all_merged_subset_labelled_new.rds')
all_merged_subset_labelled_new$Clusters <- as.character(all_merged_subset_labelled_new@active.ident)
Idents(all_merged_subset_labelled_new) <- all_merged_subset_labelled_new$Clusters

# Subset for CD3+ T cells
Tcells <- subset(all_merged_subset_labelled_new, idents = "CD3+ T cells")

# Aggregate gene expression by HIV status and sample
Tcell_cts <- Seurat::AggregateExpression(
  Tcells,
  group.by = c("HIV_Status", "sample"),
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)$RNA %>% as.data.frame()

# Filter out low-expressed genes
Tcell_cts <- Tcell_cts[rowSums(Tcell_cts) > 3, ]

# Create sample-level metadata
colData <- data.frame(samples = colnames(Tcell_cts)) %>%
  mutate(HIV_Status = case_when(
    grepl("HIV-", samples) ~ "HIV-",
    grepl("ART<3", samples) ~ "PLHIV on ART<3m",
    TRUE ~ "PLHIV on ART>1y"
  ),
  Sample_Name = str_extract(samples, "(?<=_).*")) %>%
  column_to_rownames("samples")

# Write annotation to CSV
sample_annot <- colData %>%
  rownames_to_column("SampleName") %>%
  mutate(Class = HIV_Status) %>%
  select(SampleName, Class)
write.csv(sample_annot, "scRNAseq_Results/Sample_annotation.csv", row.names = FALSE)

# Load gene set and interaction references
gmt_in <- read_gmt(system.file("extdata", "pathways.gmt", package = "CEMiTool"))
int_df <- read.delim(system.file("extdata", "interactions.tsv", package = "CEMiTool"))

# Run CEMiTool
Tcell_cem <- cemitool(
  expr = Tcell_cts,
  annot = sample_annot,
  gmt = gmt_in,
  interactions = int_df,
  filter = TRUE,
  filter_pval = 0.05,
  apply_vst = TRUE,
  cor_method = "spearman",
  cor_function = "cor",
  network_type = "signed",
  gsea_scale = TRUE,
  force_beta = TRUE,
  ora_pval = 0.05,
  min_ngen = 5,
  gsea_min_size = 3,
  gsea_max_size = 800,
  plot = TRUE,
  center_func = "median",
  plot_diagnostics = TRUE,
  verbose = TRUE
)

# Save and reload object
saveRDS(Tcell_cem, "scRNAseq_Results/Tcell_cem.rds")
Tcell_cem <- readRDS("scRNAseq_Results/Tcell_cem.rds")

# GSEA plots
gsea <- show_plot(Tcell_cem, "gsea")
ora <- show_plot(Tcell_cem, "ora")

# Enrichment matrices
es <- Tcell_cem@enrichment$es %>%
  pivot_longer(cols = c("HIV-","PLHIV on ART<3m","PLHIV on ART>1y"), names_to = "HIV_Status", values_to = "ES")
nes <- Tcell_cem@enrichment$nes %>%
  pivot_longer(cols = c("HIV-","PLHIV on ART<3m","PLHIV on ART>1y"),, names_to = "HIV_Status", values_to = "NES")
padj <- Tcell_cem@enrichment$padj %>%
  pivot_longer(cols = c("HIV-","PLHIV on ART<3m","PLHIV on ART>1y"),, names_to = "HIV_Status", values_to = "adj_pvalue")

Tcell_enrichment <- full_join(nes, es, by = c("pathway", "HIV_Status")) %>%
  full_join(padj, by = c("pathway", "HIV_Status")) %>%
  drop_na() %>%
  mutate(adj_pvalue = -log10(adj_pvalue))

# Dot Plot: T cell enriched modules
Fig7d_new <- ggplot(Tcell_enrichment, aes(x = HIV_Status, y = pathway)) +
  geom_point(aes(color = NES, size = adj_pvalue)) +
  scale_color_gradient2(low = "#0000B4", mid = "white", high = "#931500") +
  labs(title = "d.", subtitle = "T cell enriched modules") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 60, face = "bold", hjust = -0.1),
    plot.subtitle = element_text(size = 30, hjust = 0.5),
    axis.text = element_text(size = 30),
    axis.title = element_blank(),
    axis.ticks.length = unit(0.5, 'cm')
  )

ggsave("New_Figures/Figure 7/Fig7d_new.png", Fig7d_new, width = 10, height = 10, dpi = 300)
ggsave("New_Figures/Figure 7/Fig7d_new.pdf", Fig7d_new, width = 10, height = 10, dpi = 300)

# Additional figures for ORA can be added similarly
