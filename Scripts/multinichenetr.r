# Load Required Libraries
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
library(nichenetr)
library(Seurat)
library(gridGraphics)
library(patchwork)

#------------------------#
# 1. Load Input Data
#------------------------#

seurat_obj <- readRDS("Data/Single_Cell_Data/all_merged_subset_labelled_new.rds")
organism <- "human"

# Ligand-Receptor Network
lr_network <- readRDS("NicheNet/lr_network_human_allInfo_30112033.rds") %>%
  mutate(
    ligand = convert_alias_to_symbols(ligand, organism),
    receptor = convert_alias_to_symbols(receptor, organism),
    ligand = make.names(ligand),
    receptor = make.names(receptor)
  ) %>%
  distinct(ligand, receptor)

# Ligand-Target Matrix
ligand_target_matrix <- readRDS("NicheNet/ligand_target_matrix_nsga2r_final.rds")
colnames(ligand_target_matrix) <- make.names(convert_alias_to_symbols(colnames(ligand_target_matrix), organism))
rownames(ligand_target_matrix) <- make.names(convert_alias_to_symbols(rownames(ligand_target_matrix), organism))

# Filter LR Network based on ligand presence
lr_network <- lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
ligand_target_matrix <- ligand_target_matrix[, unique(lr_network$ligand)]

#------------------------#
# 2. Preprocessing Seurat Object
#------------------------#

seurat_obj$Clusters <- make.names(seurat_obj@active.ident)
seurat_obj$HIV_Status <- make.names(seurat_obj$HIV_Status)
seurat_obj$sample <- make.names(seurat_obj$sample)

sce <- as.SingleCellExperiment(seurat_obj)

#------------------------#
# 3. Experimental Metadata Setup
#------------------------#

sample_id <- "sample"
group_id <- "HIV_Status"
celltype_id <- "Clusters"
covariates <- NA
batches <- NA

contrasts <- c(
  "HIV..ART.3.Months-(HIV..ART.1.Year+HIV.)/2",
  "HIV..ART.1.Year-(HIV..ART.3.Months+HIV.)/2",
  "HIV.-(HIV..ART.1.Year+HIV..ART.3.Months)/2"
)
contrast_tbl <- tibble(contrast = contrasts, group = c("HIV..ART.3.Months", "HIV..ART.1.Year", "HIV."))

senders <- unique(colData(sce)[, celltype_id])
receivers <- senders

# Keep only selected conditions
conditions_keep <- c("HIV..ART.3.Months", "HIV..ART.1.Year", "HIV.")
sce <- sce[, colData(sce)[, group_id] %in% conditions_keep]

#------------------------#
# 4. Abundance Filtering
#------------------------#

min_cells <- 3

abundance_info <- get_abundance_info(
  sce, sample_id, group_id, celltype_id,
  min_cells, senders, receivers, batches
)

# Filter cell types by presence
abundance_df_summary <- abundance_info$abundance_data %>%
  mutate(keep = as.logical(keep)) %>%
  group_by(group_id, celltype_id) %>%
  summarise(samples_present = sum(keep), .groups = "drop")

absent_celltypes <- abundance_df_summary %>%
  filter(samples_present < 2) %>%
  group_by(celltype_id) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n == length(unique(colData(sce)[, group_id]))) %>%
  pull(celltype_id)

senders <- setdiff(senders, absent_celltypes)
receivers <- setdiff(receivers, absent_celltypes)

sce <- sce[, colData(sce)[, celltype_id] %in% c(senders, receivers)]

#------------------------#
# 5. Expression Fraction Calculation
#------------------------#

min_sample_prop <- 0.50
fraction_cutoff <- 0.05

frq_list <- get_frac_exprs(
  sce, sample_id, celltype_id, group_id,
  batches, min_cells, fraction_cutoff, min_sample_prop
)

genes_expressed <- frq_list$expressed_df %>%
  filter(expressed) %>%
  pull(gene) %>% unique()

sce <- sce[genes_expressed, ]

#------------------------#
# 6. Pseudobulk + DE Analysis
#------------------------#

ab_expr_info <- process_abundance_expression_info(
  sce, sample_id, group_id, celltype_id, min_cells,
  senders, receivers, lr_network, batches, frq_list, abundance_info
)

DE_info <- get_DE_info(
  sce, sample_id, group_id, celltype_id,
  batches, covariates, contrasts, min_cells,
  frq_list$expressed_df
)

celltype_de <- DE_info$celltype_de$de_output_tidy

#------------------------#
# 7. Sender-Receiver DE Integration
#------------------------#

sender_receiver_de <- combine_sender_receiver_de(
  sender_de = celltype_de, receiver_de = celltype_de,
  senders_oi = senders, receivers_oi = receivers,
  lr_network = lr_network
)

#------------------------#
# 8. Ligand Activity & Prioritization
#------------------------#

logFC_threshold <- 0.5
p_val_threshold <- 0.05
p_val_adj <- FALSE
top_n_target <- 250
verbose <- TRUE
cores_system <- 1
n_cores <- min(cores_system, length(unique(celltype_de$cluster_id)))

geneset_assessment <- contrast_tbl$contrast %>%
  lapply(process_geneset_data, celltype_de, logFC_threshold, p_val_adj, p_val_threshold) %>%
  bind_rows()

ligand_activities <- get_ligand_activities_targets_DEgenes(
  receiver_de = celltype_de,
  receivers_oi = intersect(receivers, unique(celltype_de$cluster_id)),
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj,
  top_n_target = top_n_target,
  verbose = verbose,
  n.cores = n_cores
)

#------------------------#
# 9. Prioritization Tables
#------------------------#

sender_receiver_tbl <- distinct(sender_receiver_de, sender, receiver)

grouping_tbl <- as_tibble(colData(sce)) %>%
  select(all_of(c(sample_id, group_id))) %>%
  distinct() %>%
  rename(sample = !!sample_id, group = !!group_id)

prioritization_tbls <- generate_prioritization_tables(
  sender_receiver_info = ab_expr_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular",
  fraction_cutoff = fraction_cutoff,
  abundance_data_receiver = ab_expr_info$abundance_data_receiver,
  abundance_data_sender = ab_expr_info$abundance_data_sender,
  ligand_activity_down = FALSE
)

#------------------------#
# 10. Ligand-Receptor-Target Correlation
#------------------------#

lr_target_prior_cor <- lr_target_prior_cor_inference(
  receivers_oi = unique(prioritization_tbls$group_prioritization_tbl$receiver),
  abundance_expression_info = ab_expr_info,
  celltype_de = celltype_de,
  grouping_tbl = grouping_tbl,
  prioritization_tables = prioritization_tbls,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj
)

#------------------------#
# 11. Save Output
#------------------------#

output_path <- "scRNAseq_Results/"
output_file <- paste0(output_path, "multinichenet_output.rds")

multinichenet_output <- list(
  celltype_info = ab_expr_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = ab_expr_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities,
  prioritization_tables = prioritization_tbls,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
) %>% make_lite_output()

saveRDS(multinichenet_output, output_file)

#------------------------#
# 12. Visualization
#------------------------#

# Circos plots
plot_top_lr <- function(tbl, top_n = 250) {
  top_tbl <- get_top_n_lr_pairs(multinichenet_output$prioritization_tables, top_n = top_n, rank_per_group = TRUE)
  tbl <- left_join(tbl, top_tbl) %>% mutate(prioritization_score = coalesce(prioritization_score, 0))
  senders_receivers <- union(tbl$sender, tbl$receiver) %>% sort()
  col_s <- RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% set_names(senders_receivers)
  circos_list <- make_circos_group_comparison(tbl, col_s, col_s)
  return(circos_list)
}

prioritized_tbl <- multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  distinct(id, sender, receiver, ligand, receptor, group)

circos_plots <- plot_top_lr(prioritized_tbl)

# Save as ggplot objects
save_circos_ggplot <- function(recorded_plot, file) {
  grob <- recordedplot_to_grob(recorded_plot)
  p <- ggplot() + annotation_custom(grob) + theme_void()
  ggsave(p, filename = file, width = 11, height = 8, dpi = 300)
}

save_circos_ggplot(circos_plots$HIV., "Figures/Figure 4/HIV_negative.png")
save_circos_ggplot(circos_plots$HIV..ART.3.Months, "Figures/Figure 4/HIV_ART_3M.png")
save_circos_ggplot(circos_plots$HIV..ART.1.Year, "Figures/Figure 4/HIV_ART_1Y.png")
