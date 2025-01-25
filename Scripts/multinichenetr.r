library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
library(nichenetr)
library(Seurat)

load("Immune_cells.RData")
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
Immune_cells$Clusters <- paste0(Immune_cells@active.ident)
Immune_cells$HIV_Status <- paste0(Immune_cells$HIV_Status)
Immune_cells$sample <- paste0(Immune_cells$sample)
# Make names syntactically valid
Immune_cells$HIV_Status <- make.names(Immune_cells$HIV_Status)
Immune_cells$sample <- make.names(Immune_cells$sample)
Immune_cells$Clusters <- make.names(Immune_cells$Clusters)
# Convert to a single cell experiment
sce <- Seurat::as.SingleCellExperiment(Immune_cells)
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

absent_celltypes <- abundance_df_summarized %>%
  dplyr::filter(samples_present < 2) %>%
  dplyr::group_by(celltype_id) %>%
  dplyr::summarize(n = n(), .groups = "drop") %>%
  dplyr::filter(n == total_nr_conditions) %>%
  pull(celltype_id)

analyse_condition_specific_celltypes = TRUE
if(analyse_condition_specific_celltypes == TRUE){
  senders_oi = senders_oi %>% dplyr::setdiff(absent_celltypes)
  receivers_oi = receivers_oi %>% dplyr::setdiff(absent_celltypes)
} else {
  senders_oi = senders_oi %>%
    dplyr::setdiff(union(absent_celltypes, condition_specific_celltypes))
  receivers_oi = receivers_oi %>%
    dplyr::setdiff(union(absent_celltypes, condition_specific_celltypes))
}

sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in%
            c(senders_oi, receivers_oi)
]

min_sample_prop = 0.50
fraction_cutoff = 0.05

frq_list = get_frac_exprs(
  sce = sce,
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id,
  batches = batches,
  min_cells = min_cells,
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)

genes_oi = frq_list$expressed_df %>%
  filter(expressed == TRUE) %>% pull(gene) %>% unique()
sce = sce[genes_oi, ]

# Pseudobulk expression calculation
abundance_expression_info = process_abundance_expression_info(
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
  batches = batches, covariates = covariates,
  contrasts_oi = contrasts_oi,
  min_cells = min_cells,
  expressed_df = frq_list$expressed_df)
# Check DE results
DE_info$celltype_de$de_output_tidy %>% head()
DE_info$hist_pvals
DE_info$celltype_de
DE_info$hist_pvals_findmarkers
DE_info$celltype_de_findmarkers


empirical_pval = FALSE
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
  celltype_de = DE_info_emp$de_output_tidy_emp %>% select(-p_val, -p_adj) %>%
    rename(p_val = p_emp, p_adj = p_adj_emp)
} else {
  celltype_de = DE_info$celltype_de$de_output_tidy
}

# Combine DE information for ligand-senders and receptors-receivers
sender_receiver_de = multinichenetr::combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network)
sender_receiver_de %>% head(20)

# Ligand activity prediction: use DE analysis output to predict the activity
# of ligands in receiver cell types and infer their potential targer genes
logFC_threshold = 0.250
p_val_threshold = 0.1
p_val_adj = FALSE
geneset_assessment = contrast_tbl$contrast %>%
  lapply(
    process_geneset_data,
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>%
  bind_rows()
geneset_assessment

# In case i want to use adjusted p values
geneset_assessment_adjustedPval = contrast_tbl$contrast %>%
  lapply(
    process_geneset_data,
    celltype_de, logFC_threshold, p_val_adj = TRUE, p_val_threshold
  ) %>%
  bind_rows()
geneset_assessment_adjustedPval

# Perform the ligand activity analysis and ligand-target inference
top_n_target = 250

verbose = TRUE
cores_system = 1
n.cores = min(cores_system, celltype_de$cluster_id %>% unique() %>% length())

# Running the ligand activity prediction
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

# Prioritization: Rank cell-cell communication patters through multi-criteria prioritization

ligand_activity_down = FALSE
sender_receiver_tbl = sender_receiver_de %>% distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>%
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>%
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(multinichenetr::generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular", # all prioritization criteria will be weighted equally
  fraction_cutoff = fraction_cutoff,
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down))
# Check output tables
prioritization_tables$group_prioritization_tbl %>% head(20)

# Calculate cross-samples expression correlation between ligand-receptor pairs and target genes
lr_target_prior_cor = lr_target_prior_cor_inference(
  receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(),
  abundance_expression_info = abundance_expression_info,
  celltype_de = celltype_de,
  grouping_tbl = grouping_tbl,
  prioritization_tables = prioritization_tables,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj)

# Save the output of multinichenetr
Path = "scRNAseq_Results/"
Immune_multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor)
Immune_multinichenet_output = multinichenetr::make_lite_output(Immune_multinichenet_output)

save = TRUE
if(save == TRUE){
  saveRDS(Immune_multinichenet_output, paste0(Path, "Immune_multinichenet_output.rds"))
  
}

Immune_multinichenet_output <- readRDS("scRNAseq_Results/Immune_multinichenet_output.rds")

# Visualization of differential cell-cell interactions
prioritized_tbl_oi_all = get_top_n_lr_pairs(
  Immune_multinichenet_output$prioritization_tables,
  top_n = 50,
  rank_per_group = T)

prioritized_tbl_oi =
  Immune_multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  dplyr::filter(id %in% prioritized_tbl_oi_all$id) %>%
  dplyr::distinct(id, sender, receiver, ligand, receptor, group) %>%
  dplyr::left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)

# Save the output of Circos_list
HIV_neg_circos <- circos_list$HIV.
HIV_1Y_circos <- circos_list$HIV..ART.1.Year
HIV_3M_circos <- circos_list$HIV..ART.3.Months

library(gridGraphics)
recordedplot_to_grob <- function(recorded_plot) {
  gridGraphics::grid.echo(recorded_plot)
  grid::grid.grab()
}

HIV_neg_circos <- recordedplot_to_grob(HIV_neg_circos)
HIV_3M_circos <- recordedplot_to_grob(HIV_3M_circos)
HIV_1Y_circos <- recordedplot_to_grob(HIV_1Y_circos)


library(ggplot2)

# Example: Convert grob to ggplot
ggplot_from_grob <- function(grob) {
  ggplot() +
    annotation_custom(grob) +
    theme_void()  # Remove axes, grid lines, etc.
}

# Convert your grob
HIV_neg_circos <- ggplot_from_grob(HIV_neg_circos)
HIV_3M_circos <- ggplot_from_grob(HIV_3M_circos)
HIV_1Y_circos <- ggplot_from_grob(HIV_1Y_circos)

# Combine Multiple ggplot Objects
library(patchwork)
combined_plot <- HIV_neg_circos + HIV_3M_circos + HIV_1Y_circos + plot_layout(ncol = 3)
print(combined_plot)

# Save Fig4
ggsave(combined_plot,filename="Figures/Figure 4/CElltoCell.png",
width = 45,height = 15,dpi = 1080)
ggsave(combined_plot,filename="Figures/Figure 4/CElltoCell.pdf",
width = 45,height = 15,dpi = 1080)


circos_list$legend
