#installing and loading required packages----

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")

# install.packages("devtools")
devtools::install_github("saeyslab/nichenetr")
devtools::install_github("saeyslab/multinichenetr")

install.packages("dplyr")
install.packages("ggplot2")
install.packages("Seurat")

library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
library(nichenetr)
library(Seurat)

# Make gene names syntactically valid using make.names to avoid loss of genes and load the organism name----

options(timeout = 120)

if(organism == "human"){
  
  lr_network_all = readRDS(url("https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds")) %>% 
    mutate(ligand = convert_alias_to_symbols(ligand, organism = organism),
    receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
} else if(organism == "mouse"){
  lr_network_all = readRDS(url(
    "https://zenodo.org/record/10229222/files/lr_network_mouse_allInfo_30112033.rds"
    )) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"
    ))
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]  
}

# Load the seurat object and prepare it for conversion to single cells experiment object----
# Load object ----
load('HIV PAPER/Data/Single_Cell_Data/R Files/all_merged_subset_labelled.RData')
all_merged_subset_labelled_final$cell_cluster <- paste0(all_merged_subset_labelled_final@active.ident) 
DimPlot(all_merged_subset_labelled_final,reduction = 'umap.harmony',label = T,repel = T)+NoLegend()
# Make names syntactically valid----
all_merged_subset_labelled_final$HIV_Status <- make.names(all_merged_subset_labelled_final$HIV_Status)
all_merged_subset_labelled_final$cell_cluster <-make.names(all_merged_subset_labelled_final$cell_cluster)
all_merged_subset_labelled_final$sample <-make.names(all_merged_subset_labelled_final$sample)
view(all_merged_subset_labelled_final@meta.data)
# Convert to single cell experiment----
sce <- as.SingleCellExperiment(all_merged_subset_labelled_final)

# Prepare the cell-cell communication analysis----
sample_id='sample'
group_id='HIV_Status'
celltype_id='cell_cluster'
covariates=NA
batches=NA
# Define contrasts----
contrasts_oi = c("'HIV..ART.3.Months-(HIV..ART.1.Year+HIV.)/2','HIV..ART.1.Year-(HIV..ART.3.Months+HIV.)/2','HIV.-(HIV..ART.1.Year+HIV..ART.3.Months)/2'")

contrast_tbl = tibble(contrast =
                        c("HIV..ART.3.Months-(HIV..ART.1.Year+HIV.)/2","HIV..ART.1.Year-(HIV..ART.3.Months+HIV.)/2","HIV.-(HIV..ART.1.Year+HIV..ART.3.Months)/2"),
                      group = c("HIV..ART.3.Months","HIV..ART.1.Year","HIV."))

senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
]

conditions_keep = c("HIV..ART.3.Months", "HIV..ART.1.Year", "HIV.")
sce = sce[, SummarizedExperiment::colData(sce)[,group_id] %in% 
            conditions_keep
]

min_cells = 3
abundance_info = get_abundance_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  batches = batches
)
abundance_info$abund_plot_sample

abundance_df_summarized = abundance_info$abundance_data %>% 
  mutate(keep = as.logical(keep)) %>% 
  group_by(group_id, celltype_id) %>% 
  summarise(samples_present = sum((keep)))

celltypes_absent_one_condition = abundance_df_summarized %>% 
  filter(samples_present == 0) %>% pull(celltype_id) %>% unique() 
# find truly condition-specific cell types by searching for cell types 
# truely absent in at least one condition

celltypes_present_one_condition = abundance_df_summarized %>% 
  filter(samples_present >= 2) %>% pull(celltype_id) %>% unique() 
# require presence in at least 2 samples of one group so 
# it is really present in at least one condition

condition_specific_celltypes = intersect(
  celltypes_absent_one_condition, 
  celltypes_present_one_condition)

total_nr_conditions = SummarizedExperiment::colData(sce)[,group_id] %>% 
  unique() %>% length() 

absent_celltypes <- abundance_df_summarized %>%
  filter(samples_present < 2) %>% 
  group_by(celltype_id) %>%
  summarize(n = n()) %>% 
  filter(n == total_nr_conditions) %>%
  pull(celltype_id)


print("condition-specific celltypes:")
## [1] "condition-specific celltypes:"
print(condition_specific_celltypes)
## character(0)

print("absent celltypes:")
## [1] "absent celltypes:"
print(absent_celltypes)
## character(0)
analyse_condition_specific_celltypes = FALSE

if(analyse_condition_specific_celltypes == TRUE){
  senders_oi = senders_oi %>% setdiff(absent_celltypes)
  receivers_oi = receivers_oi %>% setdiff(absent_celltypes)
} else {
  senders_oi = senders_oi %>% 
    setdiff(union(absent_celltypes, condition_specific_celltypes))
  receivers_oi = receivers_oi %>% 
    setdiff(union(absent_celltypes, condition_specific_celltypes))
}

sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)]

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

abundance_expression_info = process_abundance_expression_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)

abundance_expression_info$celltype_info$pb_df %>% head()
abundance_expression_info$celltype_info$pb_df_group %>% head()
abundance_expression_info$sender_receiver_info$pb_df %>% head()
abundance_expression_info$sender_receiver_info$pb_df_group %>% head()

DE_info = get_DE_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  batches = batches, covariates = covariates, 
  contrasts_oi = contrasts_oi, 
  min_cells = min_cells, 
  expressed_df = frq_list$expressed_df)
DE_info$celltype_de$de_output_tidy %>% head()
DE_info$hist_pvals

empirical_pval = FALSE
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
  celltype_de = DE_info_emp$de_output_tidy_emp %>% select(-p_val, -p_adj) %>% 
    rename(p_val = p_emp, p_adj = p_adj_emp)
} else {
  celltype_de = DE_info$celltype_de$de_output_tidy
} 


sender_receiver_de = multinichenetr::combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)

sender_receiver_de %>% head(20)

logFC_threshold = 0.5
p_val_threshold = 0.05

p_val_adj = FALSE 

geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment

geneset_assessment_adjustedPval = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj = TRUE, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment_adjustedPval


top_n_target = 250

verbose = TRUE
cores_system = 4
n.cores = min(cores_system, celltype_de$cluster_id %>% unique() %>% length()) 

ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = verbose, 
    n.cores = n.cores
  )
))


ligand_activities_targets_DEgenes$ligand_activities %>% head(20)

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
  ligand_activity_down = ligand_activity_down
))

prioritization_tables$group_prioritization_tbl %>% head(20)
view(prioritization_tables$group_prioritization_tbl)

lr_target_prior_cor = lr_target_prior_cor_inference(
  receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(), 
  abundance_expression_info = abundance_expression_info, 
  celltype_de = celltype_de, 
  grouping_tbl = grouping_tbl, 
  prioritization_tables = prioritization_tables, 
  ligand_target_matrix = ligand_target_matrix, 
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj
)

path = "NicheNet/"

multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
) 
multinichenet_output = make_lite_output(multinichenet_output)

save = TRUE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(path, "multinichenet_output.rds"))
  
}

multinichenet_output <- readRDS("NicheNet/multinichenet_output.rds")

# Interpreting the MultiNicheNet analysis output----
prioritized_tbl_oi_all = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  top_n = 100, 
  rank_per_group = FALSE
)

prioritized_tbl_oi = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Paired') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Paired') %>% magrittr::set_names(senders_receivers)


# Create an expanded palette with 13 colors using colorRampPalette
#colors_sender = colorRampPalette(RColorBrewer::brewer.pal(12, 'Paired'))(13) %>% magrittr::set_names(senders_receivers)
#colors_receiver = colorRampPalette(RColorBrewer::brewer.pal(12, 'Paired'))(13) %>% magrittr::set_names(senders_receivers)

# Check the expanded colors
print(colors_sender)
print(colors_receiver)



circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
# Interpretable bubble plots

group_oi = "HIV."
group_oi = "HIV..ART.3.Months"
group_oi = "HIV..ART.1.Year"

prioritized_tbl_oi_M_50 = prioritized_tbl_oi_all %>% 
  filter(group == group_oi)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_M_50 %>% inner_join(lr_network)
)
plot_oi

# Intercellular regulatory network inference and visualization

lr_target_prior = prioritized_tbl_oi_all %>% inner_join(
  multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>%
    distinct(ligand, target, direction_regulation, contrast) %>% inner_join(contrast_tbl) %>% ungroup() 
) 
lr_target_df = lr_target_prior %>% distinct(group, sender, receiver, ligand, receptor, id, target, direction_regulation) 

lr_target_df %>% filter(target %in% union(lr_network$ligand, lr_network$receptor))

network = infer_intercellular_regulatory_network(lr_target_df, prioritized_tbl_oi_all)
network$links %>% head()

colors_sender["L_T_TIM3._CD38._HLADR."] = "pink" # the  original yellow background with white font is not very readable
network_graph = visualize_network(network, colors_sender)
network_graph$plot
