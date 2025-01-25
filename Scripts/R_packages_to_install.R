# R packages to install on the MLW SAPITWA server

# Seurat version 5
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)

# SingleCellExperiment
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")

# Multinichenetr
install.packages("devtools")
devtools::install_github("saeyslab/nichenetr")
devtools::install_github("saeyslab/multinichenetr")

# CellChat
devtools::install_github("jinworks/CellChat")

# CEMiTool
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("CEMiTool")

# SoupX
install.packages('SoupX')

# Organism database
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")