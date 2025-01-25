
# May you please install the following R packages.
# Thay are important for my analysis and they are memory intensive such that they cant run
# on a loal computer.


# Here are a list of the packages

## Required packages to be installed before SCENIC installation
BiocManager::install(c("AUCell", "RcisTarget"))
BiocManager::install(c("GENIE3"))

# Highly recommended packages 
BiocManager::install(c("zoo", "mixtools", "rbokeh"))
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
BiocManager::install(c("doMC", "doRNG"))

# To export/visualise 
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)


# SCENIC (This is the main package for this analysis)
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC") 
packageVersion("SCENIC")