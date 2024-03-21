################################################################################
# Script to create SpatialExperiment object of the Merscope human ovarian cancer
# patient 1 slice 2. 
# Yixing Dong, updated Mar 2024
################################################################################

# references:
#' Merscope `_cell_by_gene.csv` and `_cell_metadata.csv` file were downloaded from
#' \url{https://console.cloud.google.com/storage/browser/vz-ffpe-showcase/HumanOvarianCancerPatient2Slice1;tab=objects?prefix=&forceOnObjectsSortingFiltering=false}

# in this script we download the merscope data and reshape it into a 
# SpatialExperiment object

# First, we need to install the `SpatialExperimentIO` loader package that would 
# return Xenium as a `SpatialExperiment` object: 

# install.packages("devtools")
devtools::install_github("estellad/SpatialExperimentIO")

# or with the version of the pacakge on Bioconductor: 

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SpatialExperimentIO", version = "devel")

library(SpatialExperiment)
library(SpatialExperimentIO)

# -------------
# Download data
# -------------
# Put the downloaded unzipped file into a folder. Make sure that two mandatory 
# files do exist. 
mer_lung_p2s1_path <- here::here("raw_data/mer_patient1slice2")
mer_folder <- list.files(mer_lung_p2s1_path, pattern = ".csv")
mer_folder


# ------------------------
# Create SpatialExperiment
# ------------------------

mer_spe <- readMerscopeSXE(mer_lung_p2s1_path, return_type = "SPE")
mer_spe

# ------------------
# Saving data object
# ------------------

# We save the annotated Chromium to file for upload to `r Biocpkg("ExperimentHub")`.
repath <- here::here(file.path("raw_rds", "Merscope_humanOvarian"))
dir.create(repath, showWarnings=FALSE, recursive=TRUE)
saveRDS(mer_spe, file=file.path(here::here(repath, "Merscope_humanOvarian.rds")))

