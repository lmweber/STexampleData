################################################################################
# Script to create SpatialExperiment object of the CosMx human lung cancer
# patient 9 slice 1. 
# Yixing Dong, updated Mar 2024
################################################################################

# references:
#' CosMx `_exprMat_file.csv` and `_metadata_file.csv` file were downloaded from
#' \url{https://nanostring.com/resources/smi-ffpe-dataset-lung9-rep1-data/}

# in this script we download the cosmx data and reshape it into a 
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
cosmx_lung_p9s1_path <- here::here("raw_data/cosmx_patient9slice1")
cosmx_folder <- list.files(cosmx_lung_p9s1_path, pattern = ".csv")
cosmx_folder


# ------------------------
# Create SpatialExperiment
# ------------------------

cos_spe <- readCosmxSXE(cosmx_lung_p9s1_path, return_type = "SPE")
cos_spe

# ------------------
# Saving data object
# ------------------

# We save the annotated Chromium to file for upload to `r Biocpkg("ExperimentHub")`.
repath <- here::here(file.path("raw_rds", "Cosmx_humanLung"))
dir.create(repath, showWarnings=FALSE, recursive=TRUE)
saveRDS(cos_spe, file=file.path(here::here(repath, "Cosmx_humanLung.rds")))

