################################################################################
# Script to create SpatialExperiment object of the STARmap PLUS mouse brain
# well 05. 
# Yixing Dong, updated Mar 2024
################################################################################

# references:
#' STARmap PLUS `raw_expression_pd.csv` and `spatial.csv` file were downloaded from
#' \url{https://zenodo.org/records/8327576}

# in this script we download the STARmap PLUS data and reshape it into a 
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
star_well05_path <- here::here("raw_data/star_well05")
star_folder <- list.files(star_well05_path, pattern = ".csv")
star_folder


# ------------------------
# Create SpatialExperiment
# ------------------------

star_spe <- readStarmapplusSXE(star_well05_path, return_type = "SPE")
star_spe

# ------------------
# Saving data object
# ------------------

# We save the annotated Chromium to file for upload to `r Biocpkg("ExperimentHub")`.
repath <- here::here(file.path("raw_rds", "Starmapplus_mouseBrain"))
dir.create(repath, showWarnings=FALSE, recursive=TRUE)
saveRDS(star_spe, file=file.path(here::here(repath, "Starmapplus_mouseBrain.rds")))

