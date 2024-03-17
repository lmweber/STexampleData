################################################################################
# Script to create SpatialExperiment object of the Xenium 1 human breast cancer
# in 10x Xenium paper by Janesick et al. (2023)
# Yixing Dong, updated Mar 2024
################################################################################

# references:
#' XeniumOutputBundle .zip file was downloaded from 
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7780153}

# in this script we download the visium data and reshape it into a 
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

# In the 10X Xenium paper, the Xenium sample 1 (2 replicates) data was accompanied 
# by consecutive slices of Chromium and Visium replicates. 

# We first prepare the data folder required by reader function `SpatialExperimentIO::readXeniumSXE()`. 
# Please click to download [Xenium Sample 1 Replicate 1 outs](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7780153&format=file&file=GSM7780153%5FXenium%5FFFPE%5FHuman%5FBreast%5FCancer%5FRep1%5Fouts%2Ezip). 
# Here we placed the downloaded objects in `here::here("/raw_data/xenium_sample1rep1/")`. 

xenium_out_folder <- list.files(here::here("raw_data/xenium_sample1rep1"), pattern = ".zip")
xenium_out_folder


# We now unzip `spatial.tar.gz` to its same directory. 

unzip(file.path(here::here("raw_data/xenium_sample1rep1"), xenium_out_folder),
      exdir = here::here("raw_data/xenium_sample1rep1"))


# Have a look at all the files in the newly unzipped `/outs` folder. The only 
# files you need are `cell_feature_matrix.h5` or `cell_feature_matrix`, and 
# `cells.csv.gz`. 

list.files(here::here("raw_data/xenium_sample1rep1/outs"))

# Sanity check that the files you need are in the `/outs` folder.

all(c("cell_feature_matrix.h5", "cell_feature_matrix", "cells.csv.gz") %in% 
      list.files(here::here("raw_data/xenium_sample1rep1/outs")))


# ------------------------
# Create SpatialExperiment
# ------------------------

xe1_spe <- SpatialExperimentIO::readXeniumSXE(here::here("raw_data/xenium_sample1rep1/outs"), 
                                              return_type = "SPE")
xe1_spe <- xe1_spe[rowData(xe1_spe)$Type == "Gene Expression"]
rownames(xe1) <- rowData(xe1)$Symbol

# ------------------
# Saving data object
# ------------------

# We save the Xenium data to file for upload to `r Biocpkg("ExperimentHub")`.

repath <- here::here(file.path("raw_rds", "Xenium1_10xJanesick2022_humanBreast"))
dir.create(repath, showWarnings=FALSE, recursive=TRUE)
saveRDS(xe1_spe, file=file.path(repath, "Xenium1_10xJanesick2022_humanBreast.rds"))


