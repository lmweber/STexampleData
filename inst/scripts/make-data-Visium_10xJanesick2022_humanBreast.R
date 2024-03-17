################################################################################
# Script to create SpatialExperiment object of the Visium human breast cancer
# in 10x Xenium paper by Janesick et al. (2023)
# Yixing Dong, updated Mar 2024
################################################################################

# references:
#' Visium `Feature / cell matrix HDF5 (per-sample)` .h5 file was downloaded from
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7782699}

# in this script we download the visium data and reshape it into a 
# SpatialExperiment object

library(SpatialExperiment)
library(dplyr)


# -------------
# Download data
# -------------

# In the 10X Xenium paper, the Visium data was accompanied by consecutive 
# slices of Chromium and Xenium replicates. 


# We first prepare the `outs/` data folder required by reader function 
# `SpatialExperiment::read10xVisium()`. Please click to download 
# [filtered_feature_bc_matrix.h5](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7782699&format=file&file=GSM7782699%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5) 
# and folder [spatial.tar.gz](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7782699&format=file&file=GSM7782699%5Fspatial%2Etar%2Egz), 
# and place them in a folder named `/outs`. 
# Here we placed the downloaded objects in `here::here("/raw_data/visium/outs/")`. 
list.files(here::here("raw_data/visium/outs/"))


# We need to unzip the spatial.tar.gz into the same directory. First check what 
# is in the folder `spatial.tar.gz`.
spatial_folder <- list.files(here::here("raw_data/visium/outs/"), pattern = ".tar.gz")
untar(file.path(here::here("raw_data/visium/outs"), spatial_folder), list=TRUE)

# We now unzip `spatial.tar.gz` to its same directory. 
untar(file.path(here::here("raw_data/visium/outs"), spatial_folder),
      exdir = here::here("raw_data/visium/outs"))


# We also need to rename the count matrix `GSM7782699_filtered_feature_bc_matrix.h5` 
# to only `filtered_feature_bc_matrix.h5`. 
vis_h5 <- list.files(here::here("raw_data/visium/outs/"), pattern = ".h5")
file.rename(from = file.path(here::here("raw_data/visium/outs/"), vis_h5), 
            to = file.path(here::here("raw_data/visium/outs/"), "filtered_feature_bc_matrix.h5"))

# Now check that you have both `/spatial` folder and `filtered_feature_bc_matrix.h5` 
# count matrix in the `/outs` directory. 
list.files(here::here("raw_data/visium/outs/"))


# ------------------------
# Create SpatialExperiment
# ------------------------

# We can finally read in the Visium object with the reader function below. 
vis <- SpatialExperiment::read10xVisium(here::here("raw_data/visium/outs"), 
                                        type = "HDF5",
                                        data = "filtered",
                                        image = "lowres")
vis


# Check that there are some duplicated gene names that we should rename to avoid 
# problems with downstream analysis. 
head(sort(table(rowData(vis)$symbol), decreasing = TRUE))


# We define the following rules of renaming for the three genes with duplicated 
# names, which is also consistent with the renaming rules of Seurat's Visium 
# reader function -  `Seurat::Load10X_Spatial()`. 

rowData(vis)$geneid <- rownames(vis)
RD <- data.frame(rowData(vis))
RD <- RD %>%
  mutate(symbol = case_when(
    geneid == "ENSG00000284770" ~ "TBCE.1",
    geneid == "ENSG00000187522" ~ "HSPA14.1",
    geneid == "ENSG00000269226" ~ "TMSB15B.1",
    .default = symbol)
  )

length(unique(RD$symbol)) == length(RD$symbol)


# We have now kept the ensemble ID in `rowData()` and set the `rownames()` of 
# the object to gene symbol.

rowData(vis) <- as(RD, "DFrame")
rownames(vis) <- rowData(vis)$symbol

vis

# ------------------
# Saving data object
# ------------------

# We save the Visium to file for upload to `r Biocpkg("ExperimentHub")`.

repath <- here::here(file.path("raw_rds", "Visium_10xJanesick2022_humanBreast"))
dir.create(repath, showWarnings=FALSE, recursive=TRUE)
saveRDS(vis, file=file.path(here::here(repath, "Visium_10xJanesick2022_humanBreast.rds")))

