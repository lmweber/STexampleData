# Script to create human DLPFC data object from raw data

# TO DO: update to start from raw data
# for now: using spatialLIBD package as source


## Load data

## Load the dataset from the
## [spatialLIBD](http://bioconductor.org/packages/spatialLIBD) Bioconductor
## package. Note that we use only one sample (sample 151673) for these examples.

## note this downloads the full dataset (12 samples) from spatialLIBD

library(ExperimentHub)
library(spatialLIBD)

ehub <- ExperimentHub()

# load full dataset
sce <- fetch_data(type = "sce", eh = ehub)
# subset sample 151673
sce <- sce[, sce$sample_name == "151673"]
# create object containing raw data only
sce <- SingleCellExperiment(
  rowData = rowData(sce)[, 1:8], 
  colData = colData(sce)[, c(1:7, 9:18)], 
  assays = list(counts = assays(sce)[["counts"]])
)


## Convert the `SingleCellExperiment` object to a `SpatialExperiment`.

#library(SpatialExperiment)

## TO DO: either convert from SCE to SPE, or construct SPE directly from raw
## Visium data; possibly using raw data links from
## https://github.com/LieberInstitute/HumanPilot

# for now: use SCE object instead
spe <- sce

# for now: keep x and y spatial coordinates in colData
colData(spe)$x_coord <- colData(spe)[, "imagecol"]
colData(spe)$y_coord <- -colData(spe)[, "imagerow"]

spe


