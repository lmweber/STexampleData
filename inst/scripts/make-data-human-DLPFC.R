########################################################
# Script to create human DLPFC data object from raw data
########################################################

# ----------------------
# Download and load data
# ----------------------

# load data from spatialLIBD Bioconductor package

# note: full dataset contains 12 samples; here we use sample 151673 only

# to do: replace to load from raw data links instead (see
# https://github.com/LieberInstitute/HumanPilot)

# to do: use SpatialExperiment instead of SingleCellExperiment

library(ExperimentHub)
library(spatialLIBD)

ehub <- ExperimentHub()

# download and load full dataset (12 samples)
sce <- fetch_data(type = "sce", eh = ehub)

# subset to keep sample 151673 only
sce_sub <- sce[, sce$sample_name == "151673"]


# ------------------
# Create data object
# ------------------

# create SingleCellExperiment object containing raw counts only
spe <- SingleCellExperiment(
  assays = list(counts = counts(sce_sub)), 
  rowData = rowData(sce_sub)[, 2:7], 
  colData = colData(sce_sub)[, c(1, 4:7, 19)]
)

# add column of ground truth layer labels
colData(spe)$ground_truth <- colData(sce_sub)$layer_guess_reordered

# convert x and y coordinates and store in colData
colData(spe)$x_coord <- colData(spe)[, "imagecol"]
colData(spe)$y_coord <- -colData(spe)[, "imagerow"]


# ----------------
# Save data object
# ----------------

save(spe, file = "~/Dropbox/STdata/human_DLPFC.RData")

