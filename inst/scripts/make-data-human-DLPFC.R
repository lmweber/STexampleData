########################################################
# Script to create human DLPFC data object from raw data
########################################################

# for more details on this dataset see:
# http://spatial.libd.org/spatialLIBD/

# data can be downloaded as a SingleCellExperiment object through spatialLIBD
# package from Bioconductor (http://bioconductor.org/packages/spatialLIBD) -
# this contains all 12 samples, so here we subset sample 151673

# raw data files are also available from:
# http://spatial.libd.org/spatialLIBD/
# https://github.com/LieberInstitute/HumanPilot/tree/master/10X/151673


library(ExperimentHub)
library(spatialLIBD)
library(rjson)
library(SingleCellExperiment)
library(SpatialExperiment)


# ---------
# Load data
# ---------

# load data for all 12 samples from spatialLIBD package

ehub <- ExperimentHub()

# download and load full dataset (12 samples)
sce <- fetch_data(type = "sce", eh = ehub)

# subset to keep sample 151673 only
sce_sub <- sce[, sce$sample_name == "151673"]


# ----------------
# Load image files
# ----------------

# download image files separately from links at:
# http://spatial.libd.org/spatialLIBD/
# and:
# https://github.com/LieberInstitute/HumanPilot/tree/master/10X/151673

# saved locally in the following location
dir_local <- "~/data/HumanPilot/10X/151673"

# image file paths
# note: not including "aligned_fiducials.jpg" or "detected_tissue_image.jpg" from Loupe
img_paths <- c(
  tissue_hires_image = file.path(dir_local, "tissue_hires_image.png"), 
  tissue_lowres_image = file.path(dir_local, "tissue_lowres_image.png")
)

# spatial scale factors
file_scale_factors <- file.path(dir_local, "scalefactors_json.json")
scale_factors <- fromJSON(file = file_scale_factors)


# ------------------------
# Create SpatialExperiment
# ------------------------

# reformat SingleCellExperiment into SpatialExperiment, keeping only minimal
# columns from rowData and colData

# counts
counts <- assays(sce_sub)[["counts"]]

# row data
row_data <- rowData(sce_sub)[, c("gene_id", "gene_name", "type", "gene_source", 
                                 "gene_version", "gene_biotype")]
rownames(row_data) <- rowData(sce_sub)$gene_id

# column data
col_data <- colData(sce_sub)[, c("barcode", "row", "col", "imagerow", "imagecol", 
                                 "height", "width", "cell_count")]
col_data$ground_truth <- colData(sce_sub)[, "layer_guess_reordered"]
# add custom sample ID
# note: currently not working with custom sample ID
#col_data$sample_id <- "sample_01"
colnames(col_data)[1] <- "barcode_id"
rownames(col_data) <- colData(sce_sub)$barcode

# spatial coordinates
# add custom "x_coord" and "y_coord" with flipped/reversed coordinates for Visium platform
spatial_coords <- colData(sce_sub)[, c("barcode", "tissue")]
colnames(spatial_coords) <- c("barcode_id", "in_tissue")
spatial_coords$x_coord <- colData(sce_sub)[, "imagecol"]
spatial_coords$y_coord <- -1 * colData(sce_sub)[, "imagerow"] + max(colData(sce_sub)[, "imagerow"]) + 1
# note: column "in_tissue" must be logical
spatial_coords$in_tissue <- as.logical(as.numeric(spatial_coords$in_tissue))
rownames(spatial_coords) <- colData(sce_sub)$barcode

# image data
# both low and high resolution images from Space Ranger
img_data <- readImgData(
  path = dir_local, 
  sample_id = "Sample01", 
  imageSources = c(img_paths["tissue_lowres_image"], img_paths["tissue_hires_image"]), 
  scaleFactors = file_scale_factors, 
  load = TRUE
)

# create SpatialExperiment
spe <- SpatialExperiment(
  assays = list(counts = counts), 
  rowData = row_data, 
  colData = col_data, 
  spatialCoords = spatial_coords, 
  imgData = img_data
)

spe


# ----------------
# Save data object
# ----------------

# to do: move to ExperimentHub
# for now: saving as publicly accessible Dropbox link
save(spe, file = "~/Dropbox/STdata/human_DLPFC.RData")

