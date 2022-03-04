###############################################################
# Script to create Visium human DLPFC data object from raw data
# Lukas Weber, updated Mar 2022
###############################################################

# For more details on this dataset see:
# Maynard and Collado-Torres et al. (2020): https://www.nature.com/articles/s41593-020-00787-0
# spatialLIBD website: http://spatial.libd.org/spatialLIBD/

# The full dataset (12 samples) can also be downloaded as a SingleCellExperiment
# object through the spatialLIBD package from Bioconductor at:
# http://bioconductor.org/packages/spatialLIBD

# Raw data files are also available from:
# http://spatial.libd.org/spatialLIBD/
# https://github.com/LieberInstitute/HumanPilot/tree/master/10X/151673

# Here we build a SpatialExperiment object for 1 sample only (sample 151673),
# from raw data files (http://spatial.libd.org/spatialLIBD/) together with
# ground truth labels from the spatialLIBD Bioconductor package
# (http://bioconductor.org/packages/spatialLIBD)


library(SpatialExperiment)
library(SingleCellExperiment)
library(Matrix)
library(rjson)
library(dplyr)


# ---------
# Load data
# ---------

# data files downloaded and saved locally from links above (and/or JHPCE cluster)

dir_data <- "~/data/STexampleData/spatialLIBD"

# barcodes
file_barcodes <- file.path(dir_data, "JHPCE", "151673_raw_feature_bc_matrix__barcodes.tsv.gz")
df_barcodes <- read.csv(file_barcodes, sep = "\t", header = FALSE, 
                        col.names = c("barcode_id"))

# features
file_features <- file.path(dir_data, "JHPCE", "151673_raw_feature_bc_matrix__features.tsv.gz")
df_features <- read.csv(file_features, sep = "\t", header = FALSE, 
                        col.names = c("gene_id", "gene_name", "feature_type"))

# counts
file_counts <- file.path(dir_data, "JHPCE", "151673_raw_feature_bc_matrix__matrix.mtx.gz")
counts <- readMM(file = file_counts)

stopifnot(nrow(counts) == nrow(df_features))
stopifnot(ncol(counts) == nrow(df_barcodes))

# spatial coordinates
file_tisspos <- file.path(dir_data, "JHPCE", "tissue_positions_list.txt")
df_tisspos <- read.csv(file_tisspos, header = FALSE, 
                       col.names=c("barcode_id", "in_tissue", "array_row", "array_col", 
                                   "pxl_row_in_fullres", "pxl_col_in_fullres"))

# check dimensions
dim(df_barcodes)
dim(df_features)
dim(counts)
# note: df_tisspos contains all spots (even if using filtered data) and order of
# spots does not match df_barcodes; match and re-order below
dim(df_tisspos)

# image file paths
img_paths <- c(
  tissue_hires_image = file.path(dir_data, "JHPCE", "tissue_hires_image.png"), 
  tissue_lowres_image = file.path(dir_data, "JHPCE", "tissue_lowres_image.png")
)

# spatial scale factors
file_scale_factors <- file.path(dir_data, "JHPCE", "scalefactors_json.json")
scale_factors <- fromJSON(file = file_scale_factors)


# -------------------------------------------------
# Ground truth layer labels (and other column data)
# -------------------------------------------------

# SingleCellExperiment object from spatialLIBD Bioconductor package contains 
# ground truth layer labels and other useful barcode-level data in 'colData'

# here we extract these columns and match rows to 'df_barcodes'

# note: large download
#library(spatialLIBD)
#ehub <- ExperimentHub()
#sce <- fetch_data(type = "sce", eh = ehub)

# subset to keep sample 151673 only
#sce_sub <- sce[, sce$sample_name == "151673"]

# load from previously saved version
load(file.path(dir_data, "sce_sub.RData"))

# select columns to keep
cols_keep <- c("barcode", "cell_count", "layer_guess_reordered")
df_truth <- colData(sce_sub)[, cols_keep]
# rename columns
colnames(df_truth)[colnames(df_truth) == "barcode"] <- "barcode_id"
colnames(df_truth)[colnames(df_truth) == "layer_guess_reordered"] <- "ground_truth"
# convert labels to character
df_truth$ground_truth <- as.character(df_truth$ground_truth)
# add custom sample ID
df_truth$sample_id <- "sample_151673"

# note: contains only spots over tissue
dim(df_truth)


# --------------
# Match barcodes
# --------------

# match barcodes in df_barcodes and df_tisspos

# check order of rows (barcodes) in df_barcodes and df_tisspos
dim(df_barcodes)
dim(df_tisspos)
head(df_barcodes)
head(df_tisspos)
nrow(df_barcodes) == nrow(df_tisspos)
# order does not match
all(df_barcodes$barcode_id == df_tisspos$barcode_id)

# match and re-order rows in df_tisspos
ord <- match(df_barcodes$barcode_id, df_tisspos$barcode_id)
df_tisspos_ord <- df_tisspos[ord, ]
rownames(df_tisspos_ord) <- NULL

# check
dim(df_tisspos_ord)
stopifnot(nrow(df_barcodes) == nrow(df_tisspos_ord))
stopifnot(all(df_barcodes$barcode_id == df_tisspos_ord$barcode_id))
head(df_barcodes)
head(df_tisspos_ord)


# match barcodes in df_barcodes and df_truth

# check order of rows (barcodes) in df_barcodes and df_truth
dim(df_barcodes)
dim(df_truth)
head(df_barcodes)
head(df_truth)
# df_truth contains only spots over tissue, so number of rows does not match
nrow(df_barcodes) == nrow(df_truth)

# match rows
df_truth_matched <- left_join(df_barcodes, as.data.frame(df_truth), by = "barcode_id")
dim(df_truth_matched)
head(df_truth_matched)
stopifnot(nrow(df_barcodes) == nrow(df_truth_matched))
stopifnot(all(df_barcodes$barcode_id == df_truth_matched$barcode_id))

# replace sample ID column to include all spots (including not over tissue)
df_truth_matched$sample_id <- "sample_151673"


# ------------------------
# Create SpatialExperiment
# ------------------------

# row data
row_data <- DataFrame(df_features)
rownames(row_data) <- df_features$gene_id

# column data
stopifnot(nrow(df_truth_matched) == nrow(df_tisspos_ord))
stopifnot(all(df_truth_matched$barcode_id == df_tisspos_ord$barcode_id))
col_data <- DataFrame(cbind(df_truth_matched, df_tisspos_ord[, -1]))
col_data <- col_data[, c(1, 4, 5:7, 3, 2)]
rownames(col_data) <- col_data$barcode_id

# spatial coordinates
spatial_coords <- as.matrix(df_tisspos_ord[, c("pxl_col_in_fullres", "pxl_row_in_fullres")])
rownames(spatial_coords) <- df_tisspos_ord$barcode_id

# image data
# low and high resolution images from Space Ranger
img_data <- readImgData(
  path = file.path(dir_data, "JHPCE"), 
  sample_id = "sample_151673", 
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

saveRDS(spe, file = "Visium_humanDLPFC.rds")

