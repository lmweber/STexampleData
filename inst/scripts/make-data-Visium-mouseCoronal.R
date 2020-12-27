#################################################################
# Script to create Visium mouse coronal data object from raw data
# Lukas Weber, December 2020
#################################################################

# for more details on raw data see:
# https://support.10xgenomics.com/spatial-gene-expression/datasets
# https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain

# for more details on matrix format used by Cell Ranger see:
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices


library(SpatialExperiment)
library(Matrix)
library(rjson)


# -------------
# Download data
# -------------

dir.create("tmp")

# raw feature-barcode matrix
# note: raw matrix contains all spots; filtered matrix contains only spots over tissue
url <- "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_raw_feature_bc_matrix.tar.gz"
fn <- basename(url)
download.file(url, file.path("tmp", fn))
system(paste0("tar -C tmp -xvzf ", file.path("tmp", fn)))

# spatial imaging data
url <- "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_spatial.tar.gz"
fn <- basename(url)
download.file(url, file.path("tmp", fn))
system(paste0("tar -C tmp -xvzf ", file.path("tmp", fn)))


# ---------
# Load data
# ---------

# barcodes
file_barcodes <- file.path("tmp", "raw_feature_bc_matrix", "barcodes.tsv.gz")
df_barcodes <- read.csv(file_barcodes, sep = "\t", header = FALSE, 
                        col.names = c("barcode_id"))

# features
file_features <- file.path("tmp", "raw_feature_bc_matrix", "features.tsv.gz")
df_features <- read.csv(file_features, sep = "\t", header = FALSE, 
                        col.names = c("gene_id", "gene_name", "feature_type"))

# counts
file_counts <- file.path("tmp", "raw_feature_bc_matrix", "matrix.mtx.gz")
counts <- readMM(file = file_counts)

stopifnot(nrow(counts) == nrow(df_features))
stopifnot(ncol(counts) == nrow(df_barcodes))

# spatial coordinates
file_tisspos <- file.path("tmp", "spatial", "tissue_positions_list.csv")
df_tisspos <- read.csv(file_tisspos, header = FALSE, 
                       col.names=c("barcode_id", "in_tissue", "array_row", "array_col", 
                                   "pxl_col_in_fullres", "pxl_row_in_fullres"))

# check dimensions
dim(df_barcodes)
dim(df_features)
dim(counts)
# note: df_tisspos contains all spots (even if using filtered data) and order of
# spots does not match df_barcodes; match and re-order below
dim(df_tisspos)

# image file paths
img_paths <- c(
  aligned_fiducials = file.path("tmp", "spatial", "aligned_fiducials.jpg"), 
  detected_tissue_image = file.path("tmp", "spatial", "detected_tissue_image.jpg"), 
  tissue_hires_image = file.path("tmp", "spatial", "tissue_hires_image.png"), 
  tissue_lowres_image = file.path("tmp", "spatial", "tissue_lowres_image.png")
)

# spatial scale factors
file_scale_factors <- file.path("tmp", "spatial", "scalefactors_json.json")
scale_factors <- fromJSON(file = file_scale_factors)


# --------------
# Match barcodes
# --------------

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


# ------------------------
# Create SpatialExperiment
# ------------------------

# row data
row_data <- df_features
rownames(row_data) <- df_features$gene_id

# column data
col_data <- df_barcodes
# add custom sample ID
# note: currently not working with custom sample ID
#col_data$sample_id <- "sample_01"
rownames(col_data) <- df_barcodes$barcode_id

# spatial data
spatial_data <- df_tisspos_ord[, c("barcode_id", "in_tissue")]
# flip x and y coordinates and reverse y scale to match orientation of images
# note: for this dataset, also need to add 'min(y_coord) + max(y_coord)'
spatial_data$x_coord <- df_tisspos_ord$pxl_row_in_fullres
y_coord_tmp <- df_tisspos_ord$pxl_col_in_fullres
y_coord_tmp <- (-1 * y_coord_tmp) + min(y_coord_tmp) + max(y_coord_tmp)
spatial_data$y_coord <- y_coord_tmp
# note: column 'in_tissue' must be logical
spatial_data$in_tissue <- as.logical(spatial_data$in_tissue)
rownames(spatial_data) <- df_tisspos_ord$barcode_id

# additional column data
# keep columns with raw coordinates (may be useful for some users)
col_data_additional <- df_tisspos_ord[, c("array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres")]
rownames(col_data_additional) <- df_tisspos_ord$barcode_id
stopifnot(nrow(col_data) == nrow(col_data_additional))
stopifnot(all(rownames(col_data) == rownames(col_data_additional)))
col_data <- cbind(col_data, col_data_additional)

# checks
stopifnot(nrow(col_data) == nrow(spatial_data))
stopifnot(all(rownames(col_data) == rownames(spatial_data)))

# image data
# both low and high resolution images from Space Ranger
img_data <- readImgData(
  path = file.path("tmp", "spatial"), 
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
  spatialCoords = spatial_data, 
  imgData = img_data
)

spe


# ----------------------
# Delete temporary files
# ----------------------

unlink("tmp", recursive = TRUE)


# ----------------
# Save data object
# ----------------

save(spe, file = "~/Dropbox/STdata/Visium_mouseCoronal.RData")

