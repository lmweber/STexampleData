#################################################################
# Script to create Visium mouse coronal data object from raw data
# Lukas Weber, updated Jan 2022
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
row_data <- DataFrame(df_features)
rownames(row_data) <- df_features$gene_id

# column data
# include column of sample IDs
col_data <- DataFrame(df_tisspos_ord, sample_id = "sample01")
col_data <- col_data[, c(1, 7, 2:6)]
rownames(col_data) <- col_data$barcode_id

# spatial coordinates
# use default column names 'x' and 'y'
spatial_coords <- as.matrix(df_tisspos_ord[, c("pxl_row_in_fullres", "pxl_col_in_fullres")])
colnames(spatial_coords) <- c("x", "y")
rownames(spatial_coords) <- df_tisspos_ord$barcode_id

# image data
# low and high resolution images from Space Ranger
img_data <- readImgData(
  path = file.path("tmp", "spatial"), 
  sample_id = "sample01", 
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


# ----------------------
# Delete temporary files
# ----------------------

unlink("tmp", recursive = TRUE)


# ----------------
# Save data object
# ----------------

saveRDS(spe, file = "Visium_mouseCoronal.rds")

