##########################################################
# Script to create mouse coronal data object from raw data
##########################################################

# for more details on raw data see:
# https://support.10xgenomics.com/spatial-gene-expression/datasets
# https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain

# for more details on matrix format used by Cell Ranger see:
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices


library(Matrix)
library(rjson)
library(SpatialExperiment)


# -------------
# Download data
# -------------

dir.create("tmp")

# filtered feature-barcode matrix
# note: filtered matrix contains only spots over tissue; raw matrix contains all spots
url <- "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.tar.gz"
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
file_barcodes <- file.path("tmp", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
df_barcodes <- read.csv(file_barcodes, sep = "\t", header = FALSE, 
                        col.names = c("barcode_id"))

# features
file_features <- file.path("tmp", "filtered_feature_bc_matrix", "features.tsv.gz")
df_features <- read.csv(file_features, sep = "\t", header = FALSE, 
                        col.names = c("gene_id", "gene_name", "feature_type"))

# counts
file_counts <- file.path("tmp", "filtered_feature_bc_matrix", "matrix.mtx.gz")
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
# note: df_tisspos contains all spots (not filtered), so need to match spots below
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

# match and re-order barcode IDs (rows) in df_barcodes and df_tisspos
dim(df_barcodes)
dim(df_tisspos)
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

# column data
# note: add sample ID
col_data <- df_barcodes
# note: currently not working with custom sample ID
#col_data$sample_id <- "sample_01"

# spatial coordinates
# note: add "x_coord" and "y_coord" with flipped/reversed coordinates (standard orientation)
spatial_coords <- df_tisspos_ord
# note: including both x/y_coord and pxl_row/col_in_fullres currently not working
#spatial_coords$x_coord <- spatial_coords$pxl_row_in_fullres
#spatial_coords$y_coord <- -1 * spatial_coords$pxl_col_in_fullres + max(spatial_coords$pxl_col_in_fullres) + 1
# note: column "in_tissue" must be logical
spatial_coords$in_tissue <- as.logical(spatial_coords$in_tissue)

# image data
# note: reading in both low and high resolution image from Space Ranger
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

# to do: move to ExperimentHub

# for now: saving as publicly accessible Dropbox link

save(spe, file = "~/Dropbox/SpaData/mouse_coronal.RData")

