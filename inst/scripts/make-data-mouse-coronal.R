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

dir.create("mouse_coronal")

# filtered feature-barcode matrix
url <- "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.tar.gz"
fn <- basename(url)
download.file(url, file.path("mouse_coronal", fn))
system(paste0("tar -C mouse_coronal -xvzf ", file.path("mouse_coronal", fn)))

# spatial imaging data
url <- "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_spatial.tar.gz"
fn <- basename(url)
download.file(url, file.path("mouse_coronal", fn))
system(paste0("tar -C mouse_coronal -xvzf ", file.path("mouse_coronal", fn)))


# ---------
# Load data
# ---------

# load features and barcodes

file_features <- file.path("mouse_coronal", "filtered_feature_bc_matrix", "features.tsv.gz")
file_barcodes <- file.path("mouse_coronal", "filtered_feature_bc_matrix", "barcodes.tsv.gz")

df_features <- read.delim(file_features, header = FALSE, stringsAsFactors = FALSE)
df_barcodes <- read.delim(file_barcodes, header = FALSE, stringsAsFactors = FALSE)

colnames(df_features) <- c("gene_id", "gene_name", "feature_type")
colnames(df_barcodes) <- c("barcode_id")

rownames(df_features) <- df_features$gene_id
rownames(df_barcodes) <- df_barcodes$barcode_id

head(df_features)
head(df_barcodes)


# load feature-barcode matrix

file_matrix <- file.path("mouse_coronal", "filtered_feature_bc_matrix", "matrix.mtx.gz")

mat <- readMM(file = file_matrix)

stopifnot(nrow(mat) == nrow(df_features))
stopifnot(ncol(mat) == nrow(df_barcodes))

rownames(mat) <- df_features$gene_id
colnames(mat) <- df_barcodes$barcode_id


# load spatial coordinates

file_tissue_positions <- file.path("mouse_coronal", "spatial", "tissue_positions_list.csv")

df_tissue_positions <- read.csv(file_tissue_positions, header = FALSE)

# note array row corresponds to pixel column and vice versa
colnames(df_tissue_positions) <- c("barcode_id", "in_tissue", 
                                   "array_row", "array_col", 
                                   "pxl_col_in_fullres", "pxl_row_in_fullres")

rownames(df_tissue_positions) <- df_tissue_positions$barcode_id

# keep spatial coordinates for barcodes that overlap with tissue
df_tissue_positions_keep <- df_tissue_positions[df_barcodes$barcode_id, ]

stopifnot(nrow(df_tissue_positions_keep) == nrow(df_barcodes))
table(df_tissue_positions_keep$in_tissue)


# load scale factors for image files

file_scale_factors <- file.path("mouse_coronal", "spatial", "scalefactors_json.json")

scale_factors <- fromJSON(file = file_scale_factors)


# paths to image files

image_paths <- c(
  aligned_fiducials = file.path("mouse_coronal", "spatial", "aligned_fiducials.jpg"), 
  detected_tissue_image = file.path("mouse_coronal", "spatial", "detected_tissue_image.jpg"), 
  tissue_hires_image = file.path("mouse_coronal", "spatial", "tissue_hires_image.png"), 
  tissue_lowres_image = file.path("mouse_coronal", "spatial", "tissue_lowres_image.png")
)


# -----------------------
# Create VisiumExperiment
# -----------------------

ve <- VisiumExperiment(
  rowData = df_features, 
  colData = df_barcodes, 
  assays = list(counts = mat), 
  spatialCoords = df_tissue_positions_keep, 
  scaleFactors = scale_factors, 
  imagePaths = image_paths
)


# ----------------------
# Delete temporary files
# ----------------------

unlink("mouse_coronal", recursive = TRUE)


# ----------------
# Save data object
# ----------------

# to do: move to ExperimentHub

# for now: saving as publicly accessible Dropbox link

save(ve, file = "~/Dropbox/STdata/mouse_coronal.RData")

