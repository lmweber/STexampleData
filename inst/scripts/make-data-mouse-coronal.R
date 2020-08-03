##########################################################
# Script to create mouse coronal data object from raw data
##########################################################

# -------------
# Download data
# -------------

# for more details on raw data see:
# https://support.10xgenomics.com/spatial-gene-expression/datasets
# https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain

url <- "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.tar.gz"
fn <- basename(url)
dir.create("tmp")
download.file(url, file.path("tmp", fn))

system(paste0("tar -C tmp -xvzf ", file.path("tmp", fn)))


# ---------
# Load data
# ---------

# download data matrix

# for more details on matrix format used by Cell Ranger see:
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices

library(Matrix)

matrix_dir <- file.path("tmp", "filtered_feature_bc_matrix")

barcode_path <- file.path(matrix_dir, "barcodes.tsv.gz")
features_path <- file.path(matrix_dir, "features.tsv.gz")
matrix_path <- file.path(matrix_dir, "matrix.mtx.gz")

mat <- readMM(file = matrix_path)

feature_names <- read.delim(features_path, header = FALSE, stringsAsFactors = FALSE)
barcode_names <- read.delim(barcode_path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) <- barcode_names$V1
rownames(mat) <- feature_names$V1


# download spot coordinates

url <- "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_spatial.tar.gz"
fn <- basename(url)
download.file(url, file.path("tmp", fn))

system(paste0("tar -C tmp -xvzf ", file.path("tmp", fn)))

spatial_dir <- file.path("tmp", "spatial")
tissue_positions_path <- file.path(spatial_dir, "tissue_positions_list.csv")
tissue_positions <- read.csv(tissue_positions_path, header = FALSE)

coords <- tissue_positions[, c("V1", "V5", "V6")]
colnames(coords) <- c("barcode", "x_original", "y_original")
rownames(coords) <- coords$barcode
dim(coords)

# flip x and y coordinates to match published images
coords$x_coord <- coords$y_original
coords$y_coord <- -coords$x_original


# ------------------
# Create data object
# ------------------

# set up SingleCellExperiment

# to do: use sparse matrix format

library(SingleCellExperiment)

# match x-y coordinates for spots containing data
coords_matched <- coords[colnames(mat), ]
dim(coords_matched)
stopifnot(all(rownames(coords_matched) == colnames(mat)))

# create SingleCellExperiment
col_data <- DataFrame(coords_matched)
row_data <- DataFrame(feature_names)
colnames(row_data) <- c("gene_id", "gene_name", "type")

spe <- SingleCellExperiment(
  assays = list(counts = as.matrix(mat)),  ## to do: use sparse matrix format
  rowData = row_data, 
  colData = col_data
)


# ----------------------
# Delete temporary files
# ----------------------

unlink("tmp", recursive = TRUE)


# ----------------
# Save data object
# ----------------

# to do: store on ExperimentHub

# for now: saving as publicly accessible Dropbox link

save(spe, file = "~/Dropbox/STdata/mouse_coronal.RData")

