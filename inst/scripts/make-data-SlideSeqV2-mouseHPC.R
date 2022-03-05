##################################################################
# Script to create Slide-seqV2 mouse hippocampus (HPC) data object
# Lukas Weber, updated Mar 2022
##################################################################

# references:
# Stickels et al. (2020): https://www.nature.com/articles/s41587-020-0739-1
# Cable et al. (2021): https://www.nature.com/articles/s41587-021-00830-w

# note: requires older version of RCTD package (before renaming)
# remotes::install_github("dmcable/spacexr", ref="ad6420ce4085daed520173d5c90ce055892b019f")


library(SpatialExperiment)
library(readr)
library(RCTD)


# ---------
# Load data
# ---------

# load data object from Stickels et al. (2020)
# previously downloaded and saved locally

dir_data <- "~/Dropbox/data/STexampleData/Slide_seqV2_mouse_hippo"


# load expression matrix
# note runtime: several minutes
exprs <- read_table(file.path(dir_data, "Puck_200115_08.digital_expression.txt.gz"))

dim(exprs)
format(object.size(exprs), units = "GB")
exprs[1:6, 1:6]

# get gene IDs from first column
gene_ids <- exprs$GENE
str(gene_ids)
length(gene_ids)

# get bead IDs from column names (excluding gene IDs column)
bead_ids <- colnames(exprs)[-1]
str(bead_ids)
length(bead_ids)

# convert expression matrix to numeric matrix without gene IDs
exprs <- exprs[, -1]
exprs <- as.matrix(exprs)
stopifnot(nrow(exprs) == length(gene_ids))
rownames(exprs) <- gene_ids

dim(exprs)
format(object.size(exprs), units = "GB")
exprs[1:6, 1:6]


# load beads info
# note: contains mix of tab-delimited and comma-delimited
file <- file.path(dir_data, "Puck_200115_08_bead_locations.csv")
bead_locations_colnames <- unlist(strsplit(readLines(file, n = 1), "\t"))
bead_locations <- read_csv(file, skip = 1, col_names = FALSE)
colnames(bead_locations) <- bead_locations_colnames
dim(bead_locations)
head(bead_locations)

stopifnot(nrow(bead_locations) == ncol(exprs))
stopifnot(all(bead_ids == bead_locations$barcodes))


# ----------------
# Cell type labels
# ----------------

# load predicted cell type labels from Cable et al. (2021), Figure 5A
# from data object provided by email by authors of Cable et al. (2021)

# load file containing RCTD outputs
out_rctd <- load(file.path(dir_data, "SlideseqHippo.RData"))

# extract top predicted cell type labels
labels_rctd <- results_df$first_type

# extract barcode IDs
barcodes_rctd <- puck@counts@Dimnames[[2]]
names(labels_rctd) <- barcodes_rctd

head(labels_rctd)
str(labels_rctd)


# ------------------------
# Create SpatialExperiment
# ------------------------

# convert assay matrix to sparse format
exprs_sparse <- as(exprs, "dgCMatrix")
format(object.size(exprs_sparse), units = "GB")

# row data
row_data <- DataFrame(gene_name = gene_ids)
rownames(row_data) <- row_data$gene_name

# column data
col_data <- DataFrame(
  barcode_id = bead_ids, 
  sample_id = "sample01"
)
rownames(col_data) <- col_data$barcode_id

col_data$celltype <- as.character(NA)
col_data[names(labels_rctd), "celltype"] <- as.character(labels_rctd)

# spatial coordinates
spatial_coords <- as.matrix(bead_locations[, c("xcoord", "ycoord")])
rownames(spatial_coords) <- bead_ids


spe <- SpatialExperiment(
  assays = list(counts = exprs_sparse), 
  rowData = row_data, 
  colData = col_data, 
  spatialCoords = spatial_coords
)

spe


# ----------------
# Save data object
# ----------------

saveRDS(spe, file = "SlideSeqV2_mouseHPC.rds")

