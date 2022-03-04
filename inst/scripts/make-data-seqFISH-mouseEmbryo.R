####################################################
# Script to create seqFISH data object from raw data
# Lukas Weber, updated Mar 2022
####################################################

# link to paper (Lohoff and Ghazanfar et al. 2020):
# https://www.biorxiv.org/content/10.1101/2020.11.20.391896v1

# links to raw data:
# https://marionilab.cruk.cam.ac.uk/SpatialMouseAtlas/
# https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/
# https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/README.md


library(SpatialExperiment)
library(BumpyMatrix)
library(tidyr)


# ---------
# Load data
# ---------

# data files downloaded and saved locally from links at https://marionilab.cruk.cam.ac.uk/SpatialMouseAtlas/

dir_data <- "~/data/STexampleData/SpatialMouseAtlas2020"

# contains cell-level info
metadata <- readRDS(file.path(dir_data, "metadata.Rds"))
# contains expanded lists of values for cell segmentation vertices (this information is also in metadata)
segmentation_vertices <- readRDS(file.path(dir_data, "segmentation_vertices.Rds"))
# contains info (x-y coordinates, intensity) for each mRNA molecule in each cell
mRNA <- readRDS(file.path(dir_data, "mRNA.Rds"))
# contains counts per gene and cell
counts <- readRDS(file.path(dir_data, "counts.Rds"))

head(metadata)
dim(metadata)

head(segmentation_vertices)
dim(segmentation_vertices)

head(mRNA)
dim(mRNA)

counts[1:6, 1:6]
dim(counts)


# -----------------
# Subset one sample
# -----------------

# subset one sample for the example dataset: embryo1, z == 2

metadata_sub <- metadata[metadata$embryo == "embryo1" & metadata$z == 2, ]
dim(metadata_sub)
head(metadata_sub)

mRNA_cols <- data.frame(x = mRNA$uniqueID) %>% separate(x, c("embryo", "pos", "cell", "z"))
stopifnot(nrow(mRNA_cols) == nrow(mRNA))
mRNA_sub <- mRNA[mRNA_cols$embryo == "embryo1" & mRNA_cols$z == "z2", ]
dim(mRNA_sub)
head(mRNA_sub)

counts_sub <- counts[, rownames(metadata_sub)]
dim(counts_sub)
counts_sub[1:6, 1:6]

stopifnot(ncol(counts_sub) == nrow(metadata_sub))
stopifnot(all(colnames(counts_sub) == rownames(metadata_sub)))


# ------------------------
# Create SpatialExperiment
# ------------------------

# assay: counts
dim(counts_sub)
counts_sub[1:6, 1:6]

# assay: BumpyMatrix containing mRNA molecule positions and other information
bumpy_cols <- c("x_global_affine", "y_global_affine", "seeds", "intensity")
bumpy_assay <- splitAsBumpyMatrix(mRNA_sub[, bumpy_cols], row = mRNA_sub$geneID, column = mRNA_sub$uniqueID)
# remove cells/columns missing from counts matrix and match column order
bumpy_assay <- bumpy_assay[, colnames(counts_sub)]

dim(bumpy_assay)
stopifnot(ncol(counts_sub) == ncol(bumpy_assay))
stopifnot(all(colnames(counts_sub) == colnames(bumpy_assay)))

# row data
row_data <- DataFrame(gene_name = rownames(counts_sub))
rownames(row_data) <- rownames(counts_sub)

# column data
stopifnot(all(metadata_sub$uniqueID == colnames(counts_sub)))
col_data <- DataFrame(metadata_sub)

# store segmentation vertices in colData as SplitDataFrameList (list of data frames)
seg_verts <- SplitDataFrameList(
  x = col_data$segmentation_vertices_x_global_affine, 
  y = col_data$segmentation_vertices_y_global_affine, 
  cbindArgs = TRUE
)
# remove previous format and store SplitDataFrameList in colData
col_data <- col_data[, 1:12]
col_data <- cbind(col_data, segmentation_vertices = I(seg_verts))

# add sample IDs to colData
sample_ids <- paste(col_data$embryo, paste0("z", col_data$z), sep = "_")
stopifnot(length(unique(sample_ids)) == 1)
col_data$sample_id <- sample_ids

# rename column of cell IDs
colnames(col_data)[1] <- "cell_id"

# spatial coordinates
spatial_coords <- as.matrix(metadata_sub[, c("x_global_affine", "y_global_affine")])

# create SpatialExperiment
spe <- SpatialExperiment(
  assays = list(
    counts = counts_sub, 
    molecules = bumpy_assay), 
  rowData = row_data, 
  colData = col_data, 
  spatialCoords = spatial_coords
)

spe


# ----------------
# Save data object
# ----------------

saveRDS(spe, file = "seqFISH_mouseEmbryo.rds")

