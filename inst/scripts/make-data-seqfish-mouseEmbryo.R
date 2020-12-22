####################################################
# Script to create seqFISH data object from raw data
# Lukas Weber, December 2020
####################################################

# link to paper (Lohoff and Ghazanfar et al. 2020):
# https://www.biorxiv.org/content/10.1101/2020.11.20.391896v1

# links to raw data:
# https://marionilab.cruk.cam.ac.uk/SpatialMouseAtlas/
# https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/
# https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/README.md


library(SpatialExperiment)
library(Matrix)
library(BumpyMatrix)
library(tidyverse)


# ---------
# Load data
# ---------

# data files downloaded and saved locally from links at https://marionilab.cruk.cam.ac.uk/SpatialMouseAtlas/

dir_data <- "~/data/SpatialMouseAtlas2020"

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

# row data
row_data <- data.frame(gene_name = rownames(counts_sub))
head(row_data)

# column data
stopifnot(all(metadata_sub$uniqueID == colnames(counts_sub)))
col_data <- metadata_sub
colnames(col_data)[1] <- "cell_id"
head(col_data)

# create SpatialExperiment
spe <- SpatialExperiment(
  assays = list(
    counts = counts_sub#, 
    #molecules = bumpy_assay
  ), 
  rowData = row_data, 
  colData = col_data
)

spe


# ----------------
# Save data object
# ----------------

save(spe, file = "~/Dropbox/STdata/seqFISH_mouseEmbryo.RData")

