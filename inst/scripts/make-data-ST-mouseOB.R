################################################################################
# Script to create Spatial Transcriptomics (ST) mouse olfactory bulb (OB) object
# Lukas Weber, updated Jan 2022
################################################################################

# references:
# original data:
# Stahl et al. (2016): https://www.science.org/doi/full/10.1126/science.aaf2403
# object provided in MERINGUE package:
# Miller et al. (2021): https://genome.cshlp.org/content/early/2021/05/25/gr.271288.120

# in this script we download the object from the MERINGUE package and reshape it
# into a SpatialExperiment object


library(SpatialExperiment)
library(MERINGUE)


# ---------
# Load data
# ---------

# load data object from MERINGUE package

data(mOB)

str(mOB)
names(mOB)


# ------------------
# Extract components
# ------------------

# counts
counts <- mOB$counts

# gene names
row_data <- DataFrame(gene_name = rownames(mOB$counts))
rownames(row_data) <- row_data$gene_name

# barcode IDs
col_data <- DataFrame(
  barcode_id = colnames(mOB$counts), 
  sample_id = "sample01"
)
rownames(col_data) <- col_data$barcode_id

# manually annotated ground truth layers
# note: missing for some spots, so need to match correctly
layers <- rep(NA, nrow(col_data))
names(layers) <- rownames(col_data)
layers[names(mOB$annot)] <- as.character(mOB$annot)
# simplify names
layers <- gsub(" ", "_", gsub(" layer$", "", gsub("^[0-9]: ", "", tolower(layers))))

table(layers, useNA = "always")
col_data$layer <- layers

# spatial coordinates
spatial_coords <- as.matrix(mOB$pos)
# note: re-order to match col_data
stopifnot(nrow(col_data) == nrow(spatial_coords))
spatial_coords <- spatial_coords[rownames(col_data), ]
stopifnot(all(rownames(col_data) == rownames(spatial_coords)))


# ------------------------
# Create SpatialExperiment
# ------------------------

spe <- SpatialExperiment(
  assays = list(counts = counts), 
  rowData = row_data, 
  colData = col_data, 
  spatialCoords = spatial_coords
)

spe


# ----------------
# Save data object
# ----------------

saveRDS(spe, file = "ST_mouseOB.rds")

