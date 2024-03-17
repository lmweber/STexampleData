################################################################################
# Script to create SingleCellExperiment object of the Chromium human breast cancer
# in 10x Xenium paper by Janesick et al. (2023)
# Yixing Dong, updated Mar 2024
################################################################################

# references:
#' Chromium `Feature / cell matrix HDF5 (per-sample)` .h5 file was downloaded from
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7782698}
#' and the computational cell type annotation was downloaded from 
#' \url{https://cdn.10xgenomics.com/raw/upload/v1695234604/Xenium%20Preview%20Data/Cell_Barcode_Type_Matrices.xlsx}

# in this script we download the chromium data and reshape it into a 
# SingleCellExperiment object

library(BiocFileCache)
library(DropletUtils)
library(SingleCellExperiment)


# -------------
# Download data
# -------------

# In the 10X Xenium paper, the Chromium data was accompanied by consecutive 
# slices of Visium and Xenium replicates. 
# We download and cache the count matrix and manual cell type annotation metadata 
# derived by 10x Genomics using the `r Biocpkg("BiocFileCache")` package.

bfc <- BiocFileCache("raw_data/chromium", ask = FALSE)
countmat <- bfcadd(bfc, "raw_data/chromium/archive/chrom_count.h5", fpath = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7782698&format=file&file=GSM7782698%5Fcount%5Fraw%5Ffeature%5Fbc%5Fmatrix%2Eh5")
anno <- bfcrpath(bfc, "https://zenodo.org/records/10076046/files/Cell_Barcode_Type_Matrices.xlsx?download=1")


# ---------------
# Processing data
# ---------------

# We load the counts into memory.
sce <- DropletUtils::read10xCounts(countmat, type = "HDF5")
colnames(sce) <- sce$Barcode
rownames(sce) <- rowData(sce)$Symbol


# Processing the per-cell annotation, and store as a column `Annotation` in the 
# `colData()` of the SCE object.
# # Note that the ensemble ID are stored in `rowData()` and the `rownames()` of 
# the object are gene symbol.
library(readxl)
celltype <- readxl::read_excel(anno, sheet = "scFFPE-Seq")

cells_overlap <- intersect(colnames(sce), celltype$Barcode)
sce <- sce[, cells_overlap]
celltype <- celltype[celltype$Barcode %in% cells_overlap, ]

colData(sce) <- as(merge(data.frame(colData(sce)), celltype), "DFrame")
colnames(sce) <- sce$Barcode
sce


# ------------------
# Saving data object
# ------------------

# We save the annotated Chromium to file for upload to `r Biocpkg("ExperimentHub")`.
repath <- here::here(file.path("raw_rds", "Chromium_10xJanesick2022_humanBreast"))
dir.create(repath, showWarnings=FALSE, recursive=TRUE)
saveRDS(sce, file=file.path(here::here(repath, "Chromium_10xJanesick2022_humanBreast.rds")))

