# ---------------------------
# Create metadata spreadsheet
# ---------------------------

# metadata for all datasets

df_all <- data.frame(
  BiocVersion = "3.13", 
  Genome = NA, 
  SourceType = "FASTQ", 
  SourceVersion = NA, 
  Coordinate_1_based = NA, 
  DataProvider = NA, 
  Maintainer = "Lukas M. Weber <lukas.weber.edu@gmail.com>", 
  RDataClass = "SpatialExperiment", 
  DispatchClass = "Rds", 
  stringsAsFactors = FALSE
)


# metadata for individual datasets

df_Visium_humanDLPFC <- cbind(
  df_all, 
  Title = "Visium_humanDLPFC", 
  Description = paste0(
    "A single sample (sample 151673) of human brain dorsolateral prefrontal ", 
    "cortex (DLPFC) in the human brain, measured using the 10x Genomics Visium ", 
    "platform. This is a subset of the full dataset containing 12 samples from ", 
    "3 neurotypical donors, published by Maynard and Collado-Torres et al. ", 
    "(2021). The full dataset is available from the spatialLIBD Bioconductor ", 
    "package."), 
  SourceUrl = "http://spatial.libd.org/spatialLIBD/", 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataPath = "STexampleData/3_13/Visium_humanDLPFC.rds", 
  stringsAsFactors = FALSE
)

df_Visium_mouseCoronal <- cbind(
  df_all, 
  Title = "Visium_mouseCoronal", 
  Description = paste0(
    "A single coronal section from the mouse brain spanning one hemisphere, ", 
    "measured using the 10x Genomics Visium platform. This dataset was publicly ", 
    "released by 10x Genomics."), 
  SourceUrl = "https://support.10xgenomics.com/spatial-gene-expression/datasets", 
  Species = "Mus musculus", 
  TaxonomyId = "10090", 
  RDataPath = "STexampleData/3_13/Visium_mouseCoronal.rds", 
  stringsAsFactors = FALSE
)

df_seqFISH_mouseEmbryo <- cbind(
  df_all, 
  Title = "seqFISH_mouseEmbryo", 
  Description = paste0(
    "A subset of cells (embryo 1, z-slice 2) from a dataset investigating mouse ", 
    "embryogenesis by Lohoff and Ghazanfar et al. (2020), generated using the ", 
    "seqFISH platform. The full dataset is available from the original ", 
    "publication."), 
  SourceUrl = "https://marionilab.cruk.cam.ac.uk/SpatialMouseAtlas/", 
  Species = "Mus musculus", 
  TaxonomyId = "10090", 
  RDataPath = "STexampleData/3_13/seqFISH_mouseEmbryo.rds", 
  stringsAsFactors = FALSE
)


# combine and save as .csv spreadsheet file

df_combined <- rbind(
  df_Visium_humanDLPFC, 
  df_Visium_mouseCoronal, 
  df_seqFISH_mouseEmbryo
)

write.csv(df_combined, file = "../extdata/metadata.csv", row.names = FALSE)

