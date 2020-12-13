# STdata

[![R build status](https://github.com/lmweber/STdata/workflows/R-CMD-check/badge.svg)](https://github.com/lmweber/STdata/actions)

Collection of spatially resolved transcriptomics datasets stored in Bioconductor object formats. These datasets are designed for use in our online textbook "Orchestrating Spatial Transcriptomics Analysis with Bioconductor" (OSTA), as well as for other examples, demonstrations, and teaching purposes.

Datasets are stored as objects from the [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment) class. This class is an extension of `SingleCellExperiment`, adapted for spatially resolved transcriptomics data.


## Installation

Currently, the package is available from GitHub. This will be moved to Bioconductor's ExperimentHub in the future.

To install from GitHub including vignettes, use the following options:

```
remotes::install_github("lmweber/STdata", build_vignettes = TRUE, force = TRUE)
```


## How to view vignette

To view the vignette from your local installation, use the following R code:

```
vignette("STdata_overview", package = "STdata")
```

