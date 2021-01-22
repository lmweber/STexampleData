# STexampleData

[![R build status](https://github.com/lmweber/STexampleData/workflows/R-CMD-check/badge.svg)](https://github.com/lmweber/STexampleData/actions)

Collection of spatially resolved transcriptomics datasets in [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment) Bioconductor object format.

These datasets have been collected from various publicly available sources, and cover several technological platforms. We have formatted them into the `SpatialExperiment` Bioconductor object class to make them easier to access, to enable ourselves and others to use them for examples, demonstrations, tutorials, and other purposes.

The `SpatialExperiment` class is an extension of the widely-used `SingleCellExperiment` class, adapted for the properties of spatially resolved transcriptomics data.


## Installation

Currently, the package is available from GitHub. This will be moved to Bioconductor's ExperimentHub in the future.

To install from GitHub including vignettes, use the following options:

```
remotes::install_github("lmweber/STexampleData", build_vignettes = TRUE, force = TRUE)
```


## How to view vignette

To view the vignette from your local installation, use the following R code:

```
vignette("STexampleData_overview", package = "STexampleData")
```

