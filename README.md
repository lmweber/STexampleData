# STexampleData

Collection of spatially resolved transcriptomics datasets in [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment) Bioconductor format.

These datasets have been collected from various publicly available sources, and cover several technological platforms. We provide them in the form of `SpatialExperiment` objects to make them easier to access, so that we and others can use them for examples, demonstrations, tutorials, and other purposes.

The `SpatialExperiment` class is an extension of `SingleCellExperiment`, adapted for the properties of spatially resolved transcriptomics data. For more details, see the `SpatialExperiment` documentation.


## Installation

The package is currently available from GitHub, and will be submitted to Bioconductor's ExperimentHub.

To install the current version from GitHub:

```
remotes::install_github("lmweber/STexampleData", ref = "non_EH_load_data", build_vignettes = TRUE)
```

