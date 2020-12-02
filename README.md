# SpaData

[![R build status](https://github.com/lmweber/SpaData/workflows/R-CMD-check/badge.svg)](https://github.com/lmweber/SpaData/actions)

Collection of spatial transcriptomics datasets stored in Bioconductor object formats. These datasets are designed for use in our online textbook "Orchestrating Spatial Transcriptomics Analysis with Bioconductor" (OSTA), as well as for other examples, demonstrations, and teaching purposes.

Datasets are currently stored as `SingleCellExperiment` objects. This will be updated to use the [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment) class.


## Vignette

To view the vignette while the package is still on GitHub, install the package with the options `build_vignettes = TRUE` and `force = TRUE`, then load the vignette locally:

```
devtools::install_github("lmweber/SpaData", build_vignettes = TRUE, force = TRUE)

vignette("SpaData_overview", package = "SpaData")
```

