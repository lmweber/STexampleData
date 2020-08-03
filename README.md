# STdata

Collection of spatial transcriptomics datasets stored in Bioconductor object formats. These datasets are designed for use in our online textbook "Orchestrating Spatial Transcriptomics Analysis with Bioconductor" (OSTA), as well as for other examples, demonstrations, and teaching purposes.

Datasets are currently stored as `SingleCellExperiment` objects.


## Vignette

To view the vignette while the package is on GitHub, install the package with the options `build_vignettes = TRUE` and `force = TRUE`, then load the vignette locally:

```
devtools::install_github("lmweber/STdata", build_vignettes = TRUE, force = TRUE)

vignette("STData_overview", package = "STdata")
```

