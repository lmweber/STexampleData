# STexampleData

[![R build status](https://github.com/lmweber/STexampleData/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/lmweber/STexampleData/actions)


The `STexampleData` package contains a collection of spatially-resolved transcriptomics (SRT) datasets, which have been formatted into the [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment) Bioconductor class, for use in examples, demonstrations, and tutorials. The datasets are from several different SRT platforms and have been sourced from various publicly available sources. Several datasets include images and/or ground truth annotation labels.

The `STexampleData` package is available from [Bioconductor](https://bioconductor.org/packages/STexampleData).

A vignette containing examples and documentation is available from [Bioconductor](https://bioconductor.org/packages/STexampleData).


## Installation

The release version of the `STexampleData` package can be installed from Bioconductor:

```
install.packages("BiocManager")
BiocManager::install("STexampleData")
```

The development version can be installed from the development version of Bioconductor (see Bioconductor website for details) or from GitHub:

```
install.packages("remotes")
remotes::install_github("lmweber/STexampleData", build_vignettes = TRUE)
```


## Citation

- Righell D.\*, Weber L.M.\*, Crowell H.L.\*, Pardo B., Collado-Torres L., Ghazanfar S., Lun A.T.L., Hicks S.C.\*, and Risso D.\* (2022). *SpatialExperiment: infrastructure for spatially resolved transcriptomics data in R using Bioconductor.* [Bioinformatics](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btac299/6575443), btac299.
