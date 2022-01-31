# STexampleData

The `STexampleData` package contains a collection of spatially resolved transcriptomics (ST) datasets, which have been formatted into the [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment) Bioconductor class, for use in examples, demonstrations, and tutorials. The datasets are from several different ST platforms and have been sourced from various publicly available sources. Some of the datasets include images and/or ground truth annotation labels.

The `STexampleData` package is available from [Bioconductor](https://bioconductor.org/packages/STexampleData).

A vignette containing examples and documentation is available from [Bioconductor](https://bioconductor.org/packages/STexampleData).


## Installation

The `STexampleData` package can be installed from Bioconductor. Note that Bioconductor follows a "release" and "development" schedule, where the release version is considered to be stable and updated every 6 months, and the development version contains latest updates (and then becomes the next release version every 6 months).


### Release version

To install the stable release version, install the latest release version of R from [CRAN](https://cran.r-project.org/), then install the Bioconductor package installer and the `STexampleData` package as follows.

```
install.packages("BiocManager")
BiocManager::install("STexampleData")
```

### Development version

To install the development version, there are two options.

(i) Install the [appropriate version of R (R-release between April and October, or R-devel between October and April)](http://bioconductor.org/developers/how-to/useDevel/), then install the Bioconductor package installer and the development version of the `STexampleData` package as follows.

```
install.packages("BiocManager")
BiocManager::install("STexampleData", version = "devel")
```

(ii) Alternatively, if you do not want to install R-devel (since this may cause issues with other packages), you can use the latest release version of R from [CRAN](https://cran.r-project.org/) and install the development version of the `STexampleData` package from GitHub as follows.

```
install.packages("remotes")
remotes::install_github("lmweber/STexampleData", build_vignettes = TRUE)
```


## Citation

Righell D.\*, Weber L.M.\*, Crowell H.L.\*, Pardo B., Collado-Torres L., Ghazanfar S., Lun A.T.L., Hicks S.C.<sup>+</sup>, and Risso D.<sup>+</sup> (2021), *SpatialExperiment: infrastructure for spatially resolved transcriptomics data in R using Bioconductor*, [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.01.27.428431v2).

