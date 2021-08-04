# STexampleData

The `STexampleData` package provides access to several spatially resolved transcriptomics datasets, which have been formatted into the [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment) Bioconductor class.

These datasets have been collected from various publicly available sources, and cover several technological platforms.

We provide these datasets as `SpatialExperiment` objects to make them easier to access, so that we and others can use them for examples, demonstrations, tutorials, and other purposes.

For more details, see the [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment) documentation or our [paper](https://www.biorxiv.org/content/10.1101/2021.01.27.428431v1).


## Installation

The latest version of the package can be installed from the `devel` version of Bioconductor:

```
BiocManager::install("STexampleData", version = "devel")
```

Alternatively, you can also install the latest version from GitHub:

```
remotes::install_github("lmweber/STexampleData", ref = "no_accessors", build_vignettes = TRUE)
```


## Citation

[Righell D.\*, Weber L.M.\*, Crowell H.L.\*, Pardo B., Collado-Torres L., Ghazanfar S., Lun A.T.L., Hicks S.C.<sup>+</sup>, and Risso D.<sup>+</sup> (2021), *SpatialExperiment: infrastructure for spatially resolved transcriptomics data in R using Bioconductor*, bioRxiv](https://www.biorxiv.org/content/10.1101/2021.01.27.428431v1).

