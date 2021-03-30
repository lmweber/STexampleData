\name{STexampleData}
\docType{data}
\alias{STexampleData}
\alias{STexampleData-package}
\alias{Visium_humanDLPFC}
\alias{Visium_mouseCoronal}
\alias{seqFISH_mouseEmbryo}
\title{
Collection of spatially resolved transcriptomics datasets in SpatialExperiment 
Bioconductor format
}
\description{
Collection of spatially resolved transcriptomics datasets in SpatialExperiment 
Bioconductor format, for use in examples, demonstrations, tutorials, and other 
purposes. The datasets have been sourced from various publicly available 
sources, and cover several technological platforms.
}
\details{
The \code{STexampleData} package provides access to a collection of spatially 
resolved transcriptomics datasets, which have been formatted into the 
\code{SpatialExperiment} Bioconductor object class.

These datasets have been collected from various publicly available sources, and 
cover several technological platforms. We provide them in the form of 
\code{SpatialExperiment} objects to make them easier to access, so that we and 
others can use them for examples, demonstrations, tutorials, and other purposes.

The \code{SpatialExperiment} class is an extension of \code{SingleCellExperiment}, 
adapted for the properties of spatially resolved transcriptomics data. For more 
details, see the \code{SpatialExperiment} documentation.

\emph{Datasets}

The \code{STexampleData} package contains the following datasets:

\itemize{
\item Visium_humanDLPFC (10x Genomics Visium): A single sample (sample 151673) 
of human brain dorsolateral prefrontal cortex (DLPFC) in the human brain, 
measured using the 10x Genomics Visium platform. This is a subset of the full 
dataset containing 12 samples from 3 neurotypical donors, published by Maynard 
and Collado-Torres et al. (2021). The full dataset is available from the 
\code{spatialLIBD} Bioconductor package.

\item Visium_mouseCoronal (10x Genomics Visium): A single coronal section from 
the mouse brain spanning one hemisphere, measured using the 10x Genomics Visium 
platform. This dataset was publicly released by 10x Genomics.

\item seqFISH_mouseEmbryo (seqFISH): A subset of cells (embryo 1, z-slice 2) 
from a dataset investigating mouse embryogenesis by Lohoff and Ghazanfar et al. 
(2020), generated using the seqFISH platform. The full dataset is available from 
the original publication.
}
}
