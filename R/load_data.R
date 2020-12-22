#' Load data
#' 
#' Function to load data objects from STdata package
#' 
#' Function to load data objects from Dropbox links (until package is moved to
#' Bioconductor's ExperimentHub).
#' 
#' 
#' @param dataset_name Name of dataset. Currently accepts "human_DLPFC",
#'   "mouse_coronal", or "seqFISH_mouseEmbryo".
#' 
#' @param path Directory to save temporary files. Default = "tmp".
#' 
#' @param delete Whether to delete temporary files. Default = TRUE.
#' 
#' 
#' @return SpatialExperiment object
#' 
#' 
#' @importFrom utils download.file
#' 
#' @export
#' 
#' @examples
#' #spe <- load_data("human_DLPFC")
#' #spe <- load_data("mouse_coronal")
#' #spe <- load_data("seqFISH_mouseEmbryo")
#' 
load_data <- function(dataset_name, path = "tmp", delete = TRUE) {
  
  match.arg(dataset_name, choices = c("human_DLPFC", "mouse_coronal", "seqFISH_mouseEmbryo"))
  
  # note: change "dl=0" to "dl=1" in link from Dropbox website to enable downloading
  if (dataset_name == "human_DLPFC") {
    url <- "https://www.dropbox.com/s/f26li36ov1aebl8/human_DLPFC.RData?dl=1"
  } else if (dataset_name == "mouse_coronal") {
    url <- "https://www.dropbox.com/s/qpbjgcl4o5q29b7/mouse_coronal.RData?dl=1"
  } else if (dataset_name == "seqFISH_mouseEmbryo") {
    url <- "https://www.dropbox.com/s/bit9odxqw6ugbuf/seqFISH_mouseEmbryo.RData?dl=1"
  }
  
  fn <- file.path(path, gsub("\\?.*$", "", basename(url)))
  
  dir.create(path)
  download.file(url, fn, mode = "wb")
  load(fn)
  
  if (delete) {
    unlink(path, recursive = TRUE)
  }
  
  spe
}

