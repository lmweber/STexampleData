#' Load data
#' 
#' Temporary function to load data objects from Dropbox
#' 
#' Temporary function to load data objects from Dropbox until ExperimentHub
#' package is set up
#' 
#' 
#' @param dataset_name Name of dataset. Either "human_DLPFC" or "mouse_coronal".
#' 
#' @param path Directory to save temporary files. Default = "tmp".
#' 
#' @param delete Whether to delete temporary files. Default = TRUE.
#' 
#' 
#' @return SingleCellExperiment or SpatialExperiment object
#' 
#' 
#' @importFrom utils download.file
#' 
#' @export
#' 
#' @examples
#' #spe <- load_data("human_DLPFC")
#' #spe <- load_data("mouse_coronal")
#' 
load_data <- function(dataset_name, path = "tmp", delete = TRUE) {
  
  match.arg(dataset_name, choices = c("human_DLPFC", "mouse_coronal"))
  
  # note: change "dl=0" to "dl=1" in link copied from Dropbox website to allow 
  # downloading from command line
  if (dataset_name == "human_DLPFC") {
    url <- "https://www.dropbox.com/s/b3s04x7il31g5qy/human_DLPFC_151673.RData?dl=1"
  } else if (dataset_name == "mouse_coronal") {
    url <- "https://www.dropbox.com/s/b3cu79rk7dkupa8/mouse_coronal.RData?dl=1"
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

