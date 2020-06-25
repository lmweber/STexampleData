#' Load data
#' 
#' Temporary function to load data objects from Dropbox
#' 
#' Temporary function to load data objects from Dropbox before ExperimentHub
#' package is set up
#' 
#' 
#' @param data_name Name of dataset. Currently accepts "human_DLPFC".
#' 
#' @param path Directory to save temporary file. Default = current working
#'   directory.
#' 
#' @param delete Whether to delete temporary file. Default = TRUE.
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
#' spe <- load_data("human_DLPFC")
#' 
load_data <- function(data_name, path = ".", delete = TRUE) {
  
  match.arg(data_name, choices = c("human_DLPFC"))
  
  if (data_name == "human_DLPFC") {
    url <- "https://www.dropbox.com/s/b3s04x7il31g5qy/human_DLPFC_151673.RData?dl=1"
    filename <- file.path(path, "human_DLPFC_151673.RData")
    download.file(url, filename, mode = "wb")
    load(filename)
    if (delete) unlink(filename)
  }
  
  spe
  
}

