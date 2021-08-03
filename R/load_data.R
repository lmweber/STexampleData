#' Load data
#' 
#' Temporary function to load objects in devel version of STexampleData
#' 
#' Temporary function to load objects in devel version of STexampleData package
#' 
#' 
#' @param dataset_name Name of dataset. Options are "Visium_humanDLPFC", 
#'   "Visium_mouseCoronal", or "seqFISH_mouseEmbryo".
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
#' #spe <- load_data("Visium_humanDLPFC")
#' 
load_data <- function(dataset_name, path = "tmp", delete = TRUE) {
  
  match.arg(dataset_name, choices = c("Visium_humanDLPFC", "Visium_mouseCoronal", "seqFISH_mouseEmbryo"))
  
  # note: change "dl=0" to "dl=1" in link from Dropbox website to enable downloading
  if (dataset_name == "Visium_humanDLPFC") {
    url <- "https://www.dropbox.com/s/p7ag06lervrhqb9/Visium_humanDLPFC.rds?dl=1"
  } else if (dataset_name == "Visium_mouseCoronal") {
    url <- "https://www.dropbox.com/s/3oq67hl5zo28em6/Visium_mouseCoronal.rds?dl=1"
  } else if (dataset_name == "seqFISH_mouseEmbryo") {
    url <- "https://www.dropbox.com/s/n650cvudrfvgaic/seqFISH_mouseEmbryo.rds?dl=1"
  }
  
  fn <- file.path(path, gsub("\\?.*$", "", basename(url)))
  
  dir.create(path)
  download.file(url, fn, mode = "wb")
  spe <- readRDS(fn)
  
  if (delete) {
    unlink(path, recursive = TRUE)
  }
  
  spe
}

