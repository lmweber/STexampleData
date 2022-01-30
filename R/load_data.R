#' Functions to load data
#' 
#' Named functions to load objects from URLs
#' 
#' Named functions to load objects from URLs in development version of
#' STexampleData package
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
#' # spe <- Visium_humanDLPFC()
#' # spe <- Visium_mouseCoronal()
#' # spe <- seqFISH_mouseEmbryo()
#' 
#' @name load_data
#' 
NULL

#' @rdname load_data
Visium_humanDLPFC <- function() {
  
  # note: change "dl=0" to "dl=1" in link from Dropbox website to enable downloading
  url <- "https://www.dropbox.com/s/hkhudyddds0c4j5/Visium_humanDLPFC.rds?dl=1"
  
  path <- "tmp"
  fn <- file.path(path, basename(url))
  dir.create(path)
  download.file(url, fn, mode = "wb")
  
  spe <- readRDS(fn)
  
  unlink(path, recursive = TRUE)
  
  spe
}

#' @rdname load_data
Visium_mouseCoronal <- function() {
  
  # note: change "dl=0" to "dl=1" in link from Dropbox website to enable downloading
  url <- "https://www.dropbox.com/s/751ru9scwnppych/Visium_mouseCoronal.rds?dl=1"
  
  path <- "tmp"
  fn <- file.path(path, basename(url))
  dir.create(path)
  download.file(url, fn, mode = "wb")
  
  spe <- readRDS(fn)
  
  unlink(path, recursive = TRUE)
  
  spe
}

#' @rdname load_data
seqFISH_mouseEmbryo <- function() {
  
  # note: change "dl=0" to "dl=1" in link from Dropbox website to enable downloading
  url <- "https://www.dropbox.com/s/sfes5p3goxgr7zd/seqFISH_mouseEmbryo.rds?dl=1"
  
  path <- "tmp"
  fn <- file.path(path, basename(url))
  dir.create(path)
  download.file(url, fn, mode = "wb")
  
  spe <- readRDS(fn)
  
  unlink(path, recursive = TRUE)
  
  spe
}

