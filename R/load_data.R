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
#' @export
#' 
#' @examples
#' # spe <- Visium_humanDLPFC()
#' # spe <- Visium_mouseCoronal()
#' # spe <- seqFISH_mouseEmbryo()
#' # spe <- ST_mouseOB()
#' # spe <- SlideSeqV2_mouseHPC()
#' 
#' @name load_data
#' 
NULL

#' @rdname load_data
Visium_humanDLPFC <- function() {
  # note: change "dl=0" to "dl=1" in link from Dropbox website to enable downloading
  urlname <- "https://www.dropbox.com/s/tu12bjsbodwkabo/Visium_humanDLPFC.rds?dl=1"
  spe <- readRDS(url(urlname))
  spe
}

#' @rdname load_data
Visium_mouseCoronal <- function() {
  # note: change "dl=0" to "dl=1" in link from Dropbox website to enable downloading
  urlname <- "https://www.dropbox.com/s/3ch2vinhsvrp7m9/Visium_mouseCoronal.rds?dl=1"
  spe <- readRDS(url(urlname))
  spe
}

#' @rdname load_data
seqFISH_mouseEmbryo <- function() {
  # note: change "dl=0" to "dl=1" in link from Dropbox website to enable downloading
  urlname <- "https://www.dropbox.com/s/hahxamfkv4dkpeo/seqFISH_mouseEmbryo.rds?dl=1"
  spe <- readRDS(url(urlname))
  spe
}

#' @rdname load_data
ST_mouseOB <- function() {
  # note: change "dl=0" to "dl=1" in link from Dropbox website to enable downloading
  urlname <- "https://www.dropbox.com/s/zrohcp57sbyuyq0/ST_mouseOB.rds?dl=1"
  spe <- readRDS(url(urlname))
  spe
}

#' @rdname load_data
SlideSeqV2_mouseHPC <- function() {
  # note: change "dl=0" to "dl=1" in link from Dropbox website to enable downloading
  urlname <- "https://www.dropbox.com/s/qv3h54p1ulb4fxt/SlideSeqV2_mouseHPC.rds?dl=1"
  spe <- readRDS(url(urlname))
  spe
}

