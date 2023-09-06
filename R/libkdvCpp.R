#' @useDynLib Rlibkdv
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("Rlibkdv", libpath)
}
