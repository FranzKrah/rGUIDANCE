#' @title Reverse DNA Sequences of class DNAbin or AAbin
#' @description Creates the reverse of DNA of amino acid sequences.
#' @param x An object of class \code{\link[ape]{DNAbin}} or \code{\link[ape]{AAbin}}.
#' @importFrom ape as.AAbin as.DNAbin
#' @export

rev_DNA <- function(x){

  if (!inherits(x, c("DNAbin","AAbin")))
    stop("'x' not of class DNAbin or AAbin")

  if (is.matrix(x)){
    x_mat <- as.character(x)
    rev <- t(apply(x_mat, 1, rev))
  }
  if (is.list(x)){
    x_l <- as.character(x)
    rev <- lapply(x_l, rev)
  }
  if (inherits(x, "DNAbin"))
    rev <- as.DNAbin(rev)
  if (inherits(x, "AAbin"))
    rev <- as.AAbin(rev)
  return(rev)
}
