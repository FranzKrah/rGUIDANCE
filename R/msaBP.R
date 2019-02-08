#' @title Bootstrap a Multiple Sequence Alignnment
#' @description Create a bootstrap replicate of multiple sequence alignment.
#' @param msa An object of class \code{\link{DNAbin}} or \code{\link{AAbin}}.
#' @author Franz-Sebastian Krah
#' @noRd

msaBP <- function(msa){

  if (!inherits(msa, c("DNAbin", "AAbin")))
    stop("'msa' not of class DNAbin or AAbin (ape)")

  msa[, sample(ncol(msa), replace = TRUE)]
}

