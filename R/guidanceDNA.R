#' @title Create Objects of Class "guidanceDNA"
#' @description Create objects of class
#'   \code{"\link[=guidanceDNA-class]{guidanceDNA}"} from objects of class
#'   \code{"\link{DNAbin}"} and a matrix derived from \code{\link{guidance}} or
#'   \code{\link{HoT}}.
#' @param msa An object of class \code{\link{DNAbin}}.
#' @param scores A matrix of quality scores.
#' @param method A characters string giving the method to derive the quality
#'   scores.
#' @param msa.method A character string giving the alignment method
#' @include guidanceDNA-class.R
#' @importFrom methods new
#' @author Christoph Heibl
#' @export

"guidanceDNA" <- function(msa, scores, method, msa.method){

  new("guidanceDNA",
      msa = msa,
      scores = scores,
      method = method,
      msa.method = msa.method
  )
}

## setMethod: indexing, extracting scores

setMethod("show", signature = "guidanceDNA",
          function(object){
            cat(nrow(object@msa), "DNA sequences with", object@method, "scores and ", object@msa.method, "MSA program")
          })
