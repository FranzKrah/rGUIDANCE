#' @title Create Objects of Class "guidanceAA"
#' @description Create objects of class
#'   \code{"\link[=guidanceAA-class]{guidanceAA}"} from objects of class
#'   \code{"\link{AAbin}"} and a matrix derived from \code{\link{guidance}} or
#'   \code{\link{HoT}}.
#' @param msa An object of class \code{\link{AAbin}}.
#' @param scores A matrix of quality scores.
#' @param method A characters string giving the method to derive the quality
#'   scores.
#' @include guidanceAA-class.R
#' @importFrom methods new
#' @author Christoph Heibl
#' @export

"guidanceAA" <- function(msa, scores, method){

  new("guidanceAA",
    msa = msa,
    scores = scores,
    method = method
  )
}

## setMethod: indexing, extracting scores,

setMethod("show", signature = "guidanceAA",
  function(object){
    cat(nrow(object@msa), "AA sequences with", object@method, "scores")
  })
