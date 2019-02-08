setOldClass("AAbin")

#' @title An S4 Class to represent Amino Acid Alignments with Quality Scores
#' @description \code{"guidanceAA"} holds a multiple sequence alignment of alino
#'   acids together with quality scores for each cell in the alignment.
#' @slot msa An object of class \code{\link{AAbin}}.
#' @slot scores A matrix of quality scores.
#' @slot method A characters string giving the method to derive the quality
#'   scores.
#' @slot msa.method A character string giving the alignment method
#' @seealso \code{"\link[=guidanceDNA-class]{guidanceDNA}"}
#' @author Christoph Heibl


setClass("guidanceAA",
  representation = list(
    msa = "AAbin",
    scores = "matrix",
    method = "character", 
    msa.method = "character")
)
