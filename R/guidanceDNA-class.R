setOldClass("DNAbin")

#' @title An S4 Class to represent DNA Alignments with Quality Scores
#' @description \code{guidanceDNA} holds a multiple sequence alignment of DNA
#'   together with quality scores for each cell in the alignment.
#' @slot msa An object of class \code{\link{DNAbin}}.
#' @slot scores A matrix of quality scores.
#' @slot method A characters string giving the method to derive the quality
#'   scores.
#' @slot msa.method A character string giving the alignment method
#' @seealso \code{"\link[=guidanceAA-class]{guidanceAA}"}
#' @author Christoph Heibl

setClass("guidanceDNA",
         representation = list(
           msa = "DNAbin",
           scores = "matrix",
           method = "character",
           msa.method = "character")
)

