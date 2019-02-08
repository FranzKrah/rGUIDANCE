#' @title Neighbor joining BP Tree from MSA
#' @description Estimate a neighbor joing tree from a multiple sequence alignment.
#' @param msa An object of class \code{\link{DNAbin}} or \code{\link{AAbin}}
#'   containing \strong{aligned} sequences of DNA or amino acids.
#' @param outgroup A character string giving the taxon name to be used for
#'   rooting; the default (\code{"auto"}) chooses the first taxon.
#' @return tree of class phylo
#' @importFrom ape compute.brlen nj multi2di root
#' @importFrom phangorn as.phyDat dist.ml
#' @author Franz-Sebastian Krah
#' @noRd

msaBP_nj_tree <- function(msa, outgroup = "auto"){

  base.msa.bp <- msaBP(msa)
  # convert to class phyDAT
  base.msa.ml <- as.phyDat(as.character(base.msa.bp))
  # find ML distance as input to nj tree search
  ml.dist.msa <- dist.ml(base.msa.ml)
  # NJ
  tr <- nj(ml.dist.msa)
  # root
  if (outgroup == "auto")
    outgroup <- tr$tip.label[1]
  tr <- root(tr, outgroup = outgroup)
  # resolve polytomies
  tr <- multi2di(tr)
  ## Rescale branch lengths
  tr <- compute.brlen(tr)

  return(tr)
}
