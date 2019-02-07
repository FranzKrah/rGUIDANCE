#' @title Create N-3 Partitions from Phylogeny
#' @description A wrapper to \code{\link{treePart}} from the \strong{adephylo}
#'   package.
#' @param tree A phylogenetic tree of class \code{\link{phylo}}.
#' @importFrom adephylo treePart
#' @author Franz-Sebastian Krah
#' @export

partitions <- function(tree){

  # get all partitions of the tree
  partition <- treePart(tree, "dummy")

  # remove duplicated partition (first bifurcation from node)
  # if the first split includes a single sequences, then do not take it as
  # first comparision

  partition <- partition[!(colSums(partition) %in% c(nrow(partition) - 1, 1))]

  # find dupliate and remove from partitions
  ## [FK - 2017-03-27]
  ## this can be done faster if only the combinations are checked.
  dup <- list()
  res <- list()
  for (j in 1:ncol(partition)){
    part <- partition[,-j]
    for (i in 1:ncol(part)){
      dup[[i]] <- (partition[, j] > 0) == !(part[, i] > 0)
    }
    res[[j]] <- dup
  }
  dup <- lapply(res, function(x) do.call(cbind, x))
  del <- unlist(lapply(dup, function(x) which(apply(x, 2, all))))
  if (length(del))
    partition <- partition[-del[2]]
  return(partition)
}


