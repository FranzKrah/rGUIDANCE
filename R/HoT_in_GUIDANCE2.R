#' @title HoT within guidance2
#' @description Internally Heads or Tails Alignment Reliability within guidance2
#' @param msa An object of class \code{\link{DNAbin}} or
#'   \code{\link{AAbin}} containing aligned sequences of DNA or amino acids.
#' @param n.coopt An integer giving the number of sampled co-optimal MSAa.
#' @param raw_seq An object of class \code{\link{DNAbin}} or
#'   \code{\link{AAbin}} containing unaligned sequences of DNA or amino acids.
#' @param method A character string containing further arguments passed to
#'   MAFFT; default is \code{"auto"}.
#' @param msa.exec A character string giving the path to the executable of the
#'   alignment program (e.g. \code{/usr/local/bin/mafft}); possible programs are
#'   \code{MAFFT}, \code{MUSCLE}, and \code{ClustalW}.
#' @author Franz-Sebastian Krah
#' @importFrom ape compute.brlen multi2di nj Ntip root
#' @import foreach
#' @importFrom ips read.fas
#' @importFrom utils globalVariables
#' @noRd

Hot_GUIDANCE2 <- function(msa, n.coopt, 
                          raw_seq, method,
                          msa.exec){
  
  ## declare i to be a global variable; this is necessary because
  ## foreach uses non-standard evaluation that codetools is not aware of
  ## [http://r.789695.n4.nabble.com/R-CMD-check-and-foreach-code-td4687660.html]
  ## Does not work [CH 2018-01-23]
  # globalVariables("i")
  
  ## look up MSA program specified
  msa.program <- str_extract(msa.exec, "mafft|muscle|clustal\\w")
  
  # create start_tree
  li <- msa
  if (is.character(msa))
    msa <- read.fas(msa)
  msa.nj <- nj(dist.ml(as.phyDat(msa)))
  msa.nj <- root(msa.nj, outgroup = msa.nj$tip.label[1])
  msa.nj <- multi2di(msa.nj)
  msa.nj <- compute.brlen(msa.nj)
  
  ## produce MSA partitions
  align_parts <- partitions(msa.nj)
  
  ## sample 4 or n co-optimal
  nt <- Ntip(msa.nj)
  n.co <- sample((nt - 3) * 8, n.coopt)
  n.co_mat <- data.frame(n = 1:((nt - 3) * 8),
                         part = rep(1:(nt - 3), each = 8),
                         n.in.part = rep(1:8, times = (nt - 3)))
  n.co.sub <- n.co_mat[n.co_mat$n %in% n.co, ]
  
  # reduce partitions to the randomly choosen co-optimal MSA number
  align_parts <- align_parts[, n.co.sub$part]
  # number of random MSA within partition (remember, each partition has 8 alignments)
  n.co.sub <- n.co.sub$n.in.part
  
  # make the 4 or n alignments
  alt_msas <- foreach(i = 1:ncol(align_parts),
                      .export = c("align_parts")) %do% {
                        align_part_set(x = raw_seq, partition_set = align_parts[, i],
                                       method = method, msa.exec = msa.exec,
                                       coopt.sub = n.co.sub[i])
                      }
  
  ## unlist
  # alt_msas <- foreach(i = seq_along(alt_msas), .combine = "c") %do% {
  #   alt_msas[[i]]
  # }
  alt_msas <- unlist(alt_msas, recursive = FALSE)
  
  ## write to files
  return(alt_msas)
}


