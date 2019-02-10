#' @title Initialize counters matrix (for developers)
#' @description Meant for developers only.
#' @param col2res Output from \code{\link{recodeMSA}}
#' @details Produce res_pair_hit object as basis for comparing residue pairs.
#' This matrix contains all residue pairs in the MSA and the information
#' if they are bases (hit = 0) or gabs (-1); the hits are later updated by
#' comparisons with the alternative MSAs (see add_msa).
#' The function 'init_counters' in 'set_msa_score' (see \code{\link{recodeMSA}}) uses 3D matrices,
#' which is not straight forewardly implemented in Rcpp, which is why I
#' work with a nxk matrix, with k many columns as the REF MSA and n rows for all
#' residue pairs.
#' This is almost equally fast, however, the subsequent calculation of derived
#' scores is not as handy. These are done in the function 'scores'.
#' @details For an example see \code{\link{addMSA}}
#' @author Franz-Sebastian Krah
#' @export

initCounter <- function(col2res){
  
  scores <- res_pair_hit(col2res$col2res)
  return(scores)
  
}