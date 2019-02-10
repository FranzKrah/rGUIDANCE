#' @title Cmatrix calculation (for developers)
#' @description Meant for developers only.
#' @param msa A MSA of class \code{\link{DNAbin}} or \code{\link{AAbin}}
#' @details Here implemented as in GUIDANCE program msa_set_score. Characters are represented by odd numbers 
#' and gaps by even numbers. 
#' For more info about the Cmatrix recoding see page 3 on the supplemantary material of: Satija et al. (2009)
#' http://www.biomedcentral.com/content/supplementary/1471-2148-9-217-s1.pdf
#' @references Satija R, Novak A., Mikls I., Lyngs R., and Hein J. (2009) BigFoot: 
#' Bayesian alignment and phylogenetic footprinting with MCMC, 
#' \emph{BMC Evolutionary Biology} \strong{9}:217
#' @details For an example see \code{\link{addMSA}}
#' @references Penn et al. 2010. An
#'   alignment confidence score capturing robustness to guide tree uncertainty.
#'   \emph{Molecular Biology and Evolution} \strong{27}:1759--1767.
#' @return a list with integers
#' @author Franz-Sebastian Krah
#' @export

recodeMSA <- function(msa){
  
  gbbin <- function(msa){
    msa <- (msa != "-") * 1
    return(msa)
  }
  
  ## Cmatrix REF
  msa <- gbbin(as.character(msa))
  cmat_msa <- msa_recode(msa)
  
  return(cmat_msa)
  
}
