#' @title Add MSA to scores matrix (for developers)
#' @description Meant for developers only.
#' @param ref_col2res output of \code{\link{recodeMSA}} for base MSA
#' @param alt_col2res output of \code{\link{recodeMSA}} for alternative MSA
#' @param hit_mat output of \code{\link{initCounter}} for base MSA Cmatrix (\code{\link{recodeMSA}})
#' @details Based on the output of recodeMSA and initCounter the MSAs are compared
#' and scores are added to the scores matrix (initCounter matrix)
#' @author Franz-Sebastian Krah
#' @examples 
#' \dontrun{
#' ## Minimal (unrealistic) example for the score computation
#' # Base MSA
#' s <- matrix(c("-", "a", "c", "t", "g", 
#' "-", "a", "c", "t", "g", 
#' "a", "-", "c", "t", "g"), 3,5, byrow = T)
#' s <- as.DNAbin(s)
#' 
#' # Alternative MSA
#' s.alt <- matrix(c("-", "a", "c", "t", "g", "-", 
#' "a", "c", "t", "g", "-", 
#' "a", "c", "t", "g"), 3,5, byrow = T)
#' s.alt <- as.DNAbin(s.alt)
#' 
#' # Recode base MSA (Cmatrix)
#' s.recode <- recodeMSA(msa = s)
#' 
#' # Recode alt MSA (Cmatrix)
#' s.alt.recode <- recodeMSA(msa = s.alt)
#' 
#' # Initialize counter (basis for residue pair score)
#' s.counter <- initCounter(s.recode)
#' s.counter
#' # Update counter
#' s.counter <- addMSA(ref_col2res = s.recode, alt_col2res = s.alt.recode, hit_mat = s.counter)
#' }
#' @export


addMSA <- function(ref_col2res, alt_col2res,
                    hit_mat){
  
  hit_mat <- add_msa(ref_col2res$col2res, alt_col2res$col2res, 
                    alt_col2res$res2col, hit_mat)
  return(hit_mat)
}