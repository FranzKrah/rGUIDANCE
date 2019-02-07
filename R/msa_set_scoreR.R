#' @title Compare reference MSAs with alternative MSAs
#' @description MSA reliability scores (Penn et al. 2010)
#' @param ref An object of class \code{\link{DNAbin}} or \code{\link{AAbin}}
#'   containing \strong{aligned} sequences of DNA or amino acids.
#' @param alt single MSA or list of MSAs or path to alternative files. Single
#'   MSAs and list members should be of class \code{\link{DNAbin}} or
#'   \code{\link{AAbin}}.
#' @return matrix containing following scores:
#' @return residue_pair_score: if one alternative MSA is supplied then the score
#'   is 1 if a residue pair was identically aligned as in the reference MSA and
#'   0 otherwise. If more than one alternative MSA is used, then the average of
#'   the residue pair scores from all comparisons
#' @references Penn et al. 2010. An alignment confidence score capturing
#'   robustness to guide tree uncertainty. \emph{Molecular Biology and
#'   Evolution} \strong{27}:1759--1767.
#' @author Franz-Sebastian Krah
#' @export


msa_set_score <- function(ref, alt){

  ## functions 'msa_recode', 'res_pair_hit' and 'add_msa'
  ## are Rcpp functions
  ## => folder 'src' of the package source code

  if (!inherits(ref, c("DNAbin", "AAbin")))
    stop("MSA not of class DNAbin or AAbin (ape)")


  gbbin <- function(msa){
    msa <- (msa != "-") * 1
    return(msa)
  }

  ## Cmatrix REF
  ref <- gbbin(as.character(ref))
  cmat_msa <- msa_recode(ref)
  ref_col2res <- cmat_msa$col2res

  ## Hit matrix
  scores <- res_pair_hit(cmat_msa$col2res)

  ## read alternative MSAs if they are dir stored
  if (is.character(alt)){
    alt <- list.files(alt, full.names = TRUE)
    alt <- lapply(alt, read.fas)
    n <- length(alt)
  }
  if (!inherits(alt, c("list", "AAbin", "DNAbin")))
    stop("ALT MSAs must be a list of MSAs of class AAbin or DNAbin")

  if (inherits(alt, "list"))
    n <- length(alt)
  if (inherits(alt, c("AAbin", "DNAbin")))
    n <- 1

  ## Compare MSAs
  for (i in 1:n){

    ## Cmatrix MSA
    if (inherits(alt, "list"))
      com <- gbbin(as.character(alt[[i]]))
    if (inherits(alt, c("AAbin", "DNAbin")))
      com <- gbbin(as.character(alt))
    cmat_alt <- msa_recode(com)
    alt_col2res <- cmat_alt$col2res
    alt_res2col <- cmat_alt$res2col

    ## Compare
    ## function overrides input *scores* internally
    scores <- add_msa(ref_col2res, alt_col2res, alt_res2col, hit_mat = scores)

  }
  scores[scores == -1] <- NA
  scores / n
}
