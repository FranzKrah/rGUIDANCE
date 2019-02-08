#' @title Filtering of the base MSA using guidance scores
#'
#' @param guidanceX object of class \code{\link{guidance}}
#' @param col.cutoff numeric between 0 and 1; removes unreliable columns below
#'   the cutoff (default: 0.2); ignored if FALSE
#' @param seq.cutoff numeric between 0 and 1; removes unreliable sequences below
#'   the cutoff (default: 0.1); ignored if FALSE
#' @param mask.cutoff residues below the cutoff are masked ('N' for DNA, 'X' for
#'   AA; default: 0.5); ignored if FALSE
#' @param filter.ends logical, if TRUE trim.ends (ips) is applied to the MSA
#' @param filter.gaps logical, if TRUE trim.gabs (ips) is applied to the MSA
#' @param na.coding value with which to replace NAs
#' 
#' @details Please note that we do not recommend filtering MSAs (Tan et al. 2015). Instead the user should use the 
#' GUIDANCE column score as individual weights to each column of the alignment to RAxMLvia the –a flag 
#' (for a practical example: Krah et al. (2018)). This can be achieved by using the function scores with score = "column_raxml"
#' and passing the score as weights to \code{\link{raxml}}. For a detailled example, see Vignette.
#' 
#' @references Tan et al. (2015). Current methods for automated filtering of
#'   multiple sequence alignments frequently worsen single-gene phylogenetic
#'   inference. \emph{Systematic biology} \strong{64}:778--791.
#' @references Krah et al. (2018). Evolutionary dynamics of host specialization in wood-decay fungi. 
#' \emph{BMC Evolutionary Biology} \strong{18}:119-132
#' 
#' @return filtered MSA of class \code{AAbin} or \code{DNAbin}
#' @return if filter.ends or filter.gaps was choosen, the adjusted 
#' @seealso \code{\link{scores}}
#' @author Franz-Sebastian Krah
#' @export

filterMSA <- function(guidanceX,
                      col.cutoff = 0.2,
                      seq.cutoff = 0.1,
                      mask.cutoff = 0.5,
                      filter.ends = FALSE,
                      filter.gaps  = FALSE,
                      na.coding = 0.5) {
  
  base_msa <- guidanceX@msa
  
  if (!mask.cutoff == FALSE) {
    r_sc <- scores(guidanceX, score = "residue", na.rm = FALSE)
    
    if (inherits(base_msa, "AAbin")) {
      base_msa <- as.character(base_msa)
      base_msa[r_sc$residue < mask.cutoff & base_msa != "-"] <- "X"
      base_msa <- as.AAbin(base_msa)
    }
    if (inherits(base_msa, "DNAbin")) {
      base_msa <- as.character(base_msa)
      base_msa[r_sc$residue < mask.cutoff & base_msa != "-"] <- "N"
      base_msa <- as.DNAbin(base_msa)
    }
  }
  
  if (!col.cutoff == FALSE) {
    g_sc <- scores(guidanceX, score = "column", na.rm = FALSE)
    # base_msa <- as.character(base_msa)
    
    base_msa <- base_msa[, c(which(g_sc$column$score >= col.cutoff), which(is.na(g_sc$column$score)))]
    g_sc <- g_sc$column[c(which(g_sc$column$score >= col.cutoff), which(is.na(g_sc$column$score))), ]
  }
  
  if (!seq.cutoff == FALSE) {
    s_sc <- scores(guidanceX, score = "sequence")
    base_msa <- base_msa[s_sc$sequence$score >= seq.cutoff, ]
  }
  
  if (filter.ends) {
    keep <- trimEnds(base_msa)[[2]]
  }
  
  if (filter.gaps) {
    keep2 <- deleteGaps(base_msa)[[2]]
  }
  
  if (filter.ends | filter.gaps) {
    keep <- mget(c("keep", "keep2"), ifnotfound = list(NULL, NULL))
    keep <- Reduce(intersect, keep)
    base_msa <- base_msa[, keep]
    
    # g_sc <- scores(guidanceX, score = "column", na.rm = FALSE)
    # g_sc <- g_sc$column[keep, ]
    # g_sc[is.na(g_sc$score), ]$score <- na.coding
  }
  return(base_msa)
}
