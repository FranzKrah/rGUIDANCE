#' Filtering of the base MSA using guidance scores
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
#' @param column_score logical, if TRUE column score (e.g. for RAxML: flag -a)
#'   is in the output
#' @param flag_a character specifying a path. If path is supplied function
#'   writes the filtered MSA into a fasta file. Additionally the function
#'   produces a file with the column score ready for RAxML input (flag -a)
#' @param na.coding value with which to replace NAs
#' @return masked MSA of class \code{AAbin} or \code{DNAbin}
#' @return column_score is optional
#' @seealso \code{\link{scores}}
#' @author Franz-Sebastian Krah
#' @export

filterMSA <- function(guidanceX,
                      col.cutoff = 0.2,
                      seq.cutoff = 0.1,
                      mask.cutoff = 0.5,
                      filter.ends = FALSE,
                      filter.gaps  = FALSE,
                      column_score = FALSE,
                      flag_a = FALSE,
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

    ## 'g_r' replaced by 'guidance' [CH-2017-11-08]  
    g_sc <- scores(guidanceX, score = "column", na.rm = FALSE)
    g_sc <- g_sc$column[keep, ]
    g_sc[is.na(g_sc$score), ]$score <- na.coding
  }

  if (!flag_a == FALSE) {
    write.fas(base_msa,
              file = paste0(flag_a, "/filtMSA_", Sys.Date(), ".fas"))
    g_sc <- paste(round(g_sc$score * 10, digits = 0), collapse = " ")
    write(
      g_sc,
      file = paste(flag_a, "/filtMSA_GCSC_", Sys.Date(), ".txt", sep = ""),
      sep = "\t"
    )
  } else{
    if (column_score) {
      return(list(msa = base_msa, column_score = g_sc))
    } else{
      return(base_msa)
    }
  }
}
