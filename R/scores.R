#' @title Score calculation from guidance objects
#' @description Calculate scores based on residue pair score inferred by HoT, guidance and guidance2
#' @param guidanceX object of class, see \code{\link{guidance}}
#' @param score A character string indicating a type of score, currently
#'   available \code{"column"}, \code{"residue"}, \code{"alignment"},
#'   \code{"sequence"}, \code{"all"}.
#' @param na.rm Logical, indicating if NA should be removed.
#' @details The score 'column' is the GUIDANCE column score which is the mean of
#'   the residue pair residue score across columns. The score 'alignment' is the
#'   mean across the residue pair residue scores. The score 'sequence' is the
#'   mean of the residue pair score across rows (sequences). The score 'residue'
#'   is the mean score across the residue pairs with that residue (residue pair
#'   score).
#' @details The GUIDANCE column score can be utilized to weight characters in
#'   RAxML (flag -a). Simple removal of sites from the MSA should be done with
#'   cautions (Tan et al. 2015).
#' @references Penn et al. (2010). An alignment confidence score capturing
#'   robustness to guide tree uncertainty. Molecular Biology and Evolution
#'   27:1759--1767.
#' @references Tan et al. (2015). Current methods for automated filtering of
#'   multiple sequence alignments frequently worsen single-gene phylogenetic
#'   inference. Systematic biology 64:778--791.
#' @return data.frame or list of data.frames with scores
#' @author Franz-Sebastian Krah
#' @import foreach
#' @import parallel
#' @importFrom utils combn globalVariables
#' @export

scores <- function(guidanceX,
                   score = c("alignment", "column", "residue", "sequence", "column_raxml"),
                   na.rm = TRUE){

  if (!inherits(guidanceX, c("guidanceDNA", "guidanceAA"))){
    stop("guidance not of class 'guidance'")
  }
  
  ## declare i to be a global variable; this is necessary because
  ## foreach uses non-standard evaluation that codetools is not aware of
  ## [http://r.789695.n4.nabble.com/R-CMD-check-and-foreach-code-td4687660.html]
  ## Does not work [CH 2018-01-23]
  globalVariables('i')

  sc <- guidanceX@scores
  base_msa <- guidanceX@msa

  if (score == "all")
    score <- c("alignment", "column", "sequence", "residue")

  if ("alignment" %in% score){
    ## calculate GUIDANCE score
    alignment <- mean(sc, na.rm = TRUE)
  }

  if ("column" %in% score){
    ## calculate GUIDANCE score
    column <- colMeans(sc, na.rm = TRUE)
    column <- data.frame(col = 1:length(column), score = column)
    if (na.rm){
      column <- column[!is.na(column$score), ]
    }
  }

  if("column_raxml" %in% score){
    ## calculate GUIDANCE score
    column <- colMeans(sc, na.rm = TRUE)
    column <- data.frame(col = 1:length(column), score = column)
    if (na.rm){
      column <- column[!is.na(column$score), ]
    }
    column_raxml <- floor(column$score*10)
    column_raxml[is.na(column_raxml)] <- max(column_raxml, na.rm = TRUE)
  }
  
  if ("sequence" %in% score){
    ## calculate GUIDANCE score
    fac <- apply(combn(nrow(base_msa), 2), 2, paste, collapse = "-")
    fac_list <- foreach(i = 1:nrow(base_msa)) %do% grep(paste0(i, "\\b"), fac)

    sequence <- rowMeans(sc, na.rm = TRUE)
    seq <- foreach(i = 1:length(fac_list), .combine = "c") %do% {
      mean(sequence[fac_list[[i]]], na.rm = TRUE)
    }
    sequence <- data.frame(seq = 1:length(seq), score = seq)
  }

  if ("residue" %in% score){
    ## Calculate residue pair residue score
    fac <- apply(combn(nrow(base_msa), 2), 2, paste, collapse = "-")
    fac_list <- foreach(i = 1:nrow(base_msa)) %do% grep(paste0(i, "\\b"), fac)

    residue <- mclapply(fac_list, function(x) {
      colMeans(sc[x, ], na.rm = TRUE)
    }, mc.cores = detectCores())

    residue <- do.call(rbind, residue)
    colnames(residue) <- 1:ncol(residue)
    rownames(residue) <- rownames(base_msa)
    if (na.rm){
      residue <- residue[, !apply(residue, 2, function(x) !any(!is.na(x)))]
    }
  }

  return(mget(score))
}
