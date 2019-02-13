#' @title Score Calculation from Guidance Objects
#' @description Calculate scores based on residue pair score inferred by HoT,
#'   guidance and guidance2
#' @param guidanceX object of class, see \code{\link{guidance}}
#' @param score A character string indicating a type of score, currently
#'   available are \code{"column"}, \code{"residue"}, \code{"alignment"},
#'   \code{"sequence"}, \code{"all"}, or \code{"column_raxml"}
#' @param na.rm Logical, indicating if NAs should be removed. If
#'   \code{"column_raxml"} is used, then this should be FALSE
#' @param na.raxml Integer, specifying what should be used if column score is
#'   NA, default is 10 (see Details).
#' @param flag_a character specifying a path. If path is supplied function
#'   writes the MSA into a fasta file. Additionally the function produces a file
#'   with the column score ready for RAxML input (flag -a) outside of R
#' @details The score 'column' is the GUIDANCE column score which is the mean of
#'   the residue pair residue score across columns. The score 'alignment' is the
#'   mean across the residue pair residue scores. The score 'sequence' is the
#'   mean of the residue pair score across rows (sequences). The score 'residue'
#'   is the mean score across the residue pairs with that residue (residue pair
#'   score).
#' @details The GUIDANCE column score can be utilized to weight characters in
#'   RAxML (flag -a). Simple removal of sites from the MSA should be done with
#'   cautions (Tan et al. 2015). Using score = "column_raxml" allows to use the
#'   output directly for ips::raxml (see Vignette for a detailed example). The
#'   column score is converted to an integer from 0 to 10. NAs are set to 10 per
#'   default, but the user may choose otherwise.
#' @references Penn et al. (2010). An alignment confidence score capturing
#'   robustness to guide tree uncertainty. \emph{Molecular Biology and
#'   Evolution} \strong{27}:1759--1767.
#' @references Tan et al. (2015). Current methods for automated filtering of
#'   multiple sequence alignments frequently worsen single-gene phylogenetic
#'   inference. \emph{Systematic Biology} \strong{64}:778--791.
#' @seealso weights in \code{\link{raxml}}
#' @return data.frame or list of data.frames with scores
#' @author Franz-Sebastian Krah
#' @importFrom foreach foreach %do%
#' @import parallel
#' @importFrom utils combn
#' @export

scores <- function(guidanceX,
                   score = c("alignment", "column", "residue", "sequence", "column_raxml"),
                   na.rm = TRUE,
                   na.raxml = 10,
                   flag_a = FALSE){

  if (!inherits(guidanceX, c("guidanceDNA", "guidanceAA"))){
    stop("guidance not of class 'guidance'")
  }
  
  ## declare i to be a global variable; this is necessary because
  ## foreach uses non-standard evaluation that codetools is not aware of
  ## [http://r.789695.n4.nabble.com/R-CMD-check-and-foreach-code-td4687660.html]
  ## Does not work: globalVariables('i') [CH 2018-01-23]

  sc <- guidanceX@scores
  base_msa <- guidanceX@msa
  
  if (flag_a)
    score <- "column_raxml"

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

  if ("column_raxml" %in% score){
    ## calculate GUIDANCE score
    column <- colMeans(sc, na.rm = TRUE)
    column <- data.frame(col = 1:length(column), score = column)
    if (na.rm){
      column <- column[!is.na(column$score), ]
    }
    column_raxml <- floor(column$score * 10)
    column_raxml[is.na(column_raxml)] <- na.raxml
  }
  
  if ("sequence" %in% score){
    ## calculate GUIDANCE score
    fac <- apply(combn(nrow(base_msa), 2), 2, paste, collapse = "-")
    fac_list <- foreach(i = 1:nrow(base_msa)) %do% { grep(paste0(i, "\\b"), fac) }

    sequence <- rowMeans(sc, na.rm = TRUE)
    seq <- foreach(i = 1:length(fac_list), .combine = "c") %do% {
      mean(sequence[fac_list[[i]]], na.rm = TRUE)
    }
    sequence <- data.frame(seq = 1:length(seq), score = seq)
  }

  if ("residue" %in% score){
    ## Calculate residue pair residue score
    fac <- apply(combn(nrow(base_msa), 2), 2, paste, collapse = "-")
    fac_list <- foreach(i = 1:nrow(base_msa)) %do% { grep(paste0(i, "\\b"), fac) }

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
  
  if (!flag_a == FALSE) {
    write.fas(base_msa,
              file = paste0(flag_a, "/filtMSA_", Sys.Date(), ".fas"))
    write(
      column_raxml,
      file = paste(flag_a, "/filtMSA_GCSC_", Sys.Date(), ".txt", sep = ""),
      sep = "\t"
    )
  } 

  return(mget(score))
}
