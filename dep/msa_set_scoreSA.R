#' @title Compare Reference MSAs with Alternative MSAs
#' @description Wrapper function for program msa_set_score v2.01 of the GUIDANCE
#'   program (see reference). Copyright: To modify the code, or use parts of it
#'   for other purposes, permission should be requested. Please contact Tal
#'   Pupko: talp@post.tau.ac.il. Please note that the use of the GUIDANCE
#'   program is for academic use only. C code computing basic MSA comparision.
#'   The most basic is the residue pairs residue score, which checks if residue
#'   pairs combinations are correctly aligned in both MSAs. From this the
#'   residue score, residue column score (GUIDANCE score), residue sequence
#'   score are computed. It also calulates the column score (CS), which simply
#'   checks if a column is identically aligned in the alternative MSA.
#' @param ref of class data.frame, is the reference MSA ('BASE MSA') with
#'   sequences as columns
#' @param alt path to alternative files
#' @param bootstrap Number of bootstraps
#' @param exec XXX.
#' @return list containing following scores:
#' @return mean_scores: residue pair score and mean column score
#' @return column_score: identically aligned columns are asigned 1, otherwise 0;
#'   if more than one alternative MSA is supplied then mean of this score
#' @return residue_column_score: if one alternative MSA is supplied, then this
#'   is the SPC (sum-of-pairs column score); if more than one this is the
#'   average of the SPCs of all MSAs. The SPC is the mean over all residue pair
#'   scores (see there) in a column. This is also refered to as the GUIDANCE
#'   score.
#' @return residue_pair_score: if one alternative MSA is supplied then this is 1
#'   if a residue pair was identically aligned as in the reference MSA and 0
#'   otherwise. If more than one, then the average of the residue pair scores
#'   that result from each comparison with the reference MSA.
#' @return residue_pair_residue_score: calculated by averaging the residue
#'   column score over all pairs that have this residue. This is a confidence
#'   score for each residue.
#' @return residual_pair_sequence_pair_score_ if one alternative MSA is supplied
#'   this compares residue pairs of all sequence pairs of the base MSA with the
#'   reference MSA and returns 1 if identical and 0 if not and then averages
#'   over all pairs of a sequence pair. If more than one the the mean over those
#'   averages are computed.
#' @return residual_pair_sequence_score: like the residue_column_score but for
#'   sequences
#' @references Penn et al. (2010). An alignment confidence score capturing
#'   robustness to guide tree uncertainty. Molecular Biology and Evolution
#'   27:1759--1767.
#' @author Franz-Sebastian Krah
#' @seealso \code{\link{guidance}}, \code{\link{guidance2}}, \code{\link{HoT}}
#' @import plyr
#' @importFrom stats runif
#' @export

msa_set_scoreSA <- function(ref, alt, 
                            exec = "/Users/krah/Documents/R/pkgs/polenta/src/msa_set_score_src/msa_set_score", 
                            bootstrap){

  # create temporary dir with temporary fasta files
  fns <- vector(length = 1)
  rn <- format(runif(1, 1, 100000), digits = 0, scientific = FALSE)
  fns <- tempfile(pattern = paste0("run", rn), tmpdir = tempdir(), fileext = ".fas")
  unlink(fns[file.exists(fns)])

  if (!is.character(alt)){
    cat("list was supplied: writing MSAs to temporary files")
    dir.create(paste(tempdir(), "alt", sep = "/"))
    dir <- paste(tempdir(), "alt", sep = "/")
    msa_out <- vector(length = bootstrap)
    for (i in seq_along(msa_out))
      msa_out[i] <- tempfile(pattern = "mafft", tmpdir = dir, fileext = ".fas")
    unlink(msa_out[file.exists(msa_out)])

    ## Undefined global functions or variables: msa [CH-2017-11-08]
    # for (i in 1:bootstrap)
    #   write.fas(msa[[i]], file = msa_out[i])

    write.fas(ref, fns[1])
    system(paste(exec, fns[1], paste0(tempdir(), "/alt"), "-d", dir),
           intern = TRUE, ignore.stdout = FALSE)
  }

  ## Where is the executable? This will be obsolete when our own
  ## C code will be called.
  # exec <- "/Users/heibl/Documents/r/pkgs/polenta/src/msa_set_score_src/msa_set_score"
  else {
    if (!dir.exists(alt))
      stop("directory not found")

    write.fas(ref, fns[1])
    dir <- tempdir()
    system(paste(exec, fns[1], alt, "-d", dir),
           intern = TRUE, ignore.stdout = FALSE)
  }

  ## read program putput which is in temp dir
  files <-  list.files(tempdir(), full.names = TRUE)
  read <- files[grep("alt", files)]

  ## Scores
  # Column score
  CS <- read.table(read[1])
  names(CS) <- c("col", "CS")

  # Mean scores
  msc <- readLines(read[2])
  msc <- gsub("#", "", msc[5])
  msc <- strsplit(msc, "  ")
  msc <- unlist(msc)
  msc <- data.frame(do.call(cbind, strsplit(msc, " ")))
  msc[] <- lapply(msc, as.character)
  names(msc) <- c(msc[1, 1], msc[1, 2])
  msc <- msc[-1, ]


  # Residue pair column score (GUIDANCE Score)
  g.sc <- read.table(read[3])
  names(g.sc) <- c("col", "col_score")

  # Residue pair residue score
  rpr.sc <- read.table(read[4])
  names(rpr.sc) <- c("col", "residue", "score")

  # Residual pair sequence pair score
  rpsp.sc <- read.table(read[5])
  names(rpsp.sc) <- c("seq_row1", "seq_row2", "score")

  # Residual pair sequence score
  rps.sc <- read.table(read[6])
  names(rps.sc) <- c("seq", "score")

  # Residue pair score

  rp.sc <- read.table(read[7])
  names(rp.sc) <- c("col1", "row1", "row2", "score")

  ## delete temp files
  unlink(fns[file.exists(fns)], force = TRUE)
  unlink(fns, force = TRUE)
  files <- list.files(tempdir(), full.names = TRUE)
  files <- files[-grep("rs-graphics", files)]
  unlink(files, force = TRUE, recursive = TRUE)


  # output
  res <- list(score_means = msc,
              column_score = CS,
              residue_pair_column_score = g.sc,
              residue_pair_residue_score = rpr.sc,
              residual_pair_sequence_pair_score  = rpsp.sc,
              residual_pair_sequence_score = rps.sc,
              residue_pair_score = rp.sc)

  return(res)
}

