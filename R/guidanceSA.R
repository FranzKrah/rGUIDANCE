#' @title Interface to guidance program
#'
#' @param seqs An object of class \code{\link{DNAbin}} or
#'   \code{\link{AAbin}} containing unaligned sequences of DNA or amino acids.
#' @param bootstrap An integer giving the number of perturbated MSAs.
#' @param msa.program A charcter string giving the name of the MSA program,
#'   currelty one of c("mafft", "muscle", "clustalo", "clustalw2"); MAFFT is
#'   default
#' @param program A charcter string giving the name one of c("guidance", "guidance2", "HoT").
#' @param outorder As input or aligned. Default=aligned
#' @param msafile Path to user supplied local MSA file
#' @param cutoff Confidence cutoff between 0 to 1. Default=0.93
#' @param moreArgs More arguments passed to GUIDANCE
#' @param exec path to guidance program folder, e.g.
#'   "/Applications/guidance.v2.02/"
#' @param proc_num Integer giving the number of cores.
#' @param quiet logical if TRUE, progress is printed to console.
#' @return list containing following scores and alignments:
#' @return mean_scores residue pair score and mean column score
#' @return column_score
#' @return residue_column_score GUIDANCE score
#' @return residue_pair_residue_score
#' @return residual_pair_sequence_pair_score
#' @return residual_pair_sequence_score
#' @return residue_pair_score
#' @return base_msa
#' @return guidance_msa is the base_MSA removed from unreliable
#'   residues/columns/sequences below cutoffs
#'
#' @import stringr
#' @import useful
#' @import foreach
#' @import ips
#' @importFrom utils read.table
#' @author Franz-Sebastian Krah

guidanceSA <- function(seqs, bootstrap, msa.program, program,
                       outorder, msafile, cutoff = 0.93,
                       moreArgs, exec, proc_num, quiet = FALSE){

  perl.call <- paste("perl", paste0(exec, "www/Guidance/guidance.pl"))
  t <- system(perl.call, ignore.stderr = TRUE)
  if (t == 2){
    stop("GUIDANCE is not responding!")
  }


  if (!inherits(seqs, c("DNAbin", "AAbin")))
    stop("sequences not of class DNAbin or AAbin (ape)")

  type <- class(seqs)
  type <- gsub("bin", "", type)
  if (type == "DNA") type <- "nuc"
  if (type == "AA") type <- "aa"

  fns <- vector(length = 1)
  for (i in seq_along(fns))
    fns[i] <- tempfile(pattern = "guidance", tmpdir = tempdir(), fileext = ".fas")
  unlink(fns[file.exists(fns)])

  write.fas(seqs, fns[1])

  # necessary
  seqFile <- paste("--seqFile", fns[1])
  if (lower.case(msa.program))
    msa.program <- str_to_upper(msa.program)
  msa.program <- paste("--msaProgram", msa.program)

  seqType <- paste("--seqType", type)
  # outDir <- paste("--outDir", outdir, sep=" ")
  if (program == "hot")
    program <- "HoT"
  if (lower.case(program))
    program <- str_to_upper(program)
  program <- paste("--program", program)
  bootstraps <- paste("--bootstraps", bootstrap)
  outdir <- paste(tempdir(), "guidance", sep = "/")
  outDir <- paste("--outDir", outdir)


  if (!missing(moreArgs))
    moreArgs <- paste("--MSA_Param", moreArgs)
  if (!missing(gencode))
    gencode <- paste("--genCode", gencode)
  if (!missing(outorder))
    outorder <- paste("--outOrder", outorder)
  if (!missing(msafile))
    msafile <- paste("--msaFile", msafile)
  if (!missing(cutoff))
    cutoff <- paste("--seqCutoff", cutoff)
  if (!missing(proc_num))
    proc_num <- paste("--proc_num", proc_num)

  
  moreArgs <- mget(c("moreArgs", "gencode","outorder",
    "msafile", "cutoff", "proc_num"), ifnotfound = "")
  moreArgs <- unlist(moreArgs)
  moreArgs <- moreArgs[nchar(moreArgs)>0]
  
  guidance.call <- paste(seqFile, msa.program, seqType,
                         outDir, program, bootstraps, moreArgs)
  guidance.call <- paste(perl.call, guidance.call)
  cat("\n", guidance.call, "\n")
  ## CALL
  if (quiet){
    system(guidance.call, ignore.stdout = TRUE, ignore.stderr = TRUE)
  } else {
    system(guidance.call)
  }

  files <- list.files(paste(tempdir(), "guidance", sep = "/"), full.names = TRUE)
  # files <- list.files("../../../PhD/proj/high_priority/color_new_alignment/rpb1.fguidance/",full.names = TRUE)
  read <- files[grep("\\.scr\\b", files)]
  read <- read[-grep("csv", read)]

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

  # base.msa
  base <- files[grep("MSA.MAFFT.aln.With_Names", files)]
  base.msa <- read.fas(base)

  # guidance.msa
  guidance.msa <- files[grep("Without_low_SP_Col.With_Names", files)]
  guidance.msa <- read.fas(guidance.msa)

  ## delete temp files
  unlink(fns[file.exists(fns)], recursive = TRUE)
  unlink(files, force = TRUE)
  files <- list.files(tempdir(), full.names = TRUE)
  files <- files[-grep("rs-graphics", files)]
  unlink(files, force = TRUE, recursive = TRUE)


  # output
  res <- list(scores = list(mean_score = msc,
                            column_score = CS,
                            residue_pair_column_score = g.sc,
                            residue_pair_residue_score = rpr.sc,
                            residual_pair_sequence_pair_score  = rpsp.sc,
                            residual_pair_sequence_score = rps.sc,
                            residue_pair_score = rp.sc),
              reduced_msa = guidance.msa,
              base_msa = base.msa)
  return(res)
}
