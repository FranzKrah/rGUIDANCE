#' @title Heads or Tails Alignment Reliability
#' @description MSA reliability assessment HoT (Landan and Graur 2008)
#' @param sequences object of class \code{\link{DNAbin}} or \code{\link{AAbin}}
#'   containing unaligned sequences of DNA or amino acids.
#' @param method further argument passed to MAFFT, default is \code{"auto"}
#' @param bootstrap integer giving the number of alternative MSAs to be computed
#' @param msa.exec A character string giving the path to the executable of the
#'   alignment program (e.g. \code{/usr/local/bin/mafft}); possible programs are
#'   \code{MAFFT}, \code{MUSCLE}, and \code{ClustalW}. 
#'   For details see \code{\link{clustal}}, \code{\link{mafft}}
#' @param ncore integer specifying the number of cores; default = 1 (serial),
#'   "auto" can be used for automated usage of all detected cores.
#' @param zip.file A character string giving the name for the output zip file.
#' @return An object of class \code{\linkS4class{guidanceDNA}} or
#'   \code{\linkS4class{guidanceAA}}.
#' @details Calculates column reliability (and other scors) by comparing
#'   alternative MSAs generated by aligning guide tree partitions as described
#'   in Landan and Graur (2008). For details see \code{compareMSAs}. 8*(N-3)
#'   alternative MSAs are generated by default, where N is the number of
#'   sequences.
#' @details For an example workflow see Vignette
#' @references Landan and Graur. 2008. Local reliability measures from
#'   sets of co-optimal multiple sequence alignments. \emph{Pacific Symposium on
#'   Biocomputing} \strong{13}:15--24.
#' @seealso \code{\link{msa_set_score}}, \code{\link{guidance}},
#'   \code{\link{guidance2}}
#' @examples
#' \dontrun{
#' # run GUIDANCE on example data using MAFFT
#' fpath <- system.file("extdata", "BB30015.fasta", package="rGUIDANCE") # random example from BALiBASE
#' fas <- ape::read.FASTA(fpath)
#' g <- HoT(sequences = fas, msa.exec= "/usr/local/bin/mafft")
#' scores <- scores(g, score = "column")
#' plot(scores$column$score, xlab = "Site", 
#' ylab = "Column score", 
#' main = "HoT", type ="l")
#' }
#'
#' @author Franz-Sebastian Krah
#' @importFrom ape compute.brlen ladderize multi2di Ntip
#' @importFrom ips mafft read.fas
#' @import doSNOW
#' @import foreach
#' @importFrom graphics legend
#' @importFrom stringr str_extract
#' @import parallel
#' @import pbmcapply
#' @import plyr
#' @importFrom phangorn as.phyDat dist.ml
#' @importFrom utils setTxtProgressBar txtProgressBar zip
#' @export


HoT <- function(sequences, method = "auto",
                bootstrap,
                msa.exec = "/usr/local/bin/mafft",
                ncore = 1,
                zip.file) {
  
  ##############################################
  ## SOME CHECKS
  ##############################################
  
  if (!inherits(sequences, c("DNAbin","AAbin")))
    stop("sequences not of classes DNAbin or AAbin (ape)")
  
  nseq <- ifelse(is.matrix(sequences), nrow(sequences), length(sequences))
  if (nseq > 199)
    warning("alignments with more than 200 sequences may run into computional problems")
  
  ## look up MSA program specified
  msa.program <- str_extract(msa.exec, "mafft|muscle|clustalo|clustalw|prank")
  if (!msa.program %in% c("mafft", "muscle", "clustalw"))
    stop("Currently only MAFFT, MUSCLE or ClustalW")
  
  ## Check for MSA program
  out <- system(paste(msa.exec, "--v"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (out == 127)
    stop("please provide msa.exec path or install MSA program in root \n
         i.e. in Unix: '/usr/local/bin/mafft'")
  
  ## Check for MSA program
  ## ---------------------
  out <- system(paste(msa.exec, "--v"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (out == 127)
    stop("please provide exec path or install MSA program in root \n
         i.e. in Unix: '/usr/local/bin/mafft'")
  
  ## generate some parameters if not specified
  #---------------------------------------------
  ## number of cores
  if (ncore == "auto") {
    ncore <- detectCores(all.tests = FALSE, logical = TRUE)
  }
  
  ##############################################
  ## PART I
  ##############################################
  ## BASE and ALTERNATIVE MSAs
  ##############################################
  cat("Generating the base alignment \n")
  
  ## create loop input
  if (msa.program == "mafft") {
    
    base_msa <- mafft(x = sequences,
                      exec = msa.exec, method = method,
                      maxiterate = 0, op = 1.53, ep = 0,
                      thread = -1)
  }
  
  if (msa.program == "clustalw") {
    
    base_msa <- clustal(x = sequences,
                        exec = msa.exec,
                        pw.gapopen = 10, pw.gapext = 0.1,
                        gapopen = 10, gapext = 0.2,
                        MoreArgs = "")
  }
  
  if (msa.program == "muscle") {
    
    base_msa <- muscle(x = sequences,
                       exec = msa.exec, 
                       MoreArgs = "")
  }
  
  
  ## Calculate start guide tree
  #----------------------------------------------
  cat("Calculate start tree \n")
  base.msa.ml <- as.phyDat(base_msa)
  # find ML distance as input to nj tree search
  ml.dist.msa <- dist.ml(base.msa.ml)
  # NJ
  start_tree <- nj(ml.dist.msa)
  start_tree <- multi2di(start_tree)
  start_tree <- compute.brlen(start_tree)
  
  ## produce MSA partitions
  align_parts <- partitions(start_tree)
  
  # here could be a sampling of co-opts like in
  # guidance2. now we sample all
  n.coopt.sub <- rep("all", ncol(align_parts))
  n.coopt <- (Ntip(start_tree) - 3) * 8
  
  
  ##############################################
  ## PART II
  ##############################################
  ## Co-optimal MSAs
  ##############################################
  cat(paste("Sampling", n.coopt, "co-optimal alignments \n", sep = " "))
  
  ## Create temporary files
  #----------------------------------------------
  msa_out <- vector(length = n.coopt)
  for (i in seq_along(msa_out))
    msa_out[i] <- tempfile(pattern = "HoT", tmpdir = tempdir(), fileext = ".fas")
  unlink(msa_out[file.exists(msa_out)])
  
  # predifined file storage allocation (because it runs in batches of 8)
  start <- seq(1, n.coopt, 8)
  end <- seq(8, n.coopt, 8)
  stend <- data.frame(start, end)
  
  ## Run batch alignments
  #----------------------------------------------
  pb <- txtProgressBar(max = ncol(align_parts), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  cl <- makeCluster(ncore)
  registerDoSNOW(cl)
  
  
  alt_msa <- foreach(i = 1:ncol(align_parts),
                     .options.snow = opts,
                     .packages = "ape"
  ) %dopar% {
    align_part_set(x = sequences, partition_set = align_parts[,i], 
                   coopt.sub = n.coopt.sub[i],
                   method = method, msa.exec = msa.exec)
  }
  
  stopCluster(cl)
  close(pb)
  
  #### unlist nested list
  alt_msa <- foreach(i = 1:length(alt_msa), .combine = c) %do% {
    alt_msa[[1]]
  }
  ##############################################
  ## PART III
  ##############################################
  ## Computation of reliability scores
  ##############################################
  cat("Calculating reliability scores \n")
  
  ## HoT Score
  #----------------------------------------------
  score <- msa_set_score(ref = base_msa,
                          alt = alt_msa)
  
  ##  if wanted, store alternative MSAs into a zip file
  if (!missing(zip.file)){
    for(i in seq_along(alt_msa)){
      write.FASTA(alt_msa[[i]], file =paste0(zip.file, "alt_msa_", i, ".fas"))
    }
    
    files <- list.files(zip.file, full.names = T)
    files <- files[grep("alt_msa", files)]
    zip(zipfile = paste0(zip.file, "Hot_alt_msas_", Sys.Date(), ".zip"), files = files)
    file.remove(files)
  }
  
  ## Return guidance class
  ## -------------------------
  if (inherits(sequences, "AAbin")){
    guidanceAA(base_msa, score, "HoT", msa.program)
  } else {
    guidanceDNA(base_msa, score, "HoT", msa.program)
  }
}
