#' @title MSA Reliability Assessment with GUIDANCE
#' @name guidance
#' @description Calculate MSA reliability scores with GUIDANCE (Penn et al. 2010).
#' @param sequences An object of class \code{\link{DNAbin}} or \code{\link{AAbin}} containing unaligned sequences of DNA or amino acids.
#' @param bootstrap An integer giving the number of alternative MSAs to be computed.
#' @param method A character string containing further arguments passed to MAFFT; default is \code{"auto"}.
#' @param msa.exec A character string giving the path to the executable of the
#'   alignment program (e.g. \code{/usr/local/bin/mafft}); possible programs are
#'   \code{MAFFT}, \code{MUSCLE}, \code{ClustalW}, \code{ClustalO} or \code{PRANK}
#'   For details see \code{\link{clustal}}, \code{\link{mafft}}, \code{\link{prank}}
#' @param ncore An integer specifying the number of cores; default = 1 (i.e. serial execution); \code{"auto"} can be used for automated usage of all detected cores.
#' @param zip.file A character string giving the dir of zip-compressed file to be produced, which contains the  alternative MSAs. If left empty (default), the alternative MSA will not be stored and cannot be assessed by the user.
#' @return An object of class \code{\linkS4class{guidanceDNA}} or \code{\linkS4class{guidanceAA}}.
#' @details Calculates column confidence (and other scores) by comparing
#'   alternative MSAs generated by alternative guide trees derived from
#'   bootstrapped MSAs (Felsenstein 1985). The basic comparison between the bootstrapped MSAs
#'   and a reference MSA is if column residue pairs are identically aligned in
#'   all alternative MSAs compared with the base MSA (see \code{msa_set_score}).
#' @details For an example workflow see Vignette
#' @references Felsenstein, J. 1985. Confidence limits on phylogenies: an
#'   approach using the bootstrap. \emph{Evolution} \strong{39}:783--791.
#' @references Penn et al. 2010. An
#'   alignment confidence score capturing robustness to guide tree uncertainty.
#'   \emph{Molecular Biology and Evolution} \strong{27}:1759--1767.
#' @seealso \code{\link{msa_set_score}}, \code{\link{guidance2}}, \code{\link{HoT}}
#' @import ips
#' @importFrom doSNOW registerDoSNOW
#' @import foreach
#' @importFrom parallel detectCores makeCluster
#' @import pbmcapply
#' @import plyr
#' @importFrom phangorn as.phyDat dist.ml

#' @examples
#' \dontrun{
#' # run GUIDANCE on example data using MAFFT
#' fpath <- system.file("extdata", "BB30015.fasta", package="rGUIDANCE") # random example from BALiBASE
#' fas <- ape::read.FASTA(fpath)
#' g <- guidance(sequences = fas, msa.exec= "/usr/local/bin/mafft")
#' scores <- scores(g, score = "column")
#' plot(scores$column$score, xlab = "Site", 
#' ylab = "Column score", 
#' main = "GUIDANCE", type = "l")
#' }
#'
#' @author Franz-Sebastian Krah
#' 
#' @export

guidance <- function(sequences,
                     bootstrap = 100,
                     method = "auto",
                     msa.exec = "/usr/local/bin/mafft",
                     ncore = 1,
                     zip.file){
  
  ##############################################
  ## SOME CHECKS
  ##############################################
  if (!inherits(sequences, c("DNAbin", "AAbin")))
    stop("sequences not of class DNAbin or AAbin (ape)")
  
  nseq <- ifelse(is.matrix(sequences), nrow(sequences), length(sequences))
  
  if (nseq < 8)
    warning("GUIDANCE is not suitable for alignments of very few sequences.\n
            As a rule of thumb, use guidance2 or HoT for < 8 sequences")
  
  if (nseq > 199)
    warning("alignments with more than 200 sequences may run into computional problems")
  
  ## look up MSA program specified
  msa.program <- str_extract(msa.exec, "mafft|muscle|clustalo|clustalw|prank")
  if (!msa.program %in% c("mafft", "muscle", "clustalw", "lustalo", "prank"))
    stop("Currently only MAFFT, MUSCLE, ClustalW or ClustalO")
  
  ## Check for MSA program
  out <- system(paste(msa.exec, "--v"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (out == 127)
    stop("please provide msa.exec path or install MSA program in root \n
         i.e. in Unix: '/usr/local/bin/mafft'")
  
  ## generate some parameters if not specified
  ## -----------------------------------------
  ## number of cores
  if (ncore == "auto"){
    ncore <- detectCores(all.tests = FALSE, logical = TRUE)
  }
  
  ##############################################
  ## PART I
  ##############################################
  ## BASE and ALTERNATIVE MSAs
  ##############################################
  
  ## Make base alignment
  ## -------------------
  cat("Generating the base alignment \n")

  if (msa.program == "mafft") {
    
    base_msa <- mafft(x = sequences,
                      exec = msa.exec, method = method,
                      thread = -1)
  }
  
  if (msa.program == "clustalo") {
    
    base_msa <- clustalomega(x = sequences, 
                             exec = msa.exec, 
                             MoreArgs = "")
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
  
  if (msa.program == "prank") {
    
    base_msa <- prank(x = sequences,
                      path = msa.exec, 
                      gaprate = 0.025, 
                      gapext = 0.75)
  }
  
  ## Compute NJ guide trees
  ## ----------------------
  cat("Generating NJ guide trees\n")
  pb <- txtProgressBar(max = bootstrap, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  cl <- makeCluster(ncore)
  registerDoSNOW(cl)
  
  nj_guidetrees <- foreach(i = 1:bootstrap,
                           .options.snow = opts,
                           .packages = "phangorn", 
                           .export = 'msaBP_nj_tree') %dopar% {
                             msaBP_nj_tree(msa = base_msa, outgroup = "auto")
                           }
  stopCluster(cl)
  close(pb)
  
  ## Alignment of MSA BP times with new NJ guide trees
  ## -------------------------------------------------
  cat("Alignment of sequences using NJ guide trees\n")
  
  ## Construct alignment function
  ## ----------------------------

  ## loop
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  cl <- makeCluster(ncore)
  registerDoSNOW(cl)
  
  if (msa.program == "mafft") {
    
    alt_msa <- foreach(i = 1:bootstrap, 
                       .options.snow = opts) %dopar% {
                         
                         mafft(x = sequences, gt = nj_guidetrees[[i]], 
                               exec = msa.exec, method = method,
                               maxiterate = 0, op = 1.53, ep = 0,
                               thread = 1)
                       }
  }
  
  if (msa.program == "clustalo") {
  
    alt_msa <- foreach(i = 1:bootstrap, 
                       .options.snow = opts,
                       .export = c("sequences", "nj_guidetrees", "msa.exec")) %dopar% {
                         
                         clustalomega(x = sequences, guide.tree = nj_guidetrees[[i]], 
                                      exec = msa.exec, MoreArgs = "")
                         
                       }
  }
  
  if (msa.program == "clustalw") {
    
    alt_msa <- foreach(i = 1:bootstrap, 
                       .options.snow = opts,
                       .export = c("sequences", "nj_guidetrees", "msa.exec")) %dopar% {
                         
                         clustal(x = sequences, guide.tree = nj_guidetrees[[i]], 
                                 exec = msa.exec,
                                 pw.gapopen = 10, pw.gapext = 0.1,
                                 gapopen = 10, gapext = 0.2,
                                 MoreArgs = "")
                         
                       }
  }
  
  if (msa.program == "muscle") {

    alt_msa <- foreach(i = 1:bootstrap, 
                       .options.snow = opts,
                       .export = c("sequences", "nj_guidetrees", "msa.exec")) %dopar% {
                         
                         muscle(x = sequences, guide.tree = nj_guidetrees[[i]], 
                                exec = msa.exec, MoreArgs = "")
                       
                         }
  }
  
  if (msa.program == "prank") {
    
    alt_msa <- foreach(i = 1:bootstrap, 
                       .options.snow = opts,
                       .export = c("sequences", "nj_guidetrees", "msa.exec")) %dopar% {
                         
                         prank(x = sequences, guidetree = nj_guidetrees[[i]], 
                               path = msa.exec, 
                               gaprate = 0.025, 
                               gapext = 0.75)
                         
                       }
  }
  
  stopCluster(cl)
  close(pb)
  
  ## Delete leftover files from MAFFT
  file.remove(list.files(pattern = "tree.mafft", full.names = TRUE))
  
  ##############################################
  ## PART II
  ##############################################
  ## Computation of GUIDANCE scores
  ##############################################
  cat("\nCalculating GUIDANCE scores \n")
  
  ## Run msa_set_score
  score <- msa_set_score(ref = base_msa,
                          alt = alt_msa)
  
  ## Store alternative MSAs in a zip file (optional)
  ## -----------------------------------------------
  if (!missing(zip.file)){
    for (i in seq_along(alt_msa)){
      write.FASTA(alt_msa[[i]], file =paste0(zip.file, "alt_msa_", i, ".fas"))
    }
    
    files <- list.files(zip.file, full.names = TRUE)
    files <- files[grep("alt_msa", files)]
    zip(zipfile = paste0(zip.file, "guidance_alt_msas_", Sys.Date(), ".zip"), files = files)
    file.remove(files)
  }
  
  ## Delete temporary files
  # this deletion approach has proven best to delete everything
  files <- list.files(tempdir(), full.names = TRUE)
  files <- files[-grep("rs-graphics", files)]
  unlink(files, force = TRUE, recursive = TRUE)
  
  ## Return guidance class
  ## -------------------------
  if (inherits(sequences, "AAbin")){
    guidanceAA(base_msa, score, "guidance", msa.program)
  } else {
    guidanceDNA(base_msa, score, "guidance", msa.program)
  }
}
