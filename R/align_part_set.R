#' @title Co-optimal alignments from tree partitions
#' @description  Produce co-optimal MSAs from the N-3 bipartitions of the
#'   starting tree. This function performs one iteration i 1: Take i partition
#'   and align both sides Heads (H) and Tails (T) and we get H1,H2 and T1, T2 2:
#'   align combinations again Heads and Tails: H1H2, H1T1, H1T2, T1T2 => 8
#'   combinations
#' @param x An object of class \code{\link{DNAbin}} or
#'   \code{\link{AAbin}} containing unaligned sequences of DNA or amino acids. Here inherited from function \code{\link{guidance2}} or \code{\link{HoT}}.
#' @param partition_set partitions of the alignment based on the guide tree of the base MSA. Here inherited from the function \code{\link{partitions}}.
#' @param msa.exec A character string giving the path to the executable of the
#'   alignment program (e.g. \code{/usr/local/bin/mafft}); possible programs are
#'   \code{MAFFT}, \code{MUSCLE}, and \code{ClustalW}.
#' @param method A character string containing further arguments passed to
#'   MAFFT; default is \code{"auto"}.
#' @param coopt.sub option for selection of specific co-optimal alignments, currently not available and set to "all"
#' @import ips ape
#' @importFrom ips mafft
#' @author Franz-Sebastian Krah
#' @noRd

align_part_set <- function(x, partition_set, msa.exec,
                           method, coopt.sub = "all"){

  
  ## Look up MSA program specified
  msa.program <- str_extract(msa.exec, "mafft|muscle|clustal\\w")
  
  if (coopt.sub == "all"){
    coopt.sub <- 1:8
  }

  seq_left <- x[partition_set > 0]
  seq_right <- x[partition_set == 0]

  ## MAFFT
  ###########
  if (msa.program == "mafft"){
    ## Aligning left HEADS and TAILS
    if (length(seq_left) == 1){
      headsA <- seq_left
      tailsA <- seq_left
    } else {
      headsA <- mafft(seq_left, exec = msa.exec, method = method)
      tailsA <- mafft(rev_DNA(seq_left), exec = msa.exec, method = method)
      tailsA <- rev_DNA(tailsA)
    }
    ## Aligning right HEADS and TAILS
    if (length(seq_right) == 1){
      headsB <- seq_right
      tailsB <- seq_right
    } else {
      headsB <- mafft(seq_right, exec = msa.exec, method = method)
      tailsB <- mafft(rev_DNA(seq_right), exec = msa.exec, method = method)
      tailsB <- rev_DNA(tailsB)
    }
    # aling 4 combinations of the basic MSAs (heads) and also
    # align in revers direction (tails)
    # headsA - headsB - HEADS
    if (1 %in% coopt.sub){
      msa1 <- mafft(x = headsA, y = headsB, add = "add",
                    method = method, exec = msa.exec)
    }
    # headsA - headsB - TAILS
    if (2 %in% coopt.sub){
      msa2 <- mafft(x = rev_DNA(headsA), y = rev_DNA(headsB), add = "add",
                    method = method, exec = msa.exec)
      msa2 <- rev_DNA(msa2)
    }
    # headsA - tailsB - HEADS
    if (3 %in% coopt.sub){
      msa3 <- mafft(headsA, tailsB, add = "add",
                    method = method, exec = msa.exec)
    }
    # headsA - tailsB - TAILS
    if (4 %in% coopt.sub){
      msa4 <- mafft(rev_DNA(headsA), rev_DNA(tailsB), add = "add",
                    method = method, exec = msa.exec)
      msa4 <- rev_DNA(msa4)
    }
    # tailsA - headsB - HEADS
    if (5 %in% coopt.sub){
      msa5 <- mafft(tailsA, headsB, add = "add",
                    method = method, exec = msa.exec)
    }
    # tailsA - headsB - TAILS
    if (6 %in% coopt.sub){
      msa6 <- mafft(rev_DNA(tailsA), rev_DNA(headsB), add = "add",
                    method= method, exec = msa.exec)
      msa6 <- rev_DNA(msa6)
    }
    # tailsA - tailsB - HEADS
    if (7 %in% coopt.sub){
      msa7 <- mafft(rev_DNA(tailsA), rev_DNA(tailsB), add = "add",
                    method = method, exec = msa.exec)
    }
    if (8 %in% coopt.sub){
      msa8 <- mafft(tailsA, tailsB, add = "add",
                    method = method, exec = msa.exec)
      msa8 <- rev_DNA(msa8)
    }
    ## list of 8 combs
    list.msas <- paste(paste0("msa", coopt.sub), collapse = ",")
    comb_msas <- eval(parse(text = paste0("list(", list.msas, ")")))
  }

  ## MUSCLE
  ###########
  if (msa.program == "muscle"){
    if (length(seq_left) == 1){
      headsA <-  seq_left
      tailsA <- rev_DNA(seq_left)
    } else {
      headsA <- muscle(seq_left, exec = msa.exec)
      tailsA <- muscle(rev_DNA(seq_left), exec = msa.exec)
    }
    if (length(seq_right) == 1){
      headsB <-  seq_right
      tailsB <- rev_DNA(seq_right)
    } else {
      headsB <- muscle(seq_right, exec = msa.exec)
      tailsB <- muscle(rev_DNA(seq_right), exec = msa.exec)
    }
    if (1 %in% coopt.sub)
      msa1 <- muscle(x = headsA, y = headsB, exec = msa.exec)

    if (2 %in% coopt.sub){
      msa2 <- muscle(x = rev_DNA(headsA), y = rev_DNA(headsB),
                      exec = msa.exec)
      msa2 <- rev_DNA(msa2)
    }
    if (3 %in% coopt.sub)
      msa3 <- muscle(headsA, rev_DNA(tailsB),  exec = msa.exec)

    if (4 %in% coopt.sub){
      msa4 <- muscle(rev_DNA(headsA), tailsB,  exec = msa.exec)
      msa4 <- rev_DNA(msa4)
    }
    if (5 %in% coopt.sub)
      msa5 <- muscle(rev_DNA(tailsA), headsB, exec = msa.exec)

    if (6 %in% coopt.sub){
      msa6 <- muscle(tailsA, rev_DNA(headsB), exec = msa.exec)
      msa6 <- rev_DNA(msa6)
    }
    if (7 %in% coopt.sub)
      msa7 <- muscle(rev_DNA(tailsA), rev_DNA(tailsB), exec = msa.exec)
    if (8 %in% coopt.sub){
      msa8 <- muscle(tailsA, tailsB, exec = msa.exec)
      msa8 <- rev_DNA(msa8)
    }
    list.msas <- paste(paste0("msa", coopt.sub), collapse =",")
    comb_msas <- eval(parse(text = paste0("list(", list.msas, ")")))
  }


  ## CLUSTAL
  ###########
  if (msa.program == "clustalo"){
    if (length(seq_left) == 1){
      headsA <-  seq_left
      tailsA <- rev_DNA(seq_left)
    } else {
      headsA <- clustalomega(seq_left, exec = msa.exec)
      tailsA <- clustalomega(rev_DNA(seq_left), exec = msa.exec)
    }
    if (length(seq_right) == 1){
      headsB <-  seq_right
      tailsB <- rev_DNA(seq_right)
    } else {
      headsB <- clustalomega(seq_right,  exec = msa.exec)
      tailsB <- clustalomega(rev_DNA(seq_right), exec = msa.exec)
    }
    if (1 %in% coopt.sub)
      msa1 <- clustalomega(x = headsA, y = headsB, exec = msa.exec)
    if (2 %in% coopt.sub){
      msa2 <- clustalomega(x = rev_DNA(headsA), y = rev_DNA(headsB),
                       exec = msa.exec)
      msa2 <- rev_DNA(msa2)
    }
    if (3 %in% coopt.sub)
      msa3 <- clustalomega(headsA, rev_DNA(tailsB), exec = msa.exec)
    if (4 %in% coopt.sub){
      msa4 <- clustalomega(rev_DNA(headsA), tailsB, exec = msa.exec)
      msa4 <- rev_DNA(msa4)
    }
    if (5 %in% coopt.sub)
      msa5 <- clustalomega(rev_DNA(tailsA), headsB, exec = msa.exec)
    if (6 %in% coopt.sub){
      msa6 <- clustalomega(tailsA, rev_DNA(headsB), exec = msa.exec)
      msa6 <- rev_DNA(msa6)
    }
    if (7 %in% coopt.sub)
      msa7 <- clustalomega(rev_DNA(tailsA), rev_DNA(tailsB), exec = msa.exec)
    if (8 %in% coopt.sub){
      msa8 <- clustalomega(tailsA, tailsB,  exec = msa.exec)
      msa8 <- rev_DNA(msa8)
    }
    list.msas <- paste(paste0("msa", coopt.sub), collapse = ",")
    comb_msas <- eval(parse(text = paste0("list(", list.msas, ")")))
  }

  ## CLUSTALW2
  #############
  if (msa.program == "clustalw"){
    if (length(seq_left) == 1){
      headsA <-  seq_left
      tailsA <- rev_DNA(seq_left)
    } else {
      headsA <- clustal(seq_left, exec = msa.exec)
      tailsA <- clustal(rev_DNA(seq_left), exec = msa.exec)
    }
    if (length(seq_right) == 1){
      headsB <-  seq_right
      tailsB <- rev_DNA(seq_right)
    } else {
      headsB <- clustal(seq_right,  exec = msa.exec)
      tailsB <- clustal(rev_DNA(seq_right), exec = msa.exec)
    }
    if (1 %in% coopt.sub)
      msa1 <- clustal(x = headsA, y = headsB, exec = msa.exec)
    if (2 %in% coopt.sub){
      msa2 <- clustal(x = rev_DNA(headsA), y = rev_DNA(headsB),
                        exec = msa.exec)
      msa2 <- rev_DNA(msa2)
    }
    if (3 %in% coopt.sub){
      msa3 <- clustal(headsA, rev_DNA(tailsB), exec = msa.exec)
    }
    if (4 %in% coopt.sub){
      msa4 <- clustal(rev_DNA(headsA), tailsB, exec = msa.exec)
      msa4 <- rev_DNA(msa4)
    }
    if (5 %in% coopt.sub)
      msa5 <- clustal(rev_DNA(tailsA), headsB, exec = msa.exec)
    if (6 %in% coopt.sub){
      msa6 <- clustal(tailsA, rev_DNA(headsB), exec = msa.exec)
      msa6 <- rev_DNA(msa6)
    }
    if (7 %in% coopt.sub)
      msa7 <- clustal(rev_DNA(tailsA), rev_DNA(tailsB), exec = msa.exec)
    if (8 %in% coopt.sub){
      msa8 <- clustal(tailsA, tailsB,  exec = msa.exec)
      msa8 <- rev_DNA(msa8)
    }
    list.msas <- paste(paste("msa", coopt.sub, sep=""), collapse = ",")
    comb_msas <- eval(parse(text = paste0("list(", list.msas, ")")))
  }
  return(comb_msas)
}
