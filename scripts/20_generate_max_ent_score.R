#' Title
#' Calculates the MaxEntScan scores (http://hollywood.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html) for the 9 bp motif sequence around the 
#' 5'ss (i.e. donor splice site) of a given intron and the 23 bp motif sequence around the 3'ss (i.e. acceptor splice site) of a given intron.
#' @param db.introns 
#' @param max.ent.tool.path Local path to the MaxEntScan tool
#' @param hs.fasta.path Local path to the homo sapiens fasta file
#'
#' @return
#' @export
#'
#' @examples
GenerateMaxEntScore <- function(db.introns,
                                max.ent.tool.path,
                                bedtools.path,
                                hs.fasta.path,
                                tmp.dir){

  logger::log_info("MaxEntScan score - extracting the sequences ...")  

  if (str_detect(db.introns$seqnames[1], "chr")) { 
    db.introns <- db.introns %>% mutate(seqnames = gsub("^chr", "", chr)) 
  } # Ensembl style manually 
  
  if (any(db.introns$seqnames == "M")) { logger::log_info("Error! There's data for chr-MT!") } 
  
  #========================================================================
  # EXTRACT THE GENOMIC SEQUENCES FOR EACH SPLICE SITE
  #========================================================================

  ## 0. Prepare the object ---------------------------------------------------
  
  setDT(db.introns)
  
  db.introns[strand == "+", `:=`(
    donorSeqStart = start - 4, donorSeqStop = start + 5,
    AcceptorSeqStart = end - 20, AcceptorSeqStop = end + 3
  )]
  db.introns[strand == "-", `:=`(
    donorSeqStart = end - 6, donorSeqStop = end + 3,
    AcceptorSeqStart = start - 4, AcceptorSeqStop = start + 19
  )]

  ## Remove rows with unresolved coordinates (e.g. strand == "*") before writing BED
  n_before <- nrow(db.introns)
  db.introns <- db.introns[!is.na(donorSeqStart) & !is.na(donorSeqStop) &
                             !is.na(AcceptorSeqStart) & !is.na(AcceptorSeqStop)]
  
  n_removed <- n_before - nrow(db.introns)
  if (n_removed > 0) {
    logger::log_info(paste0(n_removed, " row(s) removed due to missing/unresolved strand or coordinates."))
  }
  
  db.introns[1,]
  
  # Fix coordinates at source - cast to integer in data.table BEFORE building to.BED
  db.introns[, donorSeqStart := as.integer(donorSeqStart)]
  db.introns[, donorSeqStop  := as.integer(donorSeqStop)]
  db.introns[, AcceptorSeqStart := as.integer(AcceptorSeqStart)]
  db.introns[, AcceptorSeqStop  := as.integer(AcceptorSeqStop)]

  options(scipen = 999)

  ## 1. Get the genomic sequences for the DONOR splice sites ---------------------------------------------------

  to.BED <- data.frame(
    seqnames = db.introns$seqnames,
    starts   = as.integer(db.introns$donorSeqStart),
    ends     = as.integer(db.introns$donorSeqStop),
    names    = as.character(db.introns$junID),
    scores   = rep(".", nrow(db.introns)),
    strands  = db.introns$strand
  )

  tmp.file <- tempfile()
  tmp.file_seq <- tempfile()

  # Use write.table - slow to write but only called once
  write.table(to.BED, file = tmp.file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

  system.time({system(command = paste0(
    bedtools.path, "/bin/bedtools getfasta -name -s -fi ", hs.fasta.path,
    " -bed ", tmp.file, " -tab -fo ", tmp.file_seq
  ))})

  donor_sequences_input <- data.table::fread(tmp.file_seq, header = FALSE) %>% 
    mutate(junID = gsub("::.*$", "", as.character(V1)))

  # Sanity check
  stopifnot(identical(as.character(donor_sequences_input$junID), as.character(db.introns$junID)))

  db.introns <- cbind(db.introns, donor_sequence = as.character(donor_sequences_input$V2))
  db.introns %>% head()
  
  ## Remove temporary files
  rm(to.BED, tmp.file, tmp.file_seq)
  

  ## 2. Get the genomic sequences for the ACCEPTOR splice sites ---------------------------------------------------

  to.BED <- data.frame(seqnames  =  db.introns$seqnames,
                       starts    =  as.integer(db.introns$AcceptorSeqStart),
                       ends      =  as.integer(db.introns$AcceptorSeqStop),
                       names     =  as.character(db.introns$junID),
                       scores    =  c(rep(".", nrow(db.introns))),
                       strands   =  db.introns$strand)
  
  
  tmp.file <- tempfile() 
  tmp.file_seq <- tempfile()
  write.table(to.BED, file = tmp.file, quote = F, sep = "\t", row.names = F, col.names = F)
  system.time({system(command = paste0(
    bedtools.path, "/bin/bedtools getfasta -name -s -fi ", hs.fasta.path,
    " -bed ", tmp.file, " -tab -fo ", tmp.file_seq
  ))})
  acceptor_sequences_input <- data.table::fread(tmp.file_seq, header = FALSE) %>% mutate(junID = gsub("::.*$", "", as.character(V1)))
  
  head(acceptor_sequences_input)
  head(donor_sequences_input)
  
  ## Replaces everything (.*) from the '::' until the end of the string '$'
  # Sanity check
  stopifnot(identical(as.character(acceptor_sequences_input$junID), as.character(db.introns$junID)))
  db.introns <- cbind(db.introns, acceptor_sequence = as.character(acceptor_sequences_input$V2))
  db.introns %>% head()
  
  
  ## Remove temporary files
  rm(to.BED, tmp.file, tmp.file_seq)
  

  #========================================================================
  # GENERATE THE MAXENTSCORE
  #========================================================================
  
  ## 2. Generate the MaxEntScore --------------------------------------------------------------------
  logger::log_info("Generating MaxEntScan score for the donor sequences...")
  
  ## get the sequences
  tmp.file <- tempfile()
  ## get the maxentscan for the 5' splice site
  
  ## check how many sequences contain "N"
  length(grep("N",as.character(db.introns$donor_sequence)))
  
  write.table(gsub("N","A",as.character(db.introns$donor_sequence)),file=tmp.file,row.names=F,col.names=F,quote=F)
  setwd(max.ent.tool.path)
  ss5score <- read.delim(pipe(paste0("perl ", max.ent.tool.path, "score5.pl ", tmp.file)), header = F)
  identical(as.character(ss5score$V1), gsub("N","A",as.character(db.introns$donor_sequence)))
  db.introns <- cbind(db.introns, ss5score = ss5score$V2)


  logger::log_info("Generating MaxEntScan score for the acceptor sequences...")
  ## get the maxentscan for the 3' splice site
  length(grep("N",as.character(db.introns$acceptor_sequence)))
  
  write.table(gsub("N","A",as.character(acceptor_sequences_input$V2)),file=tmp.file,row.names=F,col.names=F,quote=F)
  ss3score <- read.delim(pipe(paste0("perl ", max.ent.tool.path, "/score3.pl ", tmp.file)),header = F)
  identical(as.character(ss3score$V1),gsub("N","A",as.character(db.introns$acceptor_sequence)))
  db.introns <- cbind(db.introns, ss3score = ss3score$V2)
    
  rm(ss5score, ss3score, tmp.file)
    
  return(db.introns)
}
