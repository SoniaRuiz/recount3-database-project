#' Title
#' Calculates the MaxEntScan scores (http://hollywood.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html) for the 9 bp motif sequence around the 
#' 5'ss (i.e. donor splice site) of a given intron and the 23 bp motif sequence around the 3'ss (i.e. acceptor splice site) of a given intron.
#' @param junc_tidy 
#' @param max_ent_tool_path Local path to the MaxEntScan tool
#' @param homo_sapiens_fasta_path Local path to the homo sapiens fasta file
#'
#' @return
#' @export
#'
#' @examples
GenerateMaxEntScore <- function(junc_tidy,
                                max_ent_tool_path,
                                homo_sapiens_fasta_path){
  
  junc_tidy <- junc_tidy %>% dplyr::as_tibble()
  
  junc_tidy$seqnames <- junc_tidy$seqnames %>% as.character()
  
  if (any(junc_tidy$seqnames == "M")) {
    logger::log_info("Error! There's data for chr-MT!")
  } 
  
  ## 0. Prepare the object ---------------------------------------------------
  
  ## get the ranges for the donor and acceptor sequences needed for the MaxEntScan
  junc_tidy <- junc_tidy %>%  mutate(donorSeqStart = 
                                       ifelse(strand == "-",
                                              end - 6, start - 4),
                                     donorSeqStop =
                                       ifelse(strand == "-",
                                              end + 3, start + 5),
                                     AcceptorSeqStart =
                                       ifelse(strand == "-",
                                              start - 4, end - 20),
                                     AcceptorSeqStop =
                                       ifelse(strand == "-",
                                              start + 19, end + 3)) 
  
  junc_tidy[1,]
  
  to.BED <- data.frame(seqnames = junc_tidy$seqnames,
                       starts   = as.integer(junc_tidy$donorSeqStart),
                       ends     = as.integer(junc_tidy$donorSeqStop),
                       names    = as.character(junc_tidy$junID),
                       scores   = c(rep(".", nrow(junc_tidy))),
                       strands  = junc_tidy$strand)
  to.BED[1,]
  
  ## 1. Obtain the genomic sequence for splice sites ------------------------------------------------
  
  ## Get the donor genomic sequence
  
  tmp.file <- tempfile()
  
  ## get the maxentscan for the 5' splice site
  write.table(to.BED, file = tmp.file, quote = F, sep = "\t", row.names = F, col.names = F)
  tmp.file_seq <- tempfile()
  system(paste0("bedtools getfasta -name -s -fi ", homo_sapiens_fasta_path, " -bed ",
                tmp.file, " -tab -fo ", tmp.file_seq))
  donor_sequences_input <- read.delim(tmp.file_seq, header = F)
  head(donor_sequences_input)
  head(junc_tidy)
  
  
  stopifnot(identical(gsub("\\(\\+\\)", "", gsub("\\(\\*\\)", "", gsub("\\(-\\)", "", as.character(donor_sequences_input$V1)))),
                      junc_tidy$junID %>% as.character()))
  junc_tidy <- cbind(junc_tidy, 
                     donor_sequence = as.character(donor_sequences_input$V2))
  
  junc_tidy %>% head()
  
  
  ## Get the acceptor genomic sequence
  
  to.BED <- data.frame(seqnames  =  junc_tidy$seqnames,
                       starts    =  as.integer(junc_tidy$AcceptorSeqStart),
                       ends      =  as.integer(junc_tidy$AcceptorSeqStop),
                       names     =  as.character(junc_tidy$junID),
                       scores    =  c(rep(".", nrow(junc_tidy))),
                       strands   =  junc_tidy$strand)
  
  
  tmp.file <- tempfile()
  
  write.table(to.BED, file = tmp.file, quote = F, sep = "\t", row.names = F, col.names = F)
  tmp.file_seq <- tempfile()
  system(paste0("bedtools getfasta -name -s -fi ", homo_sapiens_fasta_path, " -bed ",
                tmp.file, " -tab -fo ", tmp.file_seq))
  acceptor_sequences_input <- read.delim(tmp.file_seq, header = F)
  
  head(acceptor_sequences_input)
  head(donor_sequences_input)
  
  stopifnot(identical(gsub("\\(\\+\\)", "", gsub("\\(\\*\\)", "", gsub("\\(-\\)", "", as.character(acceptor_sequences_input$V1)))),
                      junc_tidy$junID %>% as.character()))
  junc_tidy <- cbind(junc_tidy,
                     acceptor_sequence = as.character(acceptor_sequences_input$V2))
  
  junc_tidy %>% head()
  
  
  ## Remove temporary files
  rm(to.BED, tmp.file, tmp.file_seq)
  
  
  
  ## 2. Generate the MaxEntScore --------------------------------------------------------------------
  
  
  ## get the sequences
  tmp.file <- tempfile()
  ## get the maxentscan for the 5' splice site
  
  ## check how many sequences contain "N"
  length(grep("N",as.character(junc_tidy$donor_sequence)))
  
  write.table(gsub("N","A",as.character(junc_tidy$donor_sequence)),file=tmp.file,row.names=F,col.names=F,quote=F)
  setwd(max_ent_tool_path)
  ss5score <- read.delim(pipe(paste0("perl ", max_ent_tool_path, "score5.pl ", tmp.file)),header = F)
  identical(as.character(ss5score$V1),gsub("N","A",as.character(junc_tidy$donor_sequence)))
  junc_tidy <- cbind(junc_tidy, ss5score = ss5score$V2)
  
  logger::log_info("MaxEntScan score generated for the donor sequences!")
  
  
  ## get the maxentscan for the 3' splice site
  length(grep("N",as.character(junc_tidy$acceptor_sequence)))
  
  write.table(gsub("N","A",as.character(acceptor_sequences_input$V2)),file=tmp.file,row.names=F,col.names=F,quote=F)
  ss3score <- read.delim(pipe(paste0("perl ", max_ent_tool_path, "/score3.pl ", tmp.file)),header = F)
  identical(as.character(ss3score$V1),gsub("N","A",as.character(junc_tidy$acceptor_sequence)))
  junc_tidy <- cbind(junc_tidy, ss3score = ss3score$V2)
  
  logger::log_info("MaxEntScan score generated for the acceptor sequences!")
  
  rm(ss5score, ss3score, tmp.file)
  
  return(junc_tidy)
}