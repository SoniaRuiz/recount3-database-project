################################################################
## UTILS
## functions to help in the QC data, add intronic features
#################################################################


get_mode <- function(vector) {
  uniqv <- unique(vector)
  uniqv[which.max(tabulate(match(vector, uniqv)))]
}
 
remove_encode_blacklist_regions <- function(GRdata,
                                            blacklist_path) {
  
  
  if (!exists("encode_blacklist_hg38")) {
    encode_blacklist_hg38 <- rtracklayer::import(con = blacklist_path) %>% diffloop::rmchr()
  } else {
    print("'encode_blacklist_hg38' file already loaded!")
  }
  
  overlaped_junctions <- GenomicRanges::findOverlaps(query = encode_blacklist_hg38, 
                                                     subject = GRdata %>% diffloop::rmchr(),
                                                     #type = "any",
                                                     ignore.strand = F)
  
  ## JuncID indexes to be removed: they overlap with a black region
  indexes <- S4Vectors::subjectHits(overlaped_junctions)
  
  if (length(indexes) > 0) {
    print(paste0(length(unique(indexes)), " junctions overlap with a ENCODE blacklist region"))
    GRdata <- GRdata[-indexes, ]
    print(paste0(length(GRdata), " junctions after removing overlaps with ENCODE BlackList regions!"))
  }else{
    print("No junctions overlapping with an ENCODE blacklist region")
  }
  
  ## Return tidied data
  return(GRdata)
}

generate_max_ent_score <- function(junc_tidy,
                                   # bedtools_path = NULL,
                                   max_ent_tool_path,
                                   homo_sapiens_fasta_path){
  
  library(Biostrings)
  library(tidyverse)
  library(protr)
  
  junc_tidy <- junc_tidy %>% dplyr::as_tibble()
  junc_tidy$seqnames <- junc_tidy$seqnames %>% as.character()
  
  if (any(junc_tidy$seqnames == "M")) {
    print("Error! There's data for chr-MT!")
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
                       starts = as.integer(junc_tidy$donorSeqStart),
                       ends = as.integer(junc_tidy$donorSeqStop),
                       names = as.character(junc_tidy$junID),
                       scores = c(rep(".", nrow(junc_tidy))),
                       strands = junc_tidy$strand)
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
  
  #donor_sequences_tidy <- gsub(pattern = "\\::*\\)",
                              # x = as.character(donor_sequences_input$V1))
  #donor_sequences_tidy %>% head()
  
  #stopifnot(identical(gsub("\\(\\+\\)", "", gsub("\\(\\*\\)", "", gsub("\\(-\\)", "", as.character(donor_sequences_input$V1)))),
  #                    junc_tidy$junID %>% as.character()))
  junc_tidy <- cbind(junc_tidy, donor_sequence = as.character(donor_sequences_input$V2))
  
  junc_tidy %>% head()
  
  
  ## Get the acceptor genomic sequence
  
  to.BED <- data.frame(seqnames = junc_tidy$seqnames,
                       starts = as.integer(junc_tidy$AcceptorSeqStart),
                       ends = as.integer(junc_tidy$AcceptorSeqStop),
                       names = as.character(junc_tidy$junID),
                       scores = c(rep(".", nrow(junc_tidy))),
                       strands = junc_tidy$strand)
  
  
  tmp.file <- tempfile()
  
  write.table(to.BED, file = tmp.file, quote = F, sep = "\t", row.names = F, col.names = F)
  tmp.file_seq <- tempfile()
  system(paste0("bedtools getfasta -name -s -fi ", homo_sapiens_fasta_path, " -bed ",
                tmp.file, " -tab -fo ", tmp.file_seq))
  acceptor_sequences_input <- read.delim(tmp.file_seq, header = F)
  
  head(acceptor_sequences_input)
  head(donor_sequences_input)
  
  #stopifnot(identical(gsub("\\(\\+\\)", "", gsub("\\(\\*\\)", "", gsub("\\(-\\)", "", as.character(acceptor_sequences_input$V1)))),
  #                    junc_tidy$junID %>% as.character()))
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
  
  print("MaxEntScan score generated for the donor sequences!")
  
  
  ## get the maxentscan for the 3' splice site
  length(grep("N",as.character(junc_tidy$acceptor_sequence)))
  
  write.table(gsub("N","A",as.character(acceptor_sequences_input$V2)),file=tmp.file,row.names=F,col.names=F,quote=F)
  ss3score <- read.delim(pipe(paste0("perl ", max_ent_tool_path, "/score3.pl ", tmp.file)),header = F)
  identical(as.character(ss3score$V1),gsub("N","A",as.character(junc_tidy$acceptor_sequence)))
  junc_tidy <- cbind(junc_tidy, ss3score = ss3score$V2)
  
  print("MaxEntScan score generated for the acceptor sequences!")
  
  rm(ss5score, ss3score, tmp.file)
  
  junc_tidy[1,]
  
  
  
  return(junc_tidy)
}

add_cdts_cons_scores <- function(cluster = NULL,
                                 db_introns = NULL,
                                 folder_name = NULL) {
  
  
  if (is.null(db_introns) && !is.null(cluster)) {
    ## Load the IDB 
    db_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds")) %>%
      distinct(ref_junID, .keep_all = T)
  }
  
  db_introns <- db_introns %>%
    mutate(CDTS_5ss_mean = 0.0,
           CDTS_3ss_mean = 0.0,
           phastCons20way_5ss_mean = 0.0,
           phastCons20way_3ss_mean = 0.0) %>%
    GRanges()
  
  
  print(paste0(Sys.time(), " - getting 5' scores assigned to introns from IntroVerse ..."))
  
  ## https://www.nature.com/articles/nature09000
  ## Scores from the 5'ss ---------------------------------------------------------
  overlaps <- GenomicRanges::findOverlaps(query = CNC_CDTS_CONS_gr %>% diffloop::rmchr(),
                                          subject = GenomicRanges::GRanges(seqnames = db_introns %>% seqnames(),
                                                                           ranges = IRanges(start = db_introns %>% start() - 5,
                                                                                            end = db_introns %>% start() + 35),
                                                                           strand = db_introns %>% strand()),
                                          ignore.strand = FALSE,
                                          type = "any")
  
  overlaps_tidy <- overlaps %>%
    as.data.frame() %>%
    mutate(CDTS = CNC_CDTS_CONS_gr[queryHits(overlaps),]$CDTS,
           mean_phastCons20way = CNC_CDTS_CONS_gr[queryHits(overlaps),]$mean_phastCons20way) %>%
    group_by(subjectHits) %>%
    mutate(CDTS_mean = CDTS %>% mean(),
           mean_phastCons20way_mean = mean_phastCons20way %>% mean())
  
  db_introns[subjectHits(overlaps),]$CDTS_5ss_mean <- overlaps_tidy$CDTS_mean
  db_introns[subjectHits(overlaps),]$phastCons20way_5ss_mean <- overlaps_tidy$mean_phastCons20way_mean
  
  
  
  print(paste0(Sys.time(), " - getting 3' scores assigned to introns from IntroVerse ..."))
  
  ## Scores from the 3'ss ---------------------------------------------------------
  overlaps <- GenomicRanges::findOverlaps(query = CNC_CDTS_CONS_gr %>% diffloop::rmchr(),
                                          subject = GenomicRanges::GRanges(seqnames = db_introns %>% seqnames(),
                                                                           ranges = IRanges(start = db_introns %>% end() - 35,
                                                                                            end = db_introns %>% end() + 5),
                                                                           strand = db_introns %>% strand()),
                                          
                                          ignore.strand = FALSE,
                                          type = "any")
  
  overlaps_tidy <- overlaps %>% 
    as.data.frame() %>%
    mutate(CDTS = CNC_CDTS_CONS_gr[queryHits(overlaps),]$CDTS,
           mean_phastCons20way = CNC_CDTS_CONS_gr[queryHits(overlaps),]$mean_phastCons20way) %>%
    as_tibble() %>%
    group_by(subjectHits) %>%
    mutate(CDTS_mean = CDTS %>% mean(),
           mean_phastCons20way_mean = mean_phastCons20way %>% mean())
  
  db_introns[subjectHits(overlaps),]$CDTS_3ss_mean <- overlaps_tidy$CDTS_mean
  db_introns[subjectHits(overlaps),]$phastCons20way_3ss_mean <- overlaps_tidy$mean_phastCons20way_mean
  
  #####################
  ## SAVE RESULTS
  #####################
  
  if (!is.null(cluster)) {
    file_name <- paste0(folder_name, "/", cluster, "_db_introns.rds")
    saveRDS(object = db_introns %>% data.table::as.data.table(),
            file = file_name)
    
    print(paste0(Sys.time(), " - CDTS and Conservation scores added! IDB updated!"))
    
    
    rm(overlaps_tidy)
    rm(db_introns)
    rm(file_name)
    rm(overlaps)
    
  } else {
    return(db_introns)
  }
  
  
  
  
  #gc()
}

remove_MT_genes <- function(cluster,
                            folder_name) {
  
  
  
  ## Load mis-splicing ratios df
  df_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
  df_introns %>% head()
  df_introns %>% nrow()
  
  
  ## Load MT genes
  
  print(paste0(Sys.time(), " - checking the existance of MT genes..."))
  
  MT_genes <- readRDS(file = paste0(dependencies_folder,
                                    "/MT_genes/MT_genes.rds"))
  MT_geneID <- MT_genes %>% distinct(gene_id) %>% pull(gene_id)
  
  if (any(df_introns$gene_id %in% MT_geneID)) {
    
    ## Remove MT genes
    df_introns <- df_introns %>%
      dplyr::filter(!(gene_id %in% MT_geneID)) 
    
    print(paste0(Sys.time(), " - MT genes removed! ", df_introns %>% nrow(), " final number of junctions."))
    
    
    ## SAVE RESULTS
    saveRDS(object = df_introns,
            file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
    
    print(paste0(Sys.time(), " - results saved!"))
    
  } else {
    print(paste0(Sys.time(), " - the dataset doesn't contain any MT gene!"))
  }
  
  
  ## FREE SOME MEMORY
  rm(df_introns)
  rm(MT_geneID)
  rm(MT_genes)
  gc()
  
}


#################################################################
## SQL HELPER
## auxiliary functions to help in the generation of a database
#################################################################




