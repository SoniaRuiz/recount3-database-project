
#' Title
#' Creates the 'intron' and 'novel' master tables
#' @param database.path 
#' @param gtf.version 
#' @param database.folder 
#' @param results.folder 
#'
#' @return
#' @export
#'
#' @examples
SqlCreateMasterTables <- function(database.path,
                                  gtf.version,
                                  database.folder,
                                  results.folder,
                                  dependencies.folder,
                                  discard.minor.introns = F) {
  
  
  
  ####################################
  ## A) CREATE MASTER 'INTRON' TABLE
  ####################################
  
  SqlCreateMasterTableIntron(database.path,
                             gtf.version,
                             database.folder,
                             results.folder,
                             dependencies.folder,
                             discard.minor.introns = F)
  
  #############################################
  ## B) CREATE MASTER 'NOVEL' TABLE
  ## It contains novel 5' and 3' splicing events
  #############################################
  
  SqlCreateMasterTableNovel(database.path,
                            gtf.version,
                            database.folder,
                            results.folder,
                            dependencies.folder,
                            discard.minor.introns = F)
  
  
  #############################################
  ## C) CREATE MASTER 'COMBO' TABLE
  #############################################
 
  SqlCreateMasterTableCombo(database.path,
                            gtf.version,
                            database.folder,
                            results.folder,
                            dependencies.folder,
                            discard.minor.introns = F)
  
}

#' Title
#' Creates the master 'Intron' table, which stores information about all annotated introns found
#' across all samples studied
#' @param database.path 
#' @param gtf.version 
#' @param database.folder 
#' @param results.folder 
#' @param dependencies.folder 
#' @param discard.minor.introns 
#'
#' @return
#' @export
#'
#' @examples
SqlCreateMasterTableIntron <- function(database.path,
                                       gtf.version,
                                       database.folder,
                                       results.folder,
                                       dependencies.folder,
                                       discard.minor.introns = F) {
  
  logger::log_info("loading GRCh38 reference...")
  
  hg38 <- rtracklayer::import(con = paste0(dependencies.folder, "/Homo_sapiens.GRCh38.",gtf.version,".chr.gtf"))
  
  ##########################################
  ## LOAD AND TIDY THE PAIR-WISE DISTANCES
  ##########################################
  
  con <- dbConnect(RSQLite::SQLite(), database.path)
  
  if ( file.exists(paste0(database.folder, "/all_jxn_correct_pairings.rds")) ) {
    
    logger::log_info("loading the pre-generated data...")
    
    df_all_distances_pairings <- readRDS(file = paste0(database.folder, "/all_jxn_correct_pairings.rds"))
    df_ambiguous_novel <- readRDS(file = paste0(database.folder, "/all_jxn_ambiguous_pairings.rds"))
    df_introns_never <- readRDS(file = paste0(database.folder, "/all_jxn_never_misspliced.rds"))
    df_introns_not_paired <- readRDS(file = paste0(database.folder, "/all_jxn_not_paired.rds"))
    
  } else {
    
    stop("ERROR loading file dependencies!")
  }
  
  ####################################
  ## A) GET ALL INTRONS
  ####################################
  
  ## Get all annotated introns
  logger::log_info("getting the mis-spliced introns...")
  
  df_all_introns <- df_all_distances_pairings %>%
    distinct(ref_junID, .keep_all = T) %>%
    as_tibble() %>%
    mutate(misspliced = "Yes") %>%
    dplyr::select(ref_junID,
                  seqnames = ref_seq,
                  start = ref_start,
                  strand = ref_strand,
                  end = ref_end,
                  gene_id, 
                  tx_id_junction,
                  misspliced) %>%
    GRanges() %>% ## This is to get the width() as calculated by GRanges
    as_tibble()
  
  
  logger::log_info("getting the never mis-spliced introns...")
  
  ## Remove potential * in the junID of the reference introns
  if (any(str_detect(string = df_introns_never$ref_junID, pattern = "\\*"))) {
    logger::log_info("ERROR! some never mis-spliced junctions still have an *!")
    break;
  }
  
  df_introns_never_tidy <- df_introns_never %>% 
    as_tibble() %>%
    mutate(misspliced = "No") #%>%
  #dplyr::filter(!(ref_junID %in% df_all_introns$ref_junID)) #%>%
  #dplyr::filter(!(ref_junID %in% df_ambiguous_novel$ref_junID)) 
  
  
  
  logger::log_info("getting the not paired mis-spliced introns...")
  
  ## Remove potential * in the junID of the reference introns
  if (any(str_detect(string = df_introns_never$ref_junID, pattern = "\\*"))) {
    logger::log_info("ERROR! some never mis-spliced junctions still have an *!")
    break;
  }
  
  
  
  df_introns_not_paired_tidy <- df_introns_not_paired %>% 
    as_tibble() %>%
    mutate(misspliced = "Potentially") #%>%
  
  
  
  if ( intersect(df_introns_not_paired_tidy$ref_junID, df_all_introns$ref_junID) %>% length() > 0 ) {
    logger::log_info("ERROR! Non-paired introns classified as paired!")
    break;
  }
  
  
  
  df_introns_parenting_ambiguous <- df_ambiguous_novel %>% 
    as_tibble() %>%
    mutate(misspliced = "Potentially") %>%
    dplyr::filter(!(ref_junID %in% df_all_introns$ref_junID)) %>%
    distinct(ref_junID, .keep_all=T)
  
  
  
  if ( intersect(df_introns_parenting_ambiguous$ref_junID, df_introns_not_paired_tidy$ref_junID) %>% length() > 0 ) {
    logger::log_info("ERROR! Introns parenting ambiguous junctions are classified as non-paired introns!")
    break;
  }
  
  if ( intersect(df_introns_parenting_ambiguous$ref_junID, df_all_introns$ref_junID) %>% length() > 0 ) {
    logger::log_info("ERROR! Introns parenting ambiguous junctions are classified as paired introns!")
    break;
  }
  
  
  ## QC
  
  if (any(df_all_introns$width %>% abs() < 25)) {
    logger::log_info("ERROR! some mis-spliced introns are shorter than 25bp!")
    break;
  }
  if (any(df_introns_never_tidy$width %>% abs() < 25)) {
    logger::log_info("ERROR! some never mis-spliced introns are shorter than 25bp!")
    break;
  }
  if (any(df_introns_not_paired_tidy$width %>% abs() < 25)) {
    logger::log_info("ERROR! some POTENTIALLY mis-spliced introns are shorter than 25bp!")
    break;
  }
  if (any(df_introns_parenting_ambiguous$width %>% abs() < 25)) {
    logger::log_info("ERROR! some POTENTIALLY mis-spliced introns are shorter than 25bp!")
    break;
  }
  if ( intersect(df_introns_never_tidy$ref_junID, df_all_introns$ref_junID) %>% length() > 0) {
    logger::log_info("ERROR! some never mis-spliced introns are mis-spliced!")
    break;
  }
  if ( intersect(df_introns_not_paired_tidy$ref_junID, df_all_introns$ref_junID) %>% length() > 0) {
    logger::log_info("ERROR! some POTENTIALLY mis-spliced introns are CLASSIFIED AS mis-spliced!")
    break;
  }
  
  
  ## Merge mis-spliced introns and never mis-spliced introns
  df_introns_introverse <- rbind(df_introns_never_tidy, 
                                 df_introns_not_paired_tidy, 
                                 df_introns_parenting_ambiguous,
                                 df_all_introns)
  
  if (any(str_detect(string = df_introns_introverse$ref_junID, pattern = "\\*" ))) {
    logger::log_info("ERROR! Some introns still have an ambiguous '*' strand!")
    break;
  }
  if ( any(df_introns_introverse$ref_junID %>% duplicated()) ) {
    logger::log_info("ERROR! Duplicated IDs")
    break;
  }
  
  #saveRDS(object = df_introns_introverse %>%
  #          distinct(ref_junID, .keep_all = T),
  #        file = paste0(folder_database, "/df_all_introns_database.rds"))
  
  
  
  ######################################
  ## INTRONS - REMOVE AMBIGUOUS INTRONS
  ######################################
  
  logger::log_info("getting the ambiguous introns...")
  
  df_introns_introverse_tidy <- df_introns_introverse %>%
    distinct(ref_junID, .keep_all = T) %>% 
    rowwise() %>%
    mutate(ambiguous = ifelse(gene_id %>% unlist() %>% length() > 1, T, F))
  
  ## There should not be any ambiguous intron at this point
  if ( any(df_introns_introverse_tidy %>% dplyr::filter(ambiguous == T) %>% nrow() > 0) ) {
    logger::log_info("ERROR! Still there are some ambiguous introns")
    break;
  } else {
    df_introns_introverse_tidy <- df_introns_introverse_tidy %>%
      dplyr::select(-ambiguous)
  }
  
  logger::log_info(df_introns_introverse_tidy %>% distinct(ref_junID) %>% nrow(), " annotated introns TO BE STORED IN THE DATABASE...")
  
  
  # df_introns_introverse_tidy %>% dplyr::count(misspliced)
  
  ######################################
  ## GENES - CREATE GENE TABLE
  ######################################
  
  SqlCreateMasterTableGene(database.path = database.path,
                           hg38 = hg38,
                           gene.ids = df_introns_introverse_tidy %>% unnest(gene_id) %>% distinct(gene_id) )
  
  
  ######################################
  ## TX_JUNCTION - CREATE TX TABLE
  ######################################
  
  SqlCreateMasterTableTranscript(database.path = database.path,
                                 gene.ids = df_introns_introverse_tidy %>% unnest(gene_id) %>% distinct(gene_id),
                                 dependencies.folder = dependencies.folder,
                                 hg38 = hg38,
                                 tx.ids = df_introns_introverse_tidy %>% unnest(tx_id_junction) %>% distinct(tx_id_junction))
  
  
  ############################################
  ## INTRONS - ADD THE TRANSCRIPT FOREING KEY
  ############################################
  
  logger::log_info(" --> adding TRANSCRIPT foreing key to the introns...")
  
  query <- paste0("SELECT id, transcript_id FROM 'transcript'")
  df_transcripts <- dbGetQuery(con, query) %>% as_tibble()
  
  ## Add the GENE ID for the foreign key
  df_introns_introverse_tidy <- df_introns_introverse_tidy %>%
    unnest(tx_id_junction) %>% 
    left_join(df_transcripts, by = c("tx_id_junction" = "transcript_id")) %>%
    dplyr::select(-gene_id, -tx_id_junction) %>%
    dplyr::rename(transcript_id = id)  %>%
    group_by(ref_junID) %>%
    mutate(transcript_id_list = paste(transcript_id, collapse = ",")) %>%
    ungroup() %>%
    distinct(ref_junID, .keep_all = T)
  
  
  logger::log_info(df_introns_introverse_tidy %>% distinct(ref_junID) %>% nrow(), " introns to store...")
  
  ######################################
  ## INTRONS - ADD MAXENTSCAN INFO 
  ######################################
  
  logger::log_info(" --> adding the MaxEntScan info ...")
  
  wd <- getwd()
  if ( !file.exists(file.path(dependencies.folder, "Homo_sapiens.GRCh38.dna.primary_assembly.fa")) ) {
    logger::log_info(paste0("ERROR! File dependency 'Homo_sapiens.GRCh38.dna.primary_assembly.fa' does not exist within the specified dependencies folder."))
    break;
  }
  
  ## Add MaxEntScan score to the split reads
  all_split_reads_tidy <- GenerateMaxEntScore(junc_tidy = df_introns_introverse_tidy %>% dplyr::rename(junID = ref_junID) %>% distinct(junID, .keep_all = T),
                                              max_ent_tool_path = paste0(dependencies.folder, "/fordownload/"),
                                              homo_sapiens_fasta_path = paste0(dependencies.folder, "/Homo_sapiens.GRCh38.dna.primary_assembly.fa") )
  rm(df_introns_introverse_tidy)
  gc()
  
  all_split_reads_tidy <- all_split_reads_tidy %>% as_tibble()
  
  all_split_reads_tidy <- all_split_reads_tidy %>% 
    dplyr::select(-c(donorSeqStart, donorSeqStop, AcceptorSeqStart, AcceptorSeqStop)) %>%
    dplyr::rename(ref_donor_sequence = donor_sequence,
                  ref_acceptor_sequence = acceptor_sequence)
  
  setwd(wd)
  
  df_all_introns <- all_split_reads_tidy %>%
    mutate(ref_junID = paste0("chr", seqnames, ":", start, "-", end, ":", strand)) 
  
  if ((setdiff(df_all_introns$junID, df_all_introns$ref_junID) %>% length()) > 0) {
    logger::log_info("ERROR! ")
    break;
  }
  
  df_all_introns %>% as_tibble()
  
  df_all_introns <- df_all_introns %>%
    dplyr::select(-one_of("junID", "ref_ss5score", "ref_ss3score")) %>% 
    dplyr::rename(ref_mes5ss = ss5score, ref_mes3ss = ss3score) %>%
    dplyr::relocate(ref_junID) %>%
    dplyr::relocate(c(ref_mes5ss, ref_mes3ss), .before = transcript_id ) 
  
  df_all_introns %>% as_tibble()
  
  ################################################
  ## INTRONS - ADD THE CONSERVATION AND CDTS INFO
  ################################################
  
  df_all_introns %>% as_tibble()
  logger::log_info("adding CDTS and Conservation scores...")
  
  df_all_introns <- GenerateCdtsPhastconsScores(dependencies.folder = dependencies.folder,
                                                db.introns = df_all_introns %>% distinct(ref_junID, .keep_all = T) %>%
                                                  dplyr::rename("junID" = "ref_junID") %>% as_tibble(),
                                                intron.size = 100, phastcons.type = 17) %>% as_tibble()
  
  df_all_introns <- df_all_introns %>% dplyr::rename("ref_junID" = "junID") %>% mutate_if(is.numeric, ~replace_na(., 0))
  df_all_introns %>% as_tibble()
  
  ######################################
  ## INTRONS - ADD CLINVAR DATA
  ######################################
  
  logger::log_info("adding the ClinVar data...")
  
  df_all_introns_gr <- AddClinvarData(df.all.introns = df_all_introns, dependencies.folder)
  
  df_all_introns_gr %>% head()
  df_all_introns_gr %>% as_tibble() %>% filter(clinvar == T) %>% dplyr::select(ref_junID, clinvar)
  
  ######################################
  ## INTRONS - ADD THE INTRON TYPE
  ## This intron type corresponds to whether the intron is 
  ## spliced out by the minor or the major spliceosome
  ######################################
  
  logger::log_info("adding the IAOD intron data...")
  
  ## Load intron type files
  u12_introns_gr <- readRDS(file = file.path(dependencies.folder, "minor_introns_tidy.rds")) %>% GRanges() %>% diffloop::addchr()
  
  ## Add a new column to incorporate info about intron type
  elementMetadata(df_all_introns_gr)[, "u2_intron"] = T
  
  ## MINOR INTRON
  logger::log_info("Getting junctions spliced out by the minor spliceosome.")
  overlaps <- GenomicRanges::findOverlaps(query = u12_introns_gr, subject = df_all_introns_gr,
                                          ignore.strand = FALSE, type = "equal")
  
  logger::log_info(queryHits(overlaps) %>% length(), " introns spliced by the minor spliceosome!")
  df_all_introns_gr[subjectHits(overlaps),]$u2_intron <- F
  
  if ( which(df_all_introns_gr$u2_intron == F) %>% length() != 
       df_all_introns_gr[subjectHits(overlaps),] %>% length() ) {
    logger::log_info("ERROR! Minor spliceosomal introns have not been added adequately!")
    break;
  }
  
  ## Discard introns spliced out by the minor spliceosome. Comment out to avoid applying this filter.
  df_all_introns_tidy <- df_all_introns_gr %>%
    as_tibble()
  
  if (discard.minor.introns) {
    df_all_introns_tidy <- df_all_introns_tidy %>%
      filter(u2_intron == T) %>%
      dplyr::select(-u2_intron)
  }
  
  logger::log_info(df_all_introns_tidy$ref_junID %>% unique %>% length(), " introns to be stored!")
  
  #########################################
  ## INTRONS - ADD THE TRANSCRIPT BIOTYPE
  #########################################
  
  df_biotype_junID <- readRDS(file = file.path(database.folder, "all_split_reads_qc_level1_PC_biotype.rds")) %>% as_tibble()
  
  if (any(str_detect(df_biotype_junID$junID, pattern = "\\*"))){stop("Still junctions with * as strand!")}
  
  df_all_introns_tidy <- df_all_introns_tidy %>%
    inner_join(y = df_biotype_junID %>% dplyr::select(junID, protein_coding),
               by = c("ref_junID" = "junID"))
  
  logger::log_info(df_all_introns_tidy$ref_junID %>% unique %>% length(), " introns to be stored!")
  
  ######################################
  ## INTRONS - POPULATE THE TABLE
  ######################################
  
  # mean_phastCons7way5ss_35 DOUBLE NOT NULL, 
  # mean_phastCons7way3ss_35 DOUBLE NOT NULL, 
  # mean_phastCons20way5ss_35 DOUBLE NOT NULL, 
  # mean_phastCons20way3ss_35 DOUBLE NOT NULL, 
  # mean_phastCons4way5ss_70 DOUBLE NOT NULL, 
  # mean_phastCons4way3ss_70 DOUBLE NOT NULL, 
  # mean_phastCons7way5ss_70 DOUBLE NOT NULL, 
  # mean_phastCons7way3ss_70 DOUBLE NOT NULL, 
  # mean_phastCons20way5ss_70 DOUBLE NOT NULL, 
  # mean_phastCons20way3ss_70 DOUBLE NOT NULL, 
  # mean_CDTS5ss_70 DOUBLE NOT NULL, 
  # mean_CDTS3ss_70 DOUBLE NOT NULL, 
  # mean_phastCons4way5ss_100 DOUBLE NOT NULL, 
  # mean_phastCons4way3ss_100 DOUBLE NOT NULL, 
  # mean_phastCons7way5ss_100 DOUBLE NOT NULL, 
  # mean_phastCons7way3ss_100 DOUBLE NOT NULL, 
  # mean_phastCons20way5ss_100 DOUBLE NOT NULL, 
  # mean_phastCons20way3ss_100 DOUBLE NOT NULL, 
  # mean_CDTS5ss_100 DOUBLE NOT NULL, 
  # mean_CDTS3ss_100 DOUBLE NOT NULL, 
  
  df_all_introns_tidy %>% names()
  
  query <- paste0("CREATE TABLE IF NOT EXISTS 'intron'",
                  "(ref_junID NUMERIC PRIMARY KEY NOT NULL,
                  ref_coordinates TEXT NOT NULL, 
                  
                  seqnames TEXT NOT NULL,
                  start NUMERIC NOT NULL,
                  end NUMERIC NOT NULL,
                  strand TEXT NOT NULL, 
                  
                  ref_length INTEGER NOT NULL, 
                  ref_mes5ss DOUBLE NOT NULL, 
                  ref_mes3ss DOUBLE NOT NULL, 
                  
                  mean_phastCons17way5ss_100 DOUBLE NOT NULL, 
                  mean_phastCons17way3ss_100 DOUBLE NOT NULL, 
                  
                  mean_CDTS5ss_100 DOUBLE NOT NULL, 
                  mean_CDTS3ss_100 DOUBLE NOT NULL, 
                  
                  ref_donor_sequence TEXT NOT NULL,
                  ref_acceptor_sequence TEXT NOT NULL,
                  
                  u2_intron BOOL,
                  
                  clinvar BOOL NOT NULL, 
                  
                  protein_coding DOUBLE NOT NULL, 
                  
                  misspliced TEXT NOT NULL,
                  transcript_id INTEGER NOT NULL,
                  transcript_id_list TEXT NOT NULL,
                  FOREIGN KEY (transcript_id) REFERENCES 'transcript'(id))")
  
  
  con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = database.path)
  # dbListTables(con)
  # DBI::dbRemoveTable(conn = con, 'intron')
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  logger::log_info("'Intron' table created!")
  
  ## POPULATE INTRON TABLE ----------------------------------------------------
  df_all_introns_tidy_final <- df_all_introns_tidy %>% 
    dplyr::rename(ref_length = width, ref_coordinates = ref_junID) %>%
    distinct(ref_coordinates, .keep_all = T) %>%
    tibble::rowid_to_column("ref_junID")
  
  if (all(df_all_introns_tidy_final$misspliced == T)) {
    logger::log_info("ERROR! all introns classified as misspliced!")
    break;
  }
  if (any(df_all_introns_tidy_final$transcript_id %>% is.na())) {
    logger::log_info("ERROR! some introns do not have a gene assigned")
    break;
  }
  
  DBI::dbAppendTable(conn = con, name = "intron", value = df_all_introns_tidy_final)
  
  logger::log_info("'Intron' master table populated! ", df_all_introns_tidy_final %>% distinct(ref_junID) %>% nrow(), " annotated introns stored!" )
  
  ## CREATE INDEXES TO SPEED UP QUERIES ------------------------------------------
  query <- paste0("CREATE UNIQUE INDEX 'index_intron' ON 'intron'(ref_junID)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  query <- paste0("CREATE UNIQUE INDEX 'index_intron_coord' ON 'intron'(ref_junID,ref_coordinates)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  query <- paste0("CREATE UNIQUE INDEX 'index_intron_position' ON 'intron'(seqnames,start,end,strand)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
}



#' Title
#' Creates the master 'Novel' table, which stores information about all novel 5' and 3' splicing events found
#' across all samples studied
#' @param database.path 
#' @param gtf.version 
#' @param database.folder 
#' @param results.folder 
#' @param dependencies.folder 
#' @param discard.minor.introns 
#'
#' @return
#' @export
#'
#' @examples
SqlCreateMasterTableNovel <- function(database.path,
                                      gtf.version,
                                      database.folder,
                                      results.folder,
                                      dependencies.folder,
                                      discard.minor.introns = F) {
  
  logger::log_info("loading GRCh38 reference...")
  
  hg38 <- rtracklayer::import(con = paste0(dependencies.folder, "/Homo_sapiens.GRCh38.",gtf.version,".chr.gtf"))
  
  ##########################################
  ## LOAD AND TIDY THE PAIR-WISE DISTANCES
  ##########################################
  
  con <- dbConnect(RSQLite::SQLite(), database.path)
  
  if ( file.exists(paste0(database.folder, "/all_jxn_correct_pairings.rds")) ) {
    
    logger::log_info("loading the pre-generated data...")
    
    df_all_distances_pairings <- readRDS(file = paste0(database.folder, "/all_jxn_correct_pairings.rds"))
    df_ambiguous_novel <- readRDS(file = paste0(database.folder, "/all_jxn_ambiguous_pairings.rds"))
    df_introns_never <- readRDS(file = paste0(database.folder, "/all_jxn_never_misspliced.rds"))
    df_introns_not_paired <- readRDS(file = paste0(database.folder, "/all_jxn_not_paired.rds"))
    
  } else {
    stop("ERROR loading file dependencies!")
  }
  
  logger::log_info("adding MaxEntScan scores to the NOVEL JUNCTONS...")
  
  df_all_novel_raw_tidy <- df_all_distances_pairings %>%
    mutate(start = novel_start %>% as.integer(),
           end = novel_end %>% as.integer()) %>%
    dplyr::select(seqnames = novel_seq,
                  start, end, 
                  strand = novel_strand,
                  novel_junID, ref_junID, 
                  novel_type = type, distance) %>%
    distinct(novel_junID, .keep_all = T) %>%
    filter(ref_junID %in% df_all_introns_tidy_final$ref_coordinates)
  
  wd <- getwd()
  
  ## Add MaxEntScan score to the split reads
  all_split_reads_tidy <- GenerateMaxEntScore(junc_tidy = df_all_novel_raw_tidy %>% dplyr::rename(junID = novel_junID) %>% distinct(junID, .keep_all = T),
                                              max_ent_tool_path = paste0(dependencies.folder,"/fordownload/"),
                                              homo_sapiens_fasta_path = paste0(dependencies.folder, "/Homo_sapiens.GRCh38.dna.primary_assembly.fa"))
  
  all_split_reads_tidy %>% as_tibble()
  
  all_split_reads_tidy <- all_split_reads_tidy %>% 
    dplyr::select(-c(donorSeqStart, donorSeqStop, AcceptorSeqStart, AcceptorSeqStop)) %>%
    dplyr::rename(novel_donor_sequence = donor_sequence, novel_acceptor_sequence = acceptor_sequence)
  
  setwd(wd)
  
  
  df_all_novels_tidy <- all_split_reads_tidy %>%
    mutate(novel_junID = paste0("chr", seqnames, ":", start, "-", end, ":", strand)) 
  
  if ((setdiff(df_all_novels_tidy$junID, df_all_novels_tidy$novel_junID) %>% length()) > 0) {
    stop("ERROR! Novel junctions have been analysed under different sorting.")
  }
  
  rm(all_split_reads_tidy)
  gc()
  
  df_all_novels_tidy <- df_all_novels_tidy %>%
    dplyr::select(-c(any_of(c("junID","novel_ss5score","novel_ss3score")))) %>% 
    dplyr::rename(novel_mes5ss = ss5score, novel_mes3ss = ss3score) %>%
    dplyr::relocate(ref_junID, novel_junID) %>%
    dplyr::relocate(c(novel_mes5ss, novel_mes3ss), .before = novel_type ) %>% 
    as_tibble()
  
  ################################################
  ## NOVEL - ADD THE CONSERVATION AND CDTS INFO
  ################################################
  
  df_all_novels_tidy %>% as_tibble()
  logger::log_info("\t Adding CDTS and Conservation scores to the novel junctions ...")
  
  df_all_novels_tidy <- GenerateCdtsPhastconsScores(dependencies.folder = dependencies.folder,
                                                    db.introns = df_all_novels_tidy %>%
                                                      distinct(novel_junID, .keep_all = T) %>% dplyr::rename(junID = "novel_junID") %>% as_tibble(),
                                                    intron.size = 100,
                                                    phastcons.type = 17) %>% 
    as_tibble() %>% 
    dplyr::rename("novel_junID" = "junID") %>% 
    mutate_if(is.numeric, ~replace_na(., 0))
  
  df_all_novels_tidy %>% as_tibble()
  
  
  ######################################
  ## NOVEL - ADD CLINVAR DATA
  ######################################
  
  logger::log_info("adding the ClinVar data...")
  
  df_all_novels_tidy_gr <- AddClinvarData(df.all.introns = df_all_novels_tidy, dependencies.folder)
  
  df_all_novels_tidy_gr %>% head()
  df_all_novels_tidy_gr %>% as_tibble() %>% filter(clinvar == T) %>% dplyr::select(novel_junID, clinvar)
  
  ##############################################
  ## NOVEL - ADD INTRON FOREIGN KEY REFERENCE
  ##############################################
  
  logger::log_info("adding the INTRON foreign key to the NOVEL JUNCTONS...")
  
  df_all_novels_tidy_final <- df_all_novels_tidy_gr %>% 
    as_tibble() %>%
    inner_join(y = df_all_introns_tidy_final %>% dplyr::select(ref_junID, ref_coordinates),
               by = c("ref_junID" = "ref_coordinates" )) %>%
    dplyr::select(-ref_junID) %>%
    dplyr::rename(ref_junID = ref_junID.y) %>% 
    dplyr::relocate(ref_junID)
  
  rm(df_all_novels_tidy)
  rm(df_all_novels_tidy_gr)
  
  df_all_novels_tidy_final %>% filter(clinvar == T)
  
  ####################################
  ## CREATE NOVEL JUNCTION TABLE
  ####################################
  
  # dbRemoveTable(conn = con, "novel")
  query <- paste0("CREATE TABLE IF NOT EXISTS 'novel'",
                  "(novel_junID NUMERIC NOT NULL,
                  ref_junID NUMERIC NOT NULL,

                  seqnames TEXT NOT NULL,
                  start NUMERIC NOT NULL,
                  end NUMERIC NOT NULL,
                  strand TEXT NOT NULL, 

                  novel_coordinates TEXT NOT NULL, 
                  novel_mes5ss DOUBLE NOT NULL, 
                  novel_mes3ss DOUBLE NOT NULL,
                  
                  mean_phastCons17way5ss_100 DOUBLE NOT NULL, 
                  mean_phastCons17way3ss_100 DOUBLE NOT NULL, 
                  
                  mean_CDTS5ss_100 DOUBLE NOT NULL, 
                  mean_CDTS3ss_100 DOUBLE NOT NULL,

                  clinvar BOOL NOT NULL, 

                  novel_donor_sequence TEXT NOT NULL,
                  novel_acceptor_sequence TEXT NOT NULL,
              
                  novel_length INTEGER NOT NULL, 
                  novel_type TEXT NOT NULL, 
                  distance INTEGER NOT NULL,

                  PRIMARY KEY (ref_junID, novel_junID),
                  FOREIGN KEY (ref_junID) REFERENCES 'intron'(ref_junID))")
  
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  logger::log_info("'Novel' master table created!")
  
  df_all_novels_tidy_final <- df_all_novels_tidy_final %>% 
    dplyr::rename(novel_coordinates = novel_junID ) %>%
    distinct(novel_coordinates, .keep_all = T) %>%
    tibble::rowid_to_column("novel_junID") %>%
    dplyr::rename(novel_length = width)
  
  
  if (any(duplicated(df_all_novels_tidy_final$novel_coordinates))) {
    logger::log_info("ERROR! some novel junctions are duplicated")
    break;
  }
  
  DBI::dbAppendTable(conn = con, name = "novel", value = df_all_novels_tidy_final)
  
  logger::log_info("'Novel' master table populated! ", 
                   df_all_novels_tidy_final %>% distinct(novel_junID) %>% nrow(), " novel junctions stored!" )
  
  
  
  ## CREATE INDEXES TO SPEED UP QUERIES ------------------------------------------
  query <- paste0("CREATE UNIQUE INDEX 'index_novel' ON 'novel'(ref_junID,novel_junID)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  query <- paste0("CREATE UNIQUE INDEX 'index_novel_coord' ON 'novel'(novel_coordinates)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  query <- paste0("CREATE UNIQUE INDEX 'index_novel_position' ON 'novel'(seqnames,start,end,strand)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
}

#' Title
#' Creates the Master Novel Combo table, which stores information about all novel combination splicing events found
#' across the samples studied
#' @param database.path Local path to the .sqlite database
#' @param recount3.project.IDs Vector with the ID of the recount3 projects to work with
#' @param database.folder Local path to the folder that contains the database
#' @param results.folder Local path to the folder that contains the result files
#'
#' @return
#' @export
#'
#' @examples
SqlCreateMasterTableCombo <- function(database.path,
                                      database.folder,
                                      results.folder,
                                      dependencies.folder,
                                      con = NULL,
                                      recount3.project.IDs = NULL) {
  
  if (is.null(con)) {
    con <- dbConnect(RSQLite::SQLite(), database.path)
  }
  
  tables <- DBI::dbListTables(conn = con)
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  
  logger::log_info( " --> SQL connection stablished!") 
  logger::log_info("Querying master tables ...")
  
  ## GET FROM MASTER TABLE
  df_metadata <- dbGetQuery(con, paste0("SELECT * FROM 'metadata'")) 
  
  ## GET FROM GENE TABLE
  master_gene <- dbGetQuery(con, paste0("SELECT * FROM 'gene'")) %>% as_tibble()
  master_gene %>% nrow()
  
  ## GET FROM TRANSCRIPT TABLE
  master_transcript <- dbGetQuery(con, paste0("SELECT * FROM 'transcript'")) %>% as_tibble()
  master_transcript %>% nrow()
  
  DBI::dbDisconnect(conn = con) 
  
  
  if ( is.null(recount3.project.IDs) ){
    recount3.project.IDs <- (df_metadata$SRA_project %>% unique())
  }
  
  ## Get all novel combos across all tables
  all_split_reads_combos <- map_df(recount3.project.IDs, function(project_id) { 
    
    # project_id <- recount3.project.IDs[1]
    
    logger::log_info(" --> Working with '", project_id, "' ...")
    results_folder_local <- paste0(results.folder, "/", project_id, "/")
    
    clusters <- df_metadata %>% dplyr::filter(SRA_project == project_id) %>% distinct(cluster) %>% pull()
    
    map_df(clusters, function(cluster_id) { 
      
      # cluster_id <- clusters[1]
      logger::log_info(project_id, " --> ", cluster_id)
      
      if (file.exists(paste0(results_folder_local, "/base_data/", project_id, "_", cluster_id, "_all_split_reads_combos.rds"))) {
        
        ## Load all split reads
        readRDS(file = paste0(results_folder_local, "/base_data/", project_id, "_", cluster_id, "_all_split_reads_combos.rds")) %>%
          dplyr::select(-any_of(c("in_ref","n_projects","annotated"))) %>%
          unnest(gene_id)
      }
    })
  })
  
  all_split_read_combos_RBPs <- all_split_reads_combos %>% distinct(junID, .keep_all = T)
  
        
  ############################################
  ## ADD THE TRANSCRIPT FOREING KEY
  ############################################
  
  logger::log_info("\t Adding TRANSCRIPT foreing key to the novel combos...")

  
  all_split_read_combos_RBPs_tidy <- all_split_read_combos_RBPs %>% 
    dplyr::select(-gene_id) %>%
    unnest(tx_id_junction) %>% 
    left_join(master_transcript, by = c("tx_id_junction" = "transcript_id")) %>%
    dplyr::select(-gene_id, -tx_id_junction, -TSL, -MANE, -reads) %>%
    dplyr::rename(transcript_id = id)  %>%
    drop_na(transcript_id) %>%
    group_by(junID) %>%
    mutate(transcript_id_list = paste(transcript_id, collapse = ",")) %>%
    ungroup() %>%
    distinct(junID, .keep_all = T)
  
  
  logger::log_info(all_split_read_combos_RBPs_tidy %>% distinct(junID) %>% nrow(), " novel combos to store....")
  summary(all_split_read_combos_RBPs_tidy)
  
        
  ######################################
  ## ADD MAXENTSCAN INFO 
  ######################################
  
  logger::log_info(" --> adding the MaxEntScan info ...")
  
  wd <- getwd()
  if ( !file.exists(file.path(dependencies.folder, "/Homo_sapiens.GRCh38.dna.primary_assembly.fa")) ) {
    stop("ERROR! File dependency 'Homo_sapiens.GRCh38.dna.primary_assembly.fa' does not exist within the specified dependencies folder.")
  }
  ## Add MaxEntScan score to the split reads
  all_split_read_combos_RBPs_w_MES <- GenerateMaxEntScore(junc_tidy = all_split_read_combos_RBPs_tidy,
                                                          max_ent_tool_path = file.path(dependencies.folder, "fordownload/"),
                                                          homo_sapiens_fasta_path = file.path(dependencies.folder, "Homo_sapiens.GRCh38.dna.primary_assembly.fa") ) %>% as_tibble()
  
  all_split_read_combos_RBPs_w_MES <- all_split_read_combos_RBPs_w_MES %>% 
    dplyr::select(-c(donorSeqStart, donorSeqStop, AcceptorSeqStart, AcceptorSeqStop)) %>%
    dplyr::rename(ref_donor_sequence = donor_sequence, ref_acceptor_sequence = acceptor_sequence) %>%
    dplyr::rename(ref_mes5ss = ss5score, ref_mes3ss = ss3score) %>%
    dplyr::relocate(c(ref_mes5ss, ref_mes3ss), .before = transcript_id ) %>% 
    as_tibble()
  
  setwd(wd)
        
        
  ################################################
  ## ADD THE CONSERVATION AND CDTS INFO
  ################################################
  
  all_split_read_combos_RBPs_w_MES %>% as_tibble()
  logger::log_info("adding CDTS and Conservation scores...")
  
  
  all_split_read_combos_RBPs_w_scores <- GenerateCdtsPhastconsScores(db.introns = all_split_read_combos_RBPs_w_MES,
                                                                     intron.size = 100,
                                                                     phastcons.type = 17,
                                                                     dependencies.folder = dependencies.folder,
                                                                     folder.name = ) %>% as_tibble() %>% mutate_if(is.numeric, ~replace_na(., 0))
  all_split_read_combos_RBPs_w_scores %>% as_tibble()
  
  #########################################
  ## ADD THE TRANSCRIPT BIOTYPE
  #########################################
  
  df_biotype_junID <- readRDS(file = file.path(database.folder, "all_split_reads_qc_level1_PC_biotype.rds")) %>% as_tibble()
  
  if (any(str_detect(df_biotype_junID$junID, pattern = "\\*"))) {
    stop("There are junctions in the transcript biotype file with * as strand!")
  }
  
  all_split_read_combos_RBPs_w_scores <- all_split_read_combos_RBPs_w_scores %>% 
    inner_join(y = df_biotype_junID %>% dplyr::select(junID, protein_coding), by = c("junID"))
  
  logger::log_info(all_split_read_combos_RBPs_w_scores$junID %>% unique %>% length(), " novel combos to be stored!")
  
  
  #########################################################
  ## CREATE AND POPULATE MASTER 'NOVEL COMBO' TABLE
  #########################################################
  
  logger::log_info( " --> creating 'novel combo' table ... ")
  
  ## 1. CREATE THE TABLE
  
  # dbRemoveTable(conn = con, paste0(cluster_id, "_", project_id))
  query <- paste0("CREATE TABLE IF NOT EXISTS 'combo'", 
                  "(ref_junID NUMERIC PRIMARY KEY NOT NULL,
                  ref_coordinates TEXT NOT NULL, 
                  
                  seqnames TEXT NOT NULL,
                  start NUMERIC NOT NULL,
                  end NUMERIC NOT NULL,
                  strand TEXT NOT NULL, 
                  ref_length INTEGER NOT NULL,
                  
                  type TEXT NOT NULL, 
                  
                  ref_mes5ss DOUBLE NOT NULL, 
                  ref_mes3ss DOUBLE NOT NULL, 
                  
                  mean_phastCons17way5ss_100 DOUBLE NOT NULL, 
                  mean_phastCons17way3ss_100 DOUBLE NOT NULL, 
                  
                  mean_CDTS5ss_100 DOUBLE NOT NULL, 
                  mean_CDTS3ss_100 DOUBLE NOT NULL, 
                  
                  ref_donor_sequence TEXT NOT NULL,
                  ref_acceptor_sequence TEXT NOT NULL,
                  
                  protein_coding DOUBLE NOT NULL, 
                  
                  transcript_id INTEGER NOT NULL,
                  transcript_biotype TEXT NOT NULL,
                  transcript_id_list TEXT NOT NULL,
                  
                  FOREIGN KEY (transcript_id) REFERENCES 'transcript'(id))")
        
        ## Connect the database
        con <- dbConnect(RSQLite::SQLite(), database.path)
        DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
        
        ## Create the NOVEL COMBO table
        res <- DBI::dbSendQuery(conn = con, statement = query)
        DBI::dbClearResult(res)
        

  ## 2. POPULATE THE TABLE
        
  all_split_read_combos_RBPs_w_scores_final <- all_split_read_combos_RBPs_w_scores %>%
    dplyr::rename(ref_coordinates = junID, ref_length = width) %>%
    tibble::rowid_to_column("ref_junID") %>% 
    mutate(transcript_id = transcript_id %>% as.integer())
  
  summary(all_split_read_combos_RBPs_w_scores_final)
  
  DBI::dbAppendTable(conn = con, name = "combo", value = all_split_read_combos_RBPs_w_scores_final)
  
  ## Disconnect the database
  DBI::dbDisconnect(conn = con)
}

 


