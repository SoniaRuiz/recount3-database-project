
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
                                  recount3.project.IDs,
                                  discard.minor.introns = F) {
  
  
  ## A) CREATE MASTER 'METADATA' TABLE -----------------------------------------
  
  SqlCreateMasterTableMetadata(database.path,
                               recount3.project.IDs,
                               results.folder)
  
  

  ## B) CREATE MASTER 'INTRON' TABLE -------------------------------------------

  
  SqlCreateMasterTableIntron(database.path,
                             gtf.version,
                             database.folder,
                             results.folder,
                             dependencies.folder,
                             discard.minor.introns = F)
  

  ## C) CREATE MASTER 'NOVEL' TABLE
  ## It contains novel 5' and 3' splicing events -------------------------------

  
  SqlCreateMasterTableNovel(database.path,
                            gtf.version,
                            database.folder,
                            results.folder,
                            dependencies.folder,
                            discard.minor.introns = F)
  
  

  ## D) CREATE MASTER 'COMBO' TABLE --------------------------------------------

  SqlCreateMasterTableCombo(database.path,
                            gtf.version,
                            database.folder,
                            results.folder,
                            dependencies.folder,
                            discard.minor.introns = F)
  
}



#' Title
#' Create metadata table
#' @param database.path Local path to the .sqlite database
#' @param recount3.project.IDs List of recount3 projects 
#' @param results.folder Path to the local folder where the results to read from are stored
#'
#' @return
#' @export
#'
#' @examples
SqlCreateMasterTableMetadata <- function(database.path,
                                         recount3.project.IDs,
                                         results.folder)  {
  
  logger::log_info("creating metadata table ... ")
  
  df_metadata <- map_df(recount3.project.IDs, function(project_id) {
    
    # project_id <- recount3.project.IDs[1]
    # project_id <- recount3.project.IDs[2]
    # project_id <- recount3.project.IDs[5]
    # project_id <- "TARDBP"
    
    logger::log_info("getting metadata info from ", project_id, "...")
    
    if (file.exists(paste0(results.folder, "/", project_id, "/base_data/", project_id, "_clusters_used.rds"))) {
      if (str_detect(database.path, pattern = "age")) {
        ## Age stratification
        metadata_file <- paste0(results.folder, "/", project_id, "/base_data/", project_id,"_age_samples_metadata.rds") 
      } else {
        metadata_file <- paste0(results.folder, "/", project_id, "/base_data/", project_id,"_samples_metadata.rds") 
      }
      
      if (file.exists(metadata_file)) {
        readRDS(file = metadata_file) %>%
          mutate(SRA_project = project_id) %>%
          return()
        
      } else {
        return(NULL)
      }
    }
    
  })
  
  
  ## Connect the database
  con <- dbConnect(RSQLite::SQLite(), database.path)
  
  
  # Generate the SQL string to create the db
  query <- paste0("CREATE TABLE 'metadata'", 
                  "(id INTEGER PRIMARY KEY, '",
                  paste(names(df_metadata), collapse = "','"), "')")

  res <- DBI::dbSendQuery(conn = con, statement = query)
  # DBI::dbRemoveTable(conn = con, name = "metadata")
  DBI::dbClearResult(res)
  
  
  ## Populate table ------------------------------------------------------
  DBI::dbAppendTable(conn = con, 
                     name = "metadata", 
                     value = df_metadata %>% tibble::rowid_to_column("id"))
  
  
  
  DBI::dbDisconnect(conn = con)
  logger::log_info(paste0("Table: 'metadata' created!"))
  
  
  
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
  
  
  ##########################################
  ## LOAD AND TIDY THE PAIR-WISE DISTANCES
  ##########################################
  
  if (file.exists(file.path(database.folder, "all_jxn_correct_pairings.rds"))) {
    logger::log_info("Loading the pre-generated data...")
    df_all_distances_pairings <- readRDS(file = file.path(database.folder, "all_jxn_correct_pairings.rds"))
    df_ambiguous_novel <- readRDS(file = file.path(database.folder, "all_jxn_ambiguous_pairings.rds"))
    df_introns_never <- readRDS(file = file.path(database.folder, "all_jxn_never_misspliced.rds"))
    df_introns_not_paired <- readRDS(file = file.path(database.folder, "all_jxn_not_paired.rds"))
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
    mutate(misspliced = "Yes") %>%
    dplyr::select(ref_junID, seqnames = ref_seq, start = ref_start, strand = ref_strand,
                  end = ref_end, gene_id, tx_id_junction, misspliced) %>%
    GRanges() %>% ## This is to get the width() as calculated by GRanges
    as_tibble()
  
  logger::log_info("Getting never mis-spliced introns...")
  
  ## Remove potential * in the junID of the reference introns
  if (any(str_detect(string = df_introns_never$ref_junID, pattern = "\\*"))) {
    stop("ERROR! some never mis-spliced junctions still have an *!")
  }
  
  df_introns_never_tidy <- df_introns_never %>% as_tibble() %>% mutate(misspliced = "No") 
  
  logger::log_info("getting the not paired mis-spliced introns...")
  
  ## Remove potential * in the junID of the reference introns
  if (any(str_detect(string = df_introns_never$ref_junID, pattern = "\\*"))) {
    stop("ERROR! some never mis-spliced junctions still have an *!")
  }
  
  df_introns_not_paired_tidy <- df_introns_not_paired %>% 
    as_tibble() %>%
    mutate(misspliced = "Potentially") #%>%
  
  if (intersect(df_introns_not_paired_tidy$ref_junID, df_all_introns$ref_junID) %>% length() > 0) {
    stop("ERROR! Non-paired introns classified as paired!")
  }
  
  df_introns_parenting_ambiguous <- df_ambiguous_novel %>% 
    as_tibble() %>%
    mutate(misspliced = "Potentially") %>%
    dplyr::filter(!(ref_junID %in% df_all_introns$ref_junID)) %>%
    distinct(ref_junID, .keep_all=T)
  
  ## QC
  
  if (intersect(df_introns_parenting_ambiguous$ref_junID, df_introns_not_paired_tidy$ref_junID) %>% length() > 0) {
    stop("ERROR! Introns parenting ambiguous junctions are classified as non-paired introns!")
  }
  
  if (intersect(df_introns_parenting_ambiguous$ref_junID, df_all_introns$ref_junID) %>% length() > 0) {
    stop("ERROR! Introns parenting ambiguous junctions are classified as paired introns!")
  }
  
  if (any(df_all_introns$width %>% abs() < 25)) {
    stop("ERROR! some mis-spliced introns are shorter than 25bp!")
  }
  
  if (any(df_introns_never_tidy$width %>% abs() < 25)) {
    stop("ERROR! some never mis-spliced introns are shorter than 25bp!")
  }
  
  if (any(df_introns_not_paired_tidy$width %>% abs() < 25)) {
    stop("ERROR! some POTENTIALLY mis-spliced introns are shorter than 25bp!")
  }
  
  if (any(df_introns_parenting_ambiguous$width %>% abs() < 25)) {
    stop("ERROR! some POTENTIALLY mis-spliced introns are shorter than 25bp!")
  }
  
  if (intersect(df_introns_never_tidy$ref_junID, df_all_introns$ref_junID) %>% length() > 0) {
    stop("ERROR! some never mis-spliced introns are mis-spliced!")
  }
  
  if (intersect(df_introns_not_paired_tidy$ref_junID, df_all_introns$ref_junID) %>% length() > 0) {
    stop("ERROR! some POTENTIALLY mis-spliced introns are CLASSIFIED AS mis-spliced!")
  }
  
  ## END QC
  
  ## Merge the different types of annotated introns
  df_introns_introverse <- rbind(df_introns_never_tidy, df_introns_not_paired_tidy, 
                                 df_introns_parenting_ambiguous, df_all_introns)
  
  if (any(str_detect(string = df_introns_introverse$ref_junID, pattern = "\\*" ))) {
    stop("ERROR! Some introns still have an ambiguous '*' strand!")
  }
  if (any(df_introns_introverse$ref_junID %>% duplicated())) {
    stop("ERROR! Duplicated IDs")
  }
  
  
  ######################################
  ## INTRONS - REMOVE AMBIGUOUS INTRONS
  ######################################
  
  logger::log_info("getting the ambiguous introns...")
  
  df_introns_introverse_tidy <- df_introns_introverse %>%
    distinct(ref_junID, .keep_all = T) %>% 
    rowwise() %>%
    mutate(ambiguous = ifelse(gene_id %>% unlist() %>% length() > 1, T, F))
  
  ## There should not be any ambiguous intron at this point
  if (any(df_introns_introverse_tidy %>% dplyr::filter(ambiguous == T) %>% nrow() > 0)) {
    stop("ERROR! Still there are some ambiguous introns")
  } else {
    df_introns_introverse_tidy <- df_introns_introverse_tidy %>% dplyr::select(-ambiguous)
  }
  
  logger::log_info(df_introns_introverse_tidy %>% distinct(ref_junID) %>% nrow(), " annotated introns TO BE STORED IN THE DATABASE...")
  
  

  ######################################
  ## GENES - CREATE GENE TABLE
  ######################################
  
  gene_ids <- df_introns_introverse_tidy %>% unnest(gene_id) %>% distinct(gene_id)
  
  logger::log_info("loading GRCh38 reference...")
  
  hg38 <- rtracklayer::import(con = paste0(dependencies.folder, "/Homo_sapiens.GRCh38.",gtf.version,".chr.gtf")) %>% 
    as_tibble() %>%
    mutate(gene_id = str_sub(gene_id, start = 1, end = 15)) %>%
    dplyr::filter(gene_id %in% (gene_ids$gene_id))
  
  
  SqlCreateMasterTableGene(database.path = database.path,
                           hg38 = hg38,
                           gene.ids = df_introns_introverse_tidy %>% unnest(gene_id) %>% distinct(gene_id))
  
  
  ######################################
  ## TX_JUNCTION - CREATE TX TABLE
  ######################################
  
  SqlCreateMasterTableTranscript(database.path = database.path,
                                 gene.ids = df_introns_introverse_tidy %>% unnest(gene_id) %>% distinct(gene_id),
                                 dependencies.folder = dependencies.folder,
                                 hg38 = hg38,
                                 tx.ids = df_introns_introverse_tidy %>% unnest(tx_id_junction) %>% distinct(tx_id_junction))
  rm(hg38)
  gc()

  
  ######################################
  ## INTRONS - ADD MAXENTSCAN INFO 
  ######################################
  
  logger::log_info("Adding the MaxEntScan info ...")
  
  wd <- getwd()
  if ( !file.exists(file.path(dependencies.folder, "Homo_sapiens.GRCh38.dna.primary_assembly.fa")) ) {
    stop("ERROR! File dependency 'Homo_sapiens.GRCh38.dna.primary_assembly.fa' does not exist within the specified dependencies folder.")
  }
  
  ## Add MaxEntScan score to the split reads
  all_split_reads_tidy <- GenerateMaxEntScore(db.introns = df_introns_introverse_tidy %>% dplyr::rename(junID = ref_junID) %>% distinct(junID, .keep_all = T),
                                              max.ent.tool.path = paste0(dependencies.folder, "/fordownload/"),
                                              bedtools.path = file.path(dependencies.folder, "bedtools2/"),
                                              hs.fasta.path = paste0(dependencies.folder, "/Homo_sapiens.GRCh38.dna.primary_assembly.fa"))
  rm(df_introns_introverse_tidy)
  gc()
  
  all_split_reads_tidy <- all_split_reads_tidy %>% as_tibble() %>% 
    dplyr::select(-c(donorSeqStart, donorSeqStop, AcceptorSeqStart, AcceptorSeqStop)) %>%
    dplyr::rename(ref_donor_sequence = donor_sequence, ref_acceptor_sequence = acceptor_sequence)
  
  setwd(wd)
  
  
  df_all_introns <- all_split_reads_tidy %>%
    mutate(ref_junID = paste0("chr", seqnames, ":", start, "-", end, ":", strand)) 
  
  if ((setdiff(df_all_introns$junID, df_all_introns$ref_junID) %>% length()) > 0) {
    stop("ERROR! ")
  }
  
  df_all_introns_gr <- df_all_introns %>%
    dplyr::select(-one_of("junID", "ref_ss5score", "ref_ss3score")) %>% 
    dplyr::rename(ref_mes5ss = ss5score, ref_mes3ss = ss3score) %>%
    dplyr::relocate(ref_junID) %>%
    dplyr::relocate(c(ref_mes5ss, ref_mes3ss), .after = tx_id_junction ) %>% 
    distinct(ref_junID, .keep_all = T) %>% 
    dplyr::rename("junID" = "ref_junID") %>% 
    GRanges()
  
  ################################################
  ## INTRONS - ADD THE CONSERVATION AND CDTS INFO
  ################################################
  
  logger::log_info("adding CDTS and Conservation scores...")
  
  seqlevelsStyle(df_all_introns_gr) <- "UCSC"
  df_all_introns_gr <- GenerateCdtsPhastconsScores(dependencies.folder = dependencies.folder,
                                                   db.introns = df_all_introns_gr,
                                                   intron.size = 100, 
                                                   phastcons.type = 17)
  
  df_all_introns <- df_all_introns_gr %>% as_tibble() %>% dplyr::rename("ref_junID" = "junID") %>% mutate_if(is.numeric, ~replace_na(., 0))
  df_all_introns
  
  rm(df_all_introns_gr)
  gc()
  
  ######################################
  ## INTRONS - ADD CLINVAR DATA
  ######################################
  
  logger::log_info("adding the ClinVar data...")
  
  df_all_introns_gr <- AddClinvarData(db.introns = df_all_introns, dependencies.folder)
  
  df_all_introns_gr %>% head()
  df_all_introns_gr %>% as_tibble() %>% filter(clinvar == T) %>% dplyr::select(ref_junID, clinvar)
  
  rm(df_all_introns)
  gc()
  
  
  ######################################
  ## INTRONS - ADD THE INTRON TYPE
  ## This intron type corresponds to whether the intron is 
  ## spliced out by the minor or the major spliceosome
  ######################################
  
  logger::log_info("adding the IAOD intron data...")
  
  ## Load intron type files
  u12_introns_gr <- read.table(file = file.path(dependencies.folder, "GRCh38_U12.bed"), 
                               header = FALSE, sep="\t", stringsAsFactors=FALSE, quote="") %>% GRanges() 
  seqlevelsStyle(u12_introns_gr) <- "UCSC"
  
  ## Add a new column to incorporate info about intron type
  elementMetadata(df_all_introns_gr)[, "u2_intron"] = T
  
  ## MINOR INTRON
  logger::log_info("Getting junctions spliced out by the minor spliceosome.")
  overlaps <- GenomicRanges::findOverlaps(query = u12_introns_gr,
                                          subject = df_all_introns_gr,
                                          ignore.strand = FALSE, 
                                          type = "equal")
  
  logger::log_info(queryHits(overlaps) %>% length(), " introns spliced by the minor spliceosome!")
  df_all_introns_gr[subjectHits(overlaps),]$u2_intron <- F
  
  if (which(df_all_introns_gr$u2_intron == F) %>% length() != df_all_introns_gr[subjectHits(overlaps),] %>% length()) {
    stop("ERROR! Minor spliceosomal introns have not been added adequately!")
  }
  
  ## Discard introns spliced out by the minor spliceosome. Comment out to avoid applying this filter.
  df_all_introns_tidy <- df_all_introns_gr %>%
    as_tibble()
  
  if (discard.minor.introns) {
    df_all_introns_tidy <- df_all_introns_tidy %>% filter(u2_intron == T) %>% dplyr::select(-u2_intron)
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
  ## INTRON MASTER TABLE 
  ## CREATE AND POPULATE THE TABLE
  ######################################
  
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
                  
                  misspliced TEXT NOT NULL)")
  
  
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
  
  DBI::dbAppendTable(conn = con, name = "intron", value = df_all_introns_tidy_final %>% dplyr::select(-tx_id_junction))
  
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
  
  
  ############################################
  ## CREATE BRIDGE TRANSCRIPT TABLE N:N
  ############################################
  
  logger::log_info("Creating bridge table between the master tables INTRON and TRANSCRIPT...")
  
  
  SqlCreateBridgeTablewTranscript(db.intron = df_all_introns_tidy_final,
                                  bridge.table.name = "bridge_intron_transcript",
                                  origin.master.table = "intron",
                                  database.path)
  
  
  logger::log_info(df_all_introns_tidy_final %>% distinct(ref_junID) %>% nrow(), " novel combos stored!")
  summary(df_all_introns_tidy_final %>% dplyr::select(-tx_id_junction))
  
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
  all_split_reads_tidy <- GenerateMaxEntScore(db.introns = df_all_novel_raw_tidy %>% dplyr::rename(junID = novel_junID) %>% distinct(junID, .keep_all = T),
                                              max.ent.tool.path = file.path(dependencies.folder,"fordownload/"),
                                              bedtools.path = file.path(dependencies.folder, "bedtools2/"),
                                              hs.fasta.path = file.path(dependencies.folder, "Homo_sapiens.GRCh38.dna.primary_assembly.fa"))

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
  
  df_all_novels_gr <- df_all_novels_tidy %>%
    dplyr::select(-c(any_of(c("junID","novel_ss5score","novel_ss3score")))) %>% 
    dplyr::rename(novel_mes5ss = ss5score, novel_mes3ss = ss3score) %>%
    dplyr::relocate(ref_junID, novel_junID) %>%
    dplyr::relocate(c(novel_mes5ss, novel_mes3ss), .before = novel_type ) %>%
    distinct(novel_junID, .keep_all = T) %>% 
    dplyr::rename(junID = "novel_junID") %>% 
    GRanges()
  
  ################################################
  ## NOVEL - ADD THE CONSERVATION AND CDTS INFO
  ################################################
  
  
  logger::log_info("\t Adding CDTS and Conservation scores to the novel junctions ...")
  
  
  seqlevelsStyle(df_all_novels_gr) <- "UCSC"
  df_all_novels_tidy <- GenerateCdtsPhastconsScores(dependencies.folder = dependencies.folder,
                                                    db.introns = df_all_novels_gr,
                                                    intron.size = 100,
                                                    phastcons.type = 17) %>% 
    as_tibble() %>% 
    dplyr::rename("novel_junID" = "junID") %>% 
    mutate_if(is.numeric, ~replace_na(., 0))
  
  df_all_novels_tidy %>% as_tibble()
  
  rm(df_all_novels_gr)
  
  ######################################
  ## NOVEL - ADD CLINVAR DATA
  ######################################
  
  logger::log_info("adding the ClinVar data...")
  
  df_all_novels_tidy_gr <- AddClinvarData(db.introns = df_all_novels_tidy, dependencies.folder)
  
  
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
  
  logger::log_info("SQL connection stablished!") 
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
  
  
  if (is.null(recount3.project.IDs)){
    recount3.project.IDs <- (df_metadata$SRA_project %>% unique())
  }
  
  ## Get all novel combos across all tables
  all_split_reads_combos <- map_df(recount3.project.IDs, function(project_id) { 
    
    # project_id <- recount3.project.IDs[1]
    
    logger::log_info("Working with '", project_id, "' ...")
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
  

  
        
  ######################################
  ## ADD MAXENTSCAN INFO 
  ######################################
  
  logger::log_info("Adding the MaxEntScan info ...")
  
  wd <- getwd()
  if ( !file.exists(file.path(dependencies.folder, "/Homo_sapiens.GRCh38.dna.primary_assembly.fa")) ) {
    stop("ERROR! File dependency 'Homo_sapiens.GRCh38.dna.primary_assembly.fa' does not exist within the specified dependencies folder.")
  }
  
  ## Add MaxEntScan score to the split reads
  all_split_read_combos_RBPs_w_MES <- GenerateMaxEntScore(db.introns = all_split_read_combos_RBPs,
                                                          max.ent.tool.path = file.path(dependencies.folder, "fordownload/"),
                                                          bedtools.path = file.path(dependencies.folder, "bedtools2/"),
                                                          hs.fasta.path = file.path(dependencies.folder, "Homo_sapiens.GRCh38.dna.primary_assembly.fa")) %>% as_tibble()
  
  all_split_read_combos_RBPs_w_MES <- all_split_read_combos_RBPs_w_MES %>% 
    dplyr::select(-c(donorSeqStart, donorSeqStop, AcceptorSeqStart, AcceptorSeqStop)) %>%
    dplyr::rename(ref_donor_sequence = donor_sequence, ref_acceptor_sequence = acceptor_sequence) %>%
    dplyr::rename(ref_mes5ss = ss5score, ref_mes3ss = ss3score) %>%
    dplyr::relocate(c(ref_mes5ss, ref_mes3ss), .before = tx_id_junction) %>% 
    as_tibble()
  
  setwd(wd)
  
  ################################################
  ## ADD THE CONSERVATION AND CDTS INFO
  ################################################
  
  all_split_read_combos_RBPs_w_MES <- all_split_read_combos_RBPs_w_MES %>% GRanges()

  logger::log_info("adding CDTS and Conservation scores...")
  
  seqlevelsStyle(all_split_read_combos_RBPs_w_MES) <- "UCSC"
  all_split_read_combos_RBPs_w_scores <- GenerateCdtsPhastconsScores(db.introns = all_split_read_combos_RBPs_w_MES,
                                                                     intron.size = 100,
                                                                     phastcons.type = 17,
                                                                     dependencies.folder = dependencies.folder) %>% 
    as_tibble() %>% 
    mutate_if(is.numeric, ~replace_na(., 0))
  
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
  
  logger::log_info("Creating 'novel combo' table ... ")
  
  
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
                  
                  left_motif TEXT NOT NULL,
                  right_motif TEXT NOT NULL,
                  ref_donor_sequence TEXT NOT NULL,
                  ref_acceptor_sequence TEXT NOT NULL,
                  
                  protein_coding DOUBLE NOT NULL)")
        

  ## Connect the database 
  con <- dbConnect(RSQLite::SQLite(), database.path)
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
        
  ## Create the NOVEL COMBO table
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
        
  # DBI::dbRemoveTable(conn = con, name = "combo")
  
  
  ############################################
  ## POPULATE THE NOVEL COMBO TABLE
  ############################################
  
  all_split_read_combos_RBPs_w_scores_final <- all_split_read_combos_RBPs_w_scores %>%
    dplyr::rename(ref_coordinates = junID, ref_length = width) %>%
    tibble::rowid_to_column("ref_junID") %>%
    dplyr::select(-gene_id)
  
  summary(all_split_read_combos_RBPs_w_scores_final)
  
  
  DBI::dbAppendTable(conn = con, name = "combo", value = all_split_read_combos_RBPs_w_scores_final %>% dplyr::select(-tx_id_junction))
  
  ## Disconnect the database
  DBI::dbDisconnect(conn = con)
  
  ############################################
  ## CREATE BRIDGE TRANSCRIPT TABLE N:N
  ############################################
  
  logger::log_info("\t Adding TRANSCRIPT foreing key to the novel combos...")
  
  SqlCreateBridgeTablewTranscript(db.intron = all_split_read_combos_RBPs_w_scores_final,
                                  bridge.table.name = "bridge_combo_transcript",
                                  origin.master.table = "combo",
                                  database.path)
  
  
  logger::log_info(all_split_read_combos_RBPs_w_scores_final %>% distinct(ref_junID) %>% nrow(), " novel combos stored!")
  summary(all_split_read_combos_RBPs_w_scores_final %>% dplyr::select(-tx_id_junction))
 
}

 

#' Title
#' Creates the master table 'gene'
#' @param database.path Local path to the .sqlite databse
#' @param hg38 Human genome. Needed to extract info about the transcripts to be stored.
#' @param gene.ids List of Ensembl gene IDs
#'
#' @return
#' @export
#'
#' @examples
SqlCreateMasterTableGene <- function(database.path, hg38,
                                     gene.ids) {
  
  
  
  logger::log_info("Creating 'gene' master table...")
  con <- dbConnect(RSQLite::SQLite(), database.path)
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  
  
  hg38_transcripts_gene <- hg38 %>%
    dplyr::count(gene_id, type) %>%
    dplyr::filter(type == "transcript") %>%
    unnest(gene_id) %>%
    dplyr::select(-type) %>%
    dplyr::rename(n_transcripts = n)
  
  
  hg38_genes <- hg38 %>%
    dplyr::filter(type == "gene") %>%
    dplyr::select(gene_id, gene_name, gene_width = width, gene_biotype)
  
  
  hg38_tidy <- hg38_transcripts_gene %>%
    inner_join(hg38_genes, by = "gene_id" ) %>% 
    as_tibble()
  
  
  
  ## CREATE GENE_NAME TABLE
  sql_statement <- paste0("CREATE TABLE IF NOT EXISTS 'gene'", 
                          "(id INTEGER PRIMARY KEY NOT NULL,
                          gene_id TEXT NOT NULL,
                          gene_name TEXT,
                          gene_width INTEGER NOT NULL, 
                          n_transcripts INTEGER NOT NULL,
                          gene_biotype TEXT NOT NULL)")
  res <- DBI::dbSendQuery(con, sql_statement)
  DBI::dbClearResult(res)
  
  ## POPULATE GENE TABLE
  hg38_tidy_final <- hg38_tidy %>% as_tibble() %>% tibble::rowid_to_column("id")
  
  DBI::dbAppendTable(conn = con,
                     name = "gene", 
                     value = hg38_tidy_final)
  
  
  logger::log_info("Gene table created!")
  logger::log_info(hg38_tidy_final %>% nrow(), " genes stored!")
  
  
  
  ## CREATE INDEXES ---------------------------------------------------
  
  query <- paste0("CREATE UNIQUE INDEX 'index_gene_id' ON 'gene'(id)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  query <- paste0("CREATE UNIQUE INDEX 'index_gene_ensembl_id' ON 'gene'(gene_id,gene_name)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  DBI::dbDisconnect(con)
}


#' Title
#' Create master table 'transcript'
#' @param database.path Local path to the .sqlite database
#' @param gene.ids List of gene Ensembl IDs
#' @param hg38 Human genome. Needed to extract info about the transcripts to be stored.
#' @param tx.ids List of transcript Ensembl IDs
#'
#' @return
#' @export
#'
#' @examples
SqlCreateMasterTableTranscript  <- function(database.path,
                                            gene.ids,
                                            dependencies.folder,
                                            hg38,
                                            tx.ids) {
  
  ## Connect database ----------------------------------------------------------
  
  logger::log_info("Creating 'transcript' master table...")
  con <- dbConnect(RSQLite::SQLite(), database.path)
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  
  
  
  ## Load MANE transcripts -----------------------------------------------------
  
  hg_mane_transcripts <- rtracklayer::import(con = file.path(dependencies.folder, "MANE.GRCh38.v1.4.ensembl_genomic.gtf")) %>%
    as_tibble() %>%
    dplyr::select(-c(source, score, phase, gene_id, gene_type, tag, protein_id,
                     db_xref, transcript_type, exon_id, exon_number, width)) %>%
    mutate(transcript_id = transcript_id %>% str_sub(start = 1, end = 15)) %>%
    drop_na() %>%
    dplyr::filter(type == "transcript") %>%
    distinct(transcript_id) %>%
    mutate(MANE = T)
  

  
  ## GET TRANSCRIPT ID
  hg38_transcripts_gene <- hg38 %>%
    as_tibble() %>%
    dplyr::filter(transcript_id %in% tx.ids$tx_id_junction, type == "transcript") %>%
    dplyr::select(transcript_id, TSL = transcript_support_level, gene_id, transcript_biotype) %>%
    mutate(TSL = str_sub(TSL, start = 1, end = 2) %>% as.integer())
  
  
  hg38_transcripts_gene[is.na(hg38_transcripts_gene[,"TSL"]),"TSL"] <- 10
  
  ## ADD MANE INFO
  hg38_transcripts_gene_mane <- hg38_transcripts_gene %>%
    left_join(y = hg_mane_transcripts, by = "transcript_id" )
  hg38_transcripts_gene_mane[is.na(hg38_transcripts_gene_mane[,"MANE"]),"MANE"] <- F
  
  ## ADD THE GENE FOREING KEY
  logger::log_info("Adding GENE foreing key to the transcripts...")
  
  query <- paste0("SELECT id, gene_id FROM 'gene'")
  df_genes <- dbGetQuery(con, query) %>% as_tibble()
  
  ## Add the GENE ID for the foreign key
  hg38_transcripts_gene_mane <- hg38_transcripts_gene_mane %>%
    inner_join(df_genes %>% dplyr::select(id, gene_id), by = "gene_id") %>%
    dplyr::select(-gene_id) %>%
    dplyr::rename(gene_id = id) 
  
  ## Add the Id to the table
  hg38_transcripts_final <- hg38_transcripts_gene_mane %>% 
    as_tibble() %>%
    tibble::rowid_to_column("id") 
  
  
  ####################################
  ## CREATE TRANSCRIPT TABLE AND POPULATE
  ####################################
  
  query <- paste0("CREATE TABLE IF NOT EXISTS 'transcript'",
                  "(id NUMERIC NOT NULL,
                  transcript_id TEXT NOT NULL,
                  TSL NUMERIC NOT NULL, 
                  MANE BOOL NOT NULL, 
                  transcript_biotype TEXT NOT NULL, 
                  gene_id INTEGER NOT NULL,
                  PRIMARY KEY (id),
                  FOREIGN KEY (gene_id) REFERENCES 'gene'(id))")
  
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  logger::log_info("'Transcript' table created!")
  
  if ( any(duplicated(hg38_transcripts_final$transcript_id)) ) {
    logger::log_info("ERROR! some novel junctions are duplicated")
  }
  DBI::dbAppendTable(conn = con,
                     name = "transcript", 
                     value = hg38_transcripts_final)
  
  logger::log_info("'Transcript' table populated!") 
  logger::log_info(hg38_transcripts_final %>% nrow(), " transcripts stored!")
  
  
  ## CREATE INDEXES ---------------------------------------------------
  
  query <- paste0("CREATE UNIQUE INDEX 'index_transcript_id' ON 'transcript'(id)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  query <- paste0("CREATE UNIQUE INDEX 'index_transcript_ensembl_id' ON 'transcript'(transcript_id)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  
  DBI::dbDisconnect(con)
}



#' Title
#' Function to create a bridge table between the master tables 'intron' and 'transcript' due to their
#' N:N relationship
#' @param db.intron 
#' @param database.path 
#'
#' @return
#' @export
#'
#' @examples
SqlCreateBridgeTablewTranscript <- function(db.intron,
                                            database.path,
                                            origin.master.table,
                                            bridge.table.name) {
  
  
  
  ## Connect to the databse ----------------------------------------------
  con <- dbConnect(RSQLite::SQLite(), database.path)
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  
  
  
  query <- paste0("SELECT id, transcript_id FROM 'transcript'")
  df_transcripts <- dbGetQuery(con, query) %>% as_tibble()
  
  
  ## Add the GENE ID for the foreign key
  db_introns_tidy <- db.intron %>%
    unnest(tx_id_junction) %>% 
    left_join(df_transcripts, 
              by = c("tx_id_junction" = "transcript_id")) %>%
    dplyr::select(any_of(c("ref_junID", "junID")), transcript_id = id)
  
  rm(df_transcripts)
  
 
  ## Create bridge table -------------------------------------------------
  query <- paste0("CREATE TABLE ", bridge.table.name, "(
    id INTEGER PRIMARY KEY,
    intron_id INTEGER,
    transcript_id INTEGER,
    FOREIGN KEY (intron_id) REFERENCES ", origin.master.table, "(ref_junID),
    FOREIGN KEY (transcript_id) REFERENCES transcript(id));")
  
  
  
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  logger::log_info("Bridge '", bridge.table.name, "' table created!")
  
  # DBI::dbRemoveTable(con, name = bridge.table.name)
  
  
  ## Populate table ------------------------------------------------------
  db_introns_final <- db_introns_tidy %>% 
    tibble::rowid_to_column("id") %>% 
    dplyr::rename(intron_id = any_of(c("ref_junID", "junID")))
  
  DBI::dbAppendTable(conn = con,
                     name = bridge.table.name, 
                     value = db_introns_final)
  
  logger::log_info("Bridge '", bridge.table.name, "' table populated!")
  
  return(db_introns_tidy)
  
}
