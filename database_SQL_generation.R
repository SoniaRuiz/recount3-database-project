
##########################################
## LOAD DATA DEPENDENCIES
##########################################

if (!exists("CNC_CDTS_CONS_gr")) {
  print("Loading the 'CNC_CDTS_CONS_gr' file...")
  load( file = paste0(dependencies_folder, "/CNC_CDTS_CONS_gr.rda") )
  print("'CNC_CDTS_CONS_gr' file loaded!")
} else {
  print("'CNC_CDTS_CONS_gr' file already loaded!")
}

## GET INFO FROM MANE
if ( !exists("hg_mane_transcripts") ) {

  print("Loading the 'hg_MANE' file...")
  hg_MANE <- rtracklayer::import(con = paste0(dependencies_folder,
                                              "/MANE.GRCh38.v1.0.ensembl_genomic.gtf"))
  hg_MANE_tidy <- hg_MANE %>%
    as_tibble() %>%
    dplyr::select(-source, -score, -phase, -gene_id, -gene_type, -tag, -protein_id,
                  -db_xref,-transcript_type,-exon_id,-exon_number, -width ) %>%
    mutate(transcript_id = transcript_id %>% str_sub(start = 1, end = 15)) %>%
    drop_na()

  hg_mane_transcripts <- hg_MANE_tidy %>%
    dplyr::filter(type == "transcript") %>%
    distinct(transcript_id) %>%
    mutate(MANE = T)

} else {
  print("'hg_mane_transcripts' file already loaded!")
}

# GET INFO FROM THE HG38
if ( !exists("hg38_transcripts") ) {

  print("Loading the 'hg38' file...")

  if ( !exists("hg38") ) {
    hg38 <- rtracklayer::import(con = paste0(dependencies_folder,
                                             "/Homo_sapiens.GRCh38.105.chr.gtf"))
  }

  hg38_transcripts <- hg38 %>%
    as_tibble() %>%
    dplyr::filter(type == "transcript") %>%
    dplyr::select(transcript_id, transcript_support_level) %>%
    mutate(transcript_support_level = str_sub(string = transcript_support_level,
                                              start = 1,
                                              end = 1))
} else {
  print("'hg38_transcripts' file already loaded!")
}


## Remove the tables -----------------------------------------------------------

remove_tables <- function(database.path,
                          all) {
  
  con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = database.path)
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=0")
  tables <- dbListTables(con)
  tables %>% print()
  
  for (table in tables) {
    #if (str_detect(table %>% tolower(), pattern = "kidney")) {
    if (!all) {
      if (!(table %in% c("gene", "transcript", "intron", 
                         "master", "novel"))) {
        dbRemoveTable(conn = con, table)
        print(paste0("Table: '", table, "' has been removed!"))
      }
    } else {
      dbRemoveTable(conn = con, table)
      print(paste0("Table: '", table, "' has been removed!"))
    }
  }
  dbListTables(con)
}
#dbRemoveTable(conn = con, "PD_SRP058181_misspliced")


## METADATA table --------------------------------------------------------------

create_metadata_table <- function(database.path,
                                  recount3.project.IDs,
                                  gtf.version,
                                  project.name,
                                  results.folder)  {
  
  print(paste0(Sys.time(), " - creating metadata table ... "))
  
  df_metadata <- map_df(recount3.project.IDs, function(project_id) { 
    
    print(paste0(Sys.time(), " - getting info from ", project_id))
    # project_id <- recount3.project.IDs[1]
    # project_id <- recount3.project.IDs[2]
    # project_id <- recount3.project.IDs[6]
    
    metadata_file <- paste0(results.folder, "/", project_id, 
                            "/base_data/",project_id,"_samples_metadata.rds")
    
    if ( file.exists(metadata_file) ) {
      
      metadata <- readRDS(file = metadata_file)
      
      #metadata <- metadata %>% 
      #  as_tibble() %>%
      #  filter(gtex.smafrze != "EXCLUDE") %>%
      #  distinct(external_id, .keep_all = T) 
      
    #  min_samples <- 1
      
    #  if (metadata %>% nrow() >= min_samples) {
        
    #    
    #    print(metadata$gtex.smtsd %>% unique())
    #    print(metadata %>% distinct(gtex.sampid) %>% nrow())
     #   
     #   df_project_metadata <- data.frame(age = metadata$gtex.age %>% as.character(),
      #                                    rin = metadata$gtex.smrin %>% as.character(),
       #                                   gender = metadata$gtex.sex %>% as.character(),
        #                                  cluster = metadata$gtex.smtsd, 
         #                                 smnabtcht = metadata$gtex.smnabtcht, 
          #                                sample_id = metadata %>% distinct(gtex.sampid) %>% nrow(),
           #                               smafrze = metadata$gtex.smafrze,
            #                              avg_read_length = metadata$recount_seq_qc.avg_len,
             #                             mapped_read_count = metadata$recount_qc.star.all_mapped_reads,
              #                            SRA_project = metadata$recount_project.project) %>% 
          #arrange(SRA_project, cluster) %>% 
          #return()
        
    #  }else {
        
    #    return(NULL)
    #  }
      return(metadata %>%
               mutate(SRA_project = project_id))
    } else {
      
      return(NULL)
    }
    
  })
  
  
  df_metadata %>% as_tibble() %>% dplyr::count(cluster)
  
  
  
  con <- dbConnect(RSQLite::SQLite(), database.path)
  DBI::dbWriteTable(conn = con,
                    name = "master",
                    value = df_metadata %>% as_tibble(),
                    overwrite = T)
  
  DBI::dbDisconnect(conn = con)
  
  print(paste0("Table: 'metadata' created!"))
}


## Create the master tables ----------------------------------------------------

create_master_tables <- function(database.path,
                                 project.name,
                                 gtf.version,
                                 database.folder,
                                 results.folder) {
  
  
  
  
  ##########################################
  ## LOAD AND TIDY THE PAIR-WISE DISTANCES
  ##########################################
  
  con <- dbConnect(RSQLite::SQLite(), database.path)

  
  if ( file.exists(paste0(database.folder, "/all_jxn_correct_pairings.rds")) ) {
    
    print(paste0(Sys.time(), " - loading the pre-generated pair-wise distances data..."))
    
    df_all_distances_pairings_raw <- readRDS(file = paste0(database.folder, "/all_jxn_correct_pairings.rds"))
    
    df_ambiguous_novel <- readRDS(file = paste0(database.folder, "/all_jxn_ambiguous_pairings.rds"))
    
    df_introns_never <- readRDS(file = paste0(database.folder, "/all_jxn_never_misspliced.rds"))
    
  } else {
    
    print("ERROR loading dependencies!")
    break;
    
  }
  
  ####################################
  ## A) GET ALL INTRONS
  ####################################
  
  
  ## Get all annotated introns
  print(paste0(Sys.time(), " - getting the mis-spliced introns..."))
  
  df_all_introns <- df_all_distances_pairings_raw %>%
    distinct(ref_junID, .keep_all = T) %>%
    as_tibble() %>%
    mutate(misspliced = T) %>%
    dplyr::select(ref_junID,seqnames = ref_seq,
                  start = ref_start,
                  strand = ref_strand,
                  end = ref_end,
                  gene_id, 
                  tx_id_junction,
                  misspliced) %>%
    GRanges() %>% ## This is to get the width() as calculated by GRanges
    as_tibble()
  
  
  print(paste0(Sys.time(), " - getting the never mis-spliced introns..."))
  
  ## Remove potential * in the junID of the reference introns
  if (any(str_detect(string = df_introns_never$ref_junID, pattern = "\\*"))) {
    print("ERROR! some never mis-spliced junctions still have an *!")
    break;
  }
  
  df_introns_never_tidy <- df_introns_never %>% 
    as_tibble() %>%
    mutate(misspliced = F) %>%
    dplyr::filter(!(ref_junID %in% df_all_introns$ref_junID)) %>%
    dplyr::filter(!(ref_junID %in% df_ambiguous_novel$ref_junID)) 
  
  
  if (any(df_introns_never_tidy$width %>% abs() < 25)) {
    print("ERROR! some mis-spliced introns are shorter than 25bp!")
  }
  if (any(df_introns_never_tidy$width %>% abs() < 25)) {
    print("ERROR! some never mis-spliced introns are shorter than 25bp!")
  }
  if (intersect(df_introns_never_tidy$ref_junID, df_all_introns$ref_junID) %>% length() > 0) {
    print("ERROR! some never mis-spliced introns are mis-spliced!")
  }
  
  
  ## Merge mis-spliced introns and never mis-spliced introns
  df_introns_introverse <- rbind(df_introns_never_tidy, df_all_introns)
  if (any(str_detect(string = df_introns_introverse$ref_junID, pattern = "\\*" ))) {
    print("ERROR! IDs with an star!")
  }
  if (any(df_introns_introverse$ref_junID %>% duplicated())) {
    print("ERROR! Duplicated IDs")
  }
  
  #saveRDS(object = df_introns_introverse %>%
  #          distinct(ref_junID, .keep_all = T),
  #        file = paste0(folder_database, "/df_all_introns_database.rds"))
  
  
  
  ######################################
  ## INTRONS - REMOVE AMBIGUOUS INTRONS
  ######################################
  
  print(paste0(Sys.time(), " - getting the ambiguous introns..."))
  
  # df_introns_introverse <- df_introns_introverse %>%
  #   dplyr::filter(!(ref_junID %in% df_ambiguous_intron$junID)) %>%
  #   as_tibble()
  
  df_introns_introverse_tidy <- df_introns_introverse %>%
    distinct(ref_junID, .keep_all = T) %>% 
    rowwise() %>%
    mutate(ambiguous = ifelse(gene_id %>% unlist() %>% length() > 1, T, F))
  
  ## There should not be any ambiguous intron at this point
  if ( any(df_introns_introverse_tidy %>% dplyr::filter(ambiguous == T) %>% nrow() > 0) ) {
    print("ERROR! Still there are some ambiguous introns")
    break;
  } else {
    df_introns_introverse_tidy <- df_introns_introverse_tidy %>%
      dplyr::select(-ambiguous)
  }
  
  #saveRDS(object = df_introns_introverse_tidy,
  #        file = paste0(folder_database, "/df_all_introns_database_tidy.rds"))
  
  ######################################
  ## GENES - CREATE GENE TABLE
  ######################################
  
  create_gene_table(database.path = database.path,
                    gene_ids = df_introns_introverse_tidy %>% unnest(gene_id) %>% distinct(gene_id) )
  
  
  ######################################
  ## TX_JUNCTION - CREATE TX TABLE
  ######################################

    
  create_transcript_table(database.path = database.path,
                          gene_ids = df_introns_introverse_tidy %>% unnest(gene_id) %>% distinct(gene_id),
                          tx_ids = df_introns_introverse_tidy %>% unnest(tx_id_junction) %>% distinct(tx_id_junction))
  
  
  ######################################
  ## INTRONS - ADD THE TRANSCRIPT FOREING KEY
  ######################################
  
  print(paste0(Sys.time(), " --> adding TRANSCRIPT foreing key to the introns..."))
  
  query <- paste0("SELECT id, transcript_id FROM 'transcript'")
  df_transcripts <- dbGetQuery(con, query) %>% as_tibble()
  
  
  ## Add the GENE ID for the foreign key
  df_introns_introverse_tidy <- df_introns_introverse_tidy %>%
    unnest(tx_id_junction) %>%
    inner_join(df_transcripts,
               by = c("tx_id_junction" = "transcript_id")) %>%
    dplyr::select(-gene_id, -tx_id_junction) %>%
    dplyr::rename(transcript_id = id) 
  
  
  df_introns_introverse_tidy %>%
    distinct(ref_junID)
  
 
  ######################################
  ## INTRONS - ADD MAXENTSCAN INFO 
  ######################################
  
  print(paste0(Sys.time(), " --> adding the MaxEntScan info ..."))
  
  wd <- getwd()
  if ( !file.exists(paste0(dependencies_folder, 
                           "/Homo_sapiens.GRCh38.dna.primary_assembly.fa")) ) {
    print(paste0("ERROR! File dependency 'Homo_sapiens.GRCh38.dna.primary_assembly.fa'
                 does not exist within the specified dependencies folder."))
  }
  ## Add MaxEntScan score to the split reads
  all_split_reads_tidy <- generate_max_ent_score(junc_tidy = df_introns_introverse_tidy %>% dplyr::rename(junID = ref_junID),
                                                 max_ent_tool_path = paste0(dependencies_folder, "/fordownload/"),
                                                 homo_sapiens_fasta_path = paste0(dependencies_folder, 
                                                                                  "/Homo_sapiens.GRCh38.dna.primary_assembly.fa") )
  
  all_split_reads_tidy <- all_split_reads_tidy %>% as_tibble()
  
  all_split_reads_tidy <- all_split_reads_tidy %>% 
    dplyr::select(-donorSeqStart,
                  -donorSeqStop,
                  -AcceptorSeqStart,
                  -AcceptorSeqStop,
                  -donor_sequence,
                  -acceptor_sequence)
  
  setwd(wd)
  
  df_all_introns <- all_split_reads_tidy %>%
    mutate(ref_junID = paste0("chr", seqnames, ":", start, "-", end, ":", strand)) 
  
  if ((setdiff(df_all_introns$junID, df_all_introns$ref_junID) %>% length()) > 0) {
    print("ERROR! ")
    break;
  }
  
  df_all_introns %>% as_tibble()
  
  df_all_introns <- df_all_introns %>%
    dplyr::select(-one_of("junID","ref_ss5score","ref_ss3score")) %>% 
    dplyr::rename(ref_ss5score = ss5score, 
                  ref_ss3score = ss3score) %>%
    dplyr::relocate(ref_junID) %>%
    dplyr::relocate(c(ref_ss5score, ref_ss3score), .before = transcript_id ) 
  
  
  df_all_introns %>% as_tibble()
  
  
  ################################################
  ## INTRONS - ADD THE CONSERVATION AND CDTS INFO
  ################################################
  
  df_all_introns %>% as_tibble()
  
  print(paste0(Sys.time(), " - adding CDTS and Conservation scores..."))
  
  df_all_introns <- add_cdts_cons_scores(db_introns = df_all_introns %>% as_tibble()) %>%
    as_tibble()
  
  df_all_introns <- df_all_introns %>%
    dplyr::rename(ref_CDTS5score = CDTS_5ss_mean,
                  ref_CDTS3score = CDTS_3ss_mean,
                  ref_cons5score = phastCons20way_5ss_mean,
                  ref_cons3score = phastCons20way_3ss_mean)
  
  df_all_introns %>% as_tibble()
  
  ######################################
  ## INTRONS - ADD THE BIOTYPE INFO
  ######################################
  
  print(paste0(Sys.time(), " - adding the transcript biotype ... "))
  
  df_protein <- readRDS(file = paste0(results.folder, "/all_split_reads_qc_level2_PC_biotype.rds")) %>%
    as_tibble()
  
  print(paste0(Sys.time(), " - ", intersect(df_all_introns$ref_junID, df_protein$junID) %>% length(), " total number of introns found!."))
  print(paste0(Sys.time(), " - ", df_all_introns %>% distinct(ref_junID) %>% nrow(), " total number of introns."))
  
  df_all_introns <- df_all_introns %>%
    left_join(y = df_protein,
              by = c("ref_junID" = "junID")) %>%
    dplyr::rename(protein_coding = percent_PC,
                  lncRNA = percent_lncRNA)
  
  print(paste0(Sys.time(), " - ", df_all_introns %>% distinct(ref_junID) %>% nrow(), " total number of introns after adding transcript biotype."))
  
  ######################################
  ## INTRONS - ADD THE CLINVAR DATA
  ######################################
  
  print(paste0(Sys.time(), " - adding the ClinVar data...")) 
  
  clinvar_tidy <- readRDS(file = paste0(dependencies_folder, "/clinvar_intronic_tidy.rda")) %>%
    GRanges() %>%
    diffloop::rmchr()
  clinvar_tidy %>% head()
  
  df_all_introns_tidy <- df_all_introns%>%
    mutate(clinvar_type = "-") %>%
    mutate(clinvar_start = 0) %>%
    mutate(clinvar_end = 0) %>%
    GRanges() 
  
  ## Find overlaps between clinvar mutations and mis-splicing ratios
  overlaps <- GenomicRanges::findOverlaps(query = clinvar_tidy,
                                          subject = df_all_introns_tidy,
                                          ignore.strand = F,
                                          type = "any")
  
  
  elementMetadata(df_all_introns_tidy)[subjectHits(overlaps), "clinvar_type"] <- ifelse(elementMetadata(clinvar_tidy)[queryHits(overlaps), "splice_site"] == "none",
                                                                                        "-", paste0(elementMetadata(clinvar_tidy)[queryHits(overlaps), "splice_site"], "-",
                                                                                                    elementMetadata(clinvar_tidy)[queryHits(overlaps), "distance"], "bp"))
  elementMetadata(df_all_introns_tidy)[subjectHits(overlaps), "clinvar_start"] <- clinvar_tidy[queryHits(overlaps), ] %>% start()
  elementMetadata(df_all_introns_tidy)[subjectHits(overlaps), "clinvar_end"] <- clinvar_tidy[queryHits(overlaps), ] %>% end()
  
  
  df_all_introns_tidy %>% head()
  
  
  ######################################
  ## INTRONS - ADD THE INTRON TYPE
  ######################################
  
  print(paste0(Sys.time(), " - adding the IAOD intron data...")) 
  
  ## Load intron type files
  u12_introns <- readRDS(file = paste0(dependencies_folder,
                                       "/minor_introns_tidy.rds"))
  
  ## Add two new columns to incorporate info about intron type
  df_all_introns_tidy <- df_all_introns_tidy %>%
    as_tibble() %>%
    mutate(u2_intron = T)
  
  
  ## MINOR INTRON
  print(paste0(Sys.time(), " - Getting junctions spliced out by the minor spliceosome."))
  overlaps <- GenomicRanges::findOverlaps(query = GenomicRanges::GRanges(seqnames = u12_introns$seqnames,
                                                                         ranges = IRanges(start = u12_introns$start, 
                                                                                          end = u12_introns$end),
                                                                         strand = u12_introns$strand),
                                          subject = GenomicRanges::GRanges(seqnames = df_all_introns_tidy$seqnames,
                                                                           ranges = IRanges(start = df_all_introns_tidy$start, 
                                                                                            end = df_all_introns_tidy$end),
                                                                           strand = df_all_introns_tidy$strand),
                                          ignore.strand = FALSE,
                                          type = "equal")
  
  queryHits(overlaps) %>% length()
  df_all_introns_tidy[subjectHits(overlaps),]$u2_intron <- F
  
  
  
  
  ## MAJOR INTRON
  #print(paste0(Sys.time(), " - Getting junctions spliced out by the major spliceosome."))
  
  #overlaps <- GenomicRanges::findOverlaps(query = GenomicRanges::GRanges(seqnames = u2_introns$seqnames %>% as.character(),
  #                                                                       ranges = IRanges(start = u2_introns$start,
  #                                                                                        end = u2_introns$end),
  #                                                                       strand = u2_introns$strand),
  #                                        subject = GenomicRanges::GRanges(seqnames = df_all_introns_tidy$seqnames,
  #                                                                         ranges = IRanges(start = df_all_introns_tidy$start,
  #                                                                                          end = df_all_introns_tidy$end),
  #                                                                         strand = df_all_introns_tidy$strand),
  #                                        ignore.strand = FALSE,
  #                                        type = "equal")
  
  # print(overlaps)
  #subjectHits(overlaps) %>% length()
  #df_all_introns_tidy[subjectHits(overlaps),]$u2_intron <- T
  

  df_all_introns_tidy$ref_junID %>% unique %>% length()
  df_all_introns_tidy %>% distinct(ref_junID, .keep_all = T) %>% dplyr::count(misspliced)
  

  ######################################
  ## INTRONS - POPULATE THE TABLE
  ######################################
  
  
  query <- paste0("CREATE TABLE IF NOT EXISTS 'intron'",
                  "(ref_junID NUMERIC PRIMARY KEY NOT NULL,
                  ref_coordinates TEXT NOT NULL, 
                  ref_length INTEGER NOT NULL, 
                  ref_ss5score DOUBLE NOT NULL, 
                  ref_ss3score DOUBLE NOT NULL, 
                  ref_cons5score DOUBLE NOT NULL,
                  ref_cons3score DOUBLE NOT NULL,
                  ref_CDTS5score DOUBLE NOT NULL,
                  ref_CDTS3score DOUBLE NOT NULL,
                  u2_intron BOOL,
                  clinvar TEXT NOT NULL, 
                  lncRNA INTEGER NOT NULL,
                  protein_coding INTEGER NOT NULL,
                  misspliced BOOL NOT NULL,
                  transcript_id INTEGER NOT NULL,
                  FOREIGN KEY (transcript_id) REFERENCES 'transcript'(id))")
  
  con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = database.path)
  dbListTables(con)
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  print("'Intron' table created!")
  
  
  ## POPULATE INTRON TABLE ----------------------------------------------------
  df_all_introns_tidy_final <- df_all_introns_tidy %>% 
    #as_tibble() %>%
    #unnest(gene_id) %>%
    dplyr::rename(ref_length = width) %>%
    mutate(ref_cons5score = ifelse(ref_cons5score %>% is.na(), 0, ref_cons5score),
           ref_cons3score = ifelse(ref_cons3score %>% is.na(), 0, ref_cons3score)) %>%
    dplyr::rename(ref_coordinates = ref_junID,
                  clinvar = clinvar_type) %>%
    distinct(ref_coordinates, .keep_all = T) %>%
    tibble::rowid_to_column("ref_junID") %>%
    dplyr::select(-c(seqnames, start, end, strand, 
                     clinvar_start, clinvar_end) )
  
  if (all(df_all_introns_tidy_final$misspliced == T)) {
    print("ERROR! all introns classified as misspliced!")
  }
  if (any(is.na(df_all_introns_tidy_final$lncRNA))) {
    print("ERROR! some introns do not have a biotype assigned")
  }
  if (any(df_all_introns_tidy_final$transcript_id %>% is.na())) {
    print("ERROR! some introns do not have a gene assigned")
  }
  
  DBI::dbAppendTable(conn = con,
                     name = "intron", 
                     value = df_all_introns_tidy_final)
  
  
  
  print(paste0("'Intron' table populated! ", 
               df_all_introns_tidy_final %>% distinct(ref_junID) %>% nrow(), " introns stored!" ))
  
  ## CREATE INDEX TO SPEED UP QUERIES ------------------------------------------
  query <- paste0("CREATE UNIQUE INDEX 'index_intron' ON 'intron'(ref_junID)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  #setdiff(df_never_misspliced$ref_coordinates, df_all_introns_tidy_final$ref_coordinates) %>% length()
  
  
  ########################################
  ## B) NOVEL JUNCTION - ADD MAXENTSCAN INFO 
  ########################################
  
  print(paste0(Sys.time(), " - adding MaxEntScan scores to the NOVEL JUNCTONS..."))
  
  df_all_novel_raw_tidy <- df_all_distances_pairings_raw %>%
    mutate( start = novel_start %>% as.integer(),
            end = novel_end %>% as.integer()) %>%
    dplyr::select(seqnames = novel_seq,
                  start, end, 
                  strand = novel_strand,
                  novel_junID, ref_junID, 
                  novel_type = type, distance) %>%
    distinct(novel_junID, .keep_all = T)
  
  wd <- getwd()
  
  ## Add MaxEntScan score to the split reads
  all_split_reads_tidy <- generate_max_ent_score(junc_tidy = df_all_novel_raw_tidy %>% dplyr::rename(junID = novel_junID),
                                                 max_ent_tool_path = paste0(dependencies_folder,"/fordownload/"),
                                                 homo_sapiens_fasta_path = paste0(dependencies_folder,
                                                                                  "/Homo_sapiens.GRCh38.dna.primary_assembly.fa"))
  
  all_split_reads_tidy %>% as_tibble()
  
  all_split_reads_tidy <- all_split_reads_tidy %>% 
    dplyr::select(-donorSeqStart,
                  -donorSeqStop,
                  -AcceptorSeqStart,
                  -AcceptorSeqStop,
                  -donor_sequence,
                  -acceptor_sequence)
  
  setwd(wd)
  
  
  df_all_novels_tidy <- all_split_reads_tidy %>%
    mutate(novel_junID = paste0("chr", seqnames, ":", start, "-", end, ":", strand)) 
  
  
  
  df_all_novels_tidy <- df_all_novels_tidy %>%
    dplyr::select(-one_of("junID","novel_ss5score","novel_ss3score")) %>% 
    dplyr::rename(novel_ss5score = ss5score, 
                  novel_ss3score = ss3score) %>%
    dplyr::relocate(ref_junID, novel_junID) %>%
    dplyr::relocate(c(novel_ss5score, novel_ss3score), .before = novel_type ) 
  
  
  df_all_novels_tidy <- df_all_novels_tidy %>% as_tibble()
  df_all_novels_tidy
  
  ####################################
  ## NOVEL - ADD INTRON REFERENCE
  ####################################
  
  print(paste0(Sys.time(), " - adding the INTRON foreign key to the NOVEL JUNCTONS..."))
  
  df_all_novels_tidy_final <- df_all_novels_tidy %>% 
    inner_join(y = df_all_introns_tidy_final %>%
                 dplyr::select(ref_junID, ref_coordinates),
               by = c("ref_junID" = "ref_coordinates" )) %>%
    dplyr::select(-ref_junID) %>%
    dplyr::rename(ref_junID = ref_junID.y) %>% 
    as_tibble() %>%
    dplyr::relocate(ref_junID)
  
  
  
  # if (any( df_all_novels_tidy_final$novel_junID %in% novel_jxn_minor_introns$novel_junID)) {
  #   print("ERROR! Some novel junctions attached to introns spliced out by the minor spliceosome
  #         have been attached to introns spliced out by the major spliceosome!")
  #   break;
  # }
  
  ####################################
  ## CREATE NOVEL JUNCTION TABLE
  ####################################
  
  query <- paste0("CREATE TABLE IF NOT EXISTS 'novel'",
                  "(novel_junID NUMERIC NOT NULL,
                  ref_junID NUMERIC NOT NULL,
                  novel_coordinates TEXT NOT NULL, 
                  novel_ss5score DOUBLE NOT NULL, 
                  novel_ss3score DOUBLE NOT NULL,
                  novel_length INTEGER NOT NULL, 
                  novel_type TEXT NOT NULL, 
                  distance INTEGER NOT NULL,
                  PRIMARY KEY (ref_junID, novel_junID),
                  FOREIGN KEY (ref_junID) REFERENCES 'intron'(ref_junID))")
  
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  print("'Novel' table created!")
  
  df_all_novels_tidy_final <- df_all_novels_tidy_final %>% 
    as_tibble() %>%
    dplyr::rename(novel_coordinates = novel_junID ) %>%
    distinct(novel_coordinates, .keep_all = T) %>%
    GRanges() %>%
    as_tibble() %>%
    tibble::rowid_to_column("novel_junID") %>%
    dplyr::select(-seqnames, -start, -end, -strand) %>%
    dplyr::rename(novel_length = width)
  
  df_all_novels_tidy_final
  if (any(duplicated(df_all_novels_tidy_final$novel_coordinates))) {
    print("ERROR! some novel junctions are duplicated")
    break;
  }
  DBI::dbAppendTable(conn = con,
                     name = "novel", 
                     value = df_all_novels_tidy_final)
  
  

  print(paste0("'Novel' table populated! ", 
               df_all_novels_tidy_final %>% distinct(novel_junID) %>% nrow(), 
               " novel junctions stored!" ))
  
}


create_gene_table <- function(database.path,
                              gene_ids) {
  
  
  print(paste0(Sys.time(), " --> creating the gene table..."))
  con <- dbConnect(RSQLite::SQLite(), database.path)
  
  
  hg38_transcripts_gene <- hg38 %>%
    as_tibble() %>%
    #dplyr::select(gene_id, gene_name) %>%
    mutate(gene_id = str_sub(gene_id, start = 1, end = 15)) %>%
    #distinct(gene_id, .keep_all = T) %>%
    dplyr::filter(gene_id %in% gene_ids$gene_id) %>%
    dplyr::count(gene_id, type) %>%
    dplyr::filter(type == "transcript") %>%
    unnest(gene_id) %>%
    dplyr::select(-type) %>%
    dplyr::rename(n_transcripts = n)
  
  hg38_genes <- hg38 %>%
    as_tibble() %>%
    mutate(gene_id = str_sub(gene_id, start = 1,end = 15)) %>%
    dplyr::filter(gene_id %in% gene_ids$gene_id) %>%
    dplyr::filter(type == "gene") %>%
    dplyr::select(gene_id, gene_name, gene_width = width)
  
  hg38_tidy <- hg38_transcripts_gene %>%
    inner_join(hg38_genes, by = "gene_id" ) %>% 
    as_tibble()
  
  
  ## CREATE GENE_NAME TABLE
  sql_statement <- paste0("CREATE TABLE IF NOT EXISTS 'gene'", 
                          "(id INTEGER PRIMARY KEY NOT NULL,
                          gene_id TEXT NOT NULL,
                          gene_name TEXT,
                          n_transcripts INTEGER NOT NULL,
                          gene_width INTEGER NOT NULL)")
  res <- DBI::dbSendQuery(conn = con, statement = sql_statement)
  DBI::dbClearResult(res)
  
  
  ## POPULATE GENE TABLE
  hg38_tidy_final <- hg38_tidy %>% 
    as_tibble() %>%
    tibble::rowid_to_column("id")
  
  DBI::dbAppendTable(conn = con,
                     name = "gene", 
                     value = hg38_tidy_final)
  
  
  print(paste0(Sys.time(), " --> gene table created!"))
  print(paste0(Sys.time(), " --> ", hg38_tidy_final %>% nrow(), " genes stored!"))
}


create_transcript_table  <- function(database.path,
                                     gene_ids,
                                     tx_ids) {
  
  print(paste0(Sys.time(), " --> creating the 'transcript' table..."))
  con <- dbConnect(RSQLite::SQLite(), database.path)
  
  ## GET TRANSCRIPT ID
  hg38_transcripts_gene <- hg38 %>%
    as_tibble() %>%
    dplyr::filter(transcript_id %in% tx_ids$tx_id_junction,
                  type == "transcript") %>%
    mutate(gene_id = str_sub(gene_id, start = 1, end = 15)) %>%
    #distinct(gene_id, .keep_all = T) %>%
    dplyr::filter(gene_id %in% gene_ids$gene_id) %>%
    dplyr::select(transcript_id, TSL = transcript_support_level, gene_id) %>%
    mutate(TSL = str_sub(TSL, start = 1, end = 2) %>% as.integer())
  hg38_transcripts_gene[is.na(hg38_transcripts_gene[,"TSL"]),"TSL"] <- 10
  
  ## ADD MANE INFO
  hg38_transcripts_gene_mane <- hg38_transcripts_gene %>%
    left_join(y = hg_mane_transcripts,
              by = "transcript_id" )
  hg38_transcripts_gene_mane[is.na(hg38_transcripts_gene_mane[,"MANE"]),"MANE"] <- F
  
  
  
  ## ADD THE GENE FOREING KEY
  print(paste0(Sys.time(), " --> adding GENE foreing key to the transcripts..."))
  
  query <- paste0("SELECT id, gene_id FROM 'gene'")
  df_genes <- dbGetQuery(con, query) %>% as_tibble()
  
  ## Add the GENE ID for the foreign key
  hg38_transcripts_gene_mane <- hg38_transcripts_gene_mane %>%
    inner_join(df_genes %>% dplyr::select(id, gene_id),
               by = "gene_id") %>%
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
                  gene_id INTEGER NOT NULL,
                  PRIMARY KEY (id),
                  FOREIGN KEY (gene_id) REFERENCES 'gene'(id))")
  
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  print(paste0(Sys.time(), " --> 'Transcript' table created!"))
  
  
  if ( any(duplicated(hg38_transcripts_final$transcript_id)) ) {
    print("ERROR! some novel junctions are duplicated")
  }
  DBI::dbAppendTable(conn = con,
                     name = "transcript", 
                     value = hg38_transcripts_final)
  
  
  print( paste0(Sys.time(), " --> 'Transcript' table populated!") )
  print(paste0(Sys.time(), " --> ", hg38_transcripts_final %>% nrow(), " transcripts stored!"))
}


## Create the cluster tables ----------------------------------------------------

create_cluster_tables <- function(database.path,
                                  project.name,
                                  recount3.project.IDs = NULL,
                                  gtf.version,
                                  database.folder,
                                  results.folder) {
  
  all_split_reads_details_105 <- readRDS(file = paste0(database.folder, "/all_split_reads_qc_level2.rds"))
  
  con <- dbConnect(RSQLite::SQLite(), database.path)
  tables <- DBI::dbListTables(conn = con)
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  
  print( paste0(Sys.time(), " --> SQL connection stablished!") )
  
  ## GET FROM MASTER TABLE
  query = paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  
  ## GET FROM INTRON TABLE
  query = paste0("SELECT * FROM 'intron'")
  df_intron <- dbGetQuery(con, query) %>% as_tibble()
  df_intron %>% nrow()
  df_intron %>% 
    dplyr::count(misspliced)
  
  ## GET FROM NOVEL JUNCTION TABLE
  query = paste0("SELECT * FROM 'novel'")
  df_novel <- dbGetQuery(con, query) %>% as_tibble()
  df_novel %>% nrow()
  df_novel %>% 
    dplyr::count(novel_type) %>%
    print()
  
  
  ## GET FROM GENE TABLE
  query = paste0("SELECT * FROM 'gene'")
  df_gene <- dbGetQuery(con, query) %>% as_tibble()
  df_gene %>% nrow()
  
  ## GET FROM TRANSCRIPT TABLE
  query = paste0("SELECT * FROM 'transcript'")
  df_transcript <- dbGetQuery(con, query) %>% as_tibble()
  df_transcript %>% nrow()
  
  if ( is.null(recount3.project.IDs) ){
    recount3.project.IDs <- (df_metadata$SRA_project %>% unique())
  }
  
  for (project_id in recount3.project.IDs) {
    
    # project_id <- recount3.project.IDs[1]
    
    print(paste0(Sys.time(), " --> Working with '", project_id, "' DataBase..."))
    results_folder_local <- paste0(results.folder, "/", project_id, "/")
    
    clusters <- df_metadata %>%
      dplyr::filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    for (cluster_id in clusters) { 
      
      # cluster_id <- clusters[2]
      print(paste0(Sys.time(), " --> ", cluster_id))
      
      ###############################
      ## PREPARE DATA
      ###############################
      #if ( !any(tables == (paste0(cluster_id, "_", project_id, "_misspliced"))) ) {
        
      if ( file.exists(paste0(results_folder_local, "/junction_pairing/", 
                              cluster_id, "/", cluster_id, "_raw_distances_tidy.rds")) ) {
        
        
        ## BASE DATA -----------------------------------------------------------
        
        print(paste0(Sys.time(), " --> load base data ... "))
        
        ## Load split read counts
        split_read_counts <- readRDS(file = paste0(results_folder_local, "/base_data/", 
                                                   project_id, "_", cluster_id, "_split_read_counts.rds")) 
        
        if ( is.null(names(split_read_counts)) ) {
          split_read_counts <- split_read_counts %>%
            as_tibble(rownames = "junID")
        }
        
        ## Load samples
        samples <- readRDS(file = paste0(results_folder_local, "/base_data/", 
                                         project_id, "_", cluster_id, "_samples_used.rds"))
        
        
        
        ## LOAD INTRONS AND NOVEL JUNCTIONS ------------------------------------
        df_cluster_distances <- readRDS(file = paste0(results_folder_local, "/junction_pairing/", cluster_id, "/", 
                                                      cluster_id, "_raw_distances_tidy.rds")) %>% as_tibble()
        
        
        
        ## INTRONS ---------------------------------------------------
        df_introns_gr <- df_cluster_distances %>%
          distinct(ref_junID, .keep_all = T) %>%
          dplyr::select(ref_junID, seqnames = ref_seq, start = ref_start,
                        end = ref_end, strand = ref_strand)
        
        
        
        ## Add reads detected for the introns in the current tissue
        split_read_counts_intron <- get_mean_coverage(split_read_counts = split_read_counts,
                                                      samples = samples,
                                                      junIDs = df_introns_gr$ref_junID) %>%
          dplyr::rename(ref_n_individuals = n_individuals,
                        ref_sum_counts = sum_counts)
        
        split_read_counts_intron %>% head()
        if ( split_read_counts_intron$ref_sum_counts %>% unique() %>% min() == 1 ) {
          print("ERROR! Minimum number of supporting reads is 1!")
          break;
        }
        
        df_introns_merged <- df_introns_gr %>%
          left_join(y = split_read_counts_intron,
                    by = c("ref_junID" = "junID"))
        
        
        ## NOVEL JUNCTIONS -----------------------------------------------------
        df_novel_gr <- df_cluster_distances %>%
          distinct(novel_junID, .keep_all = T) %>%
          dplyr::select(novel_junID, ref_junID, seqnames = novel_seq, 
                        start = novel_start, end = novel_end, strand = novel_strand)
        
        
        split_read_counts_novel <- get_mean_coverage(split_read_counts = split_read_counts,
                                                     samples = samples,
                                                     junID = df_novel_gr$novel_junID) %>%
          dplyr::rename(novel_n_individuals = n_individuals,
                        novel_sum_counts = sum_counts)
        
        
        df_novel_merged <- df_novel_gr %>%
          left_join(y = split_read_counts_novel,
                    by = c("novel_junID" = "junID"))
        
        df_novel_merged %>% as_tibble()
        
        
        ## QC data - NO '*' within the ID ---------------------------------------------
        df_novel_merged <- df_novel_merged %>% 
          rowwise() %>%
          mutate(ref_junID = ifelse(str_detect(string = ref_junID, pattern = "\\*"), 
                                    str_replace(string = ref_junID, pattern = "\\*", strand ),
                                    ref_junID)) 
        
        df_novel_merged <- df_novel_merged %>% 
          rowwise() %>%
          mutate(novel_junID = ifelse(str_detect(string = novel_junID, pattern = "\\*"),
                                      str_replace(string = novel_junID, pattern = "\\*", strand ),
                                      novel_junID)) 
        
        df_introns_merged <- df_introns_merged %>% 
          rowwise() %>%
          mutate(strand = strand %>% as.character()) %>%
          mutate(ref_junID = ifelse(str_detect(string = ref_junID, pattern = "\\*"),
                                    str_replace(string = ref_junID, pattern = "\\*", strand ),
                                    ref_junID))
        
        
        ## Merge tidy data -----------------------------------------------------
        
        print(paste0(Sys.time(), " --> merge introns and novel junctions ... "))
        
        df_intron_novel_merged <- df_novel_merged %>%
          dplyr::select(-c(seqnames, start, end, strand)) %>%
          inner_join(y = df_introns_merged %>% 
                       dplyr::select(-c(seqnames, start, end, strand)),
                     by = "ref_junID")
        
        
        
        
        ####################################
        ## QC
        ####################################
        
        if (any(str_detect(df_intron_novel_merged$ref_junID,pattern = "//*"))) {
          print("ERROR: some IDs contain '*'")
        }
        if (any(str_detect(df_intron_novel_merged$novel_junID,pattern = "//*"))) {
          print("ERROR: some IDs contain '*'")
        }
        if (any(df_intron_novel_merged$ref_junID %>% is.na())) {
          print("ERROR!!")
        }
        
        
        
        ## JOIN data with MASTER NOVEL table
        ## Only novel junctions previously paired and QC'ed are considered
        
        print(paste0(Sys.time(), " --> merge data with MASTER NOVEL table ... "))
        
        df_all_misspliced <- df_intron_novel_merged %>% 
          inner_join(y = df_novel %>% 
                       dplyr::select(novel_junID, novel_coordinates, novel_type),
                     by = c("novel_junID" = "novel_coordinates")) %>%
          dplyr::rename(novel_coordinates = novel_junID) %>%
          dplyr::rename(novel_junID = novel_junID.y)
        
        
        if(setdiff(df_all_misspliced$ref_junID, df_intron$ref_coordinates) %>% length() > 0) {
          print("ERROR! introns detected in this tissue are not stored in the master intron table.")
          break;
        }
        
        ## JOIN data with MASTER INTRON table
        df_all_misspliced <- df_all_misspliced %>% 
          left_join(y = df_intron %>%
                      dplyr::filter(misspliced == T) %>%
                      dplyr::select(ref_junID, ref_coordinates, transcript_id),
                    by = c("ref_junID" = "ref_coordinates")) %>%
          dplyr::select(-ref_junID) %>%
          dplyr::rename(ref_junID = ref_junID.y)
        
        
        if (which(str_detect(df_all_misspliced$novel_coordinates,pattern = "//*")) %>% length() > 0) {
          print("ERROR: some IDs contain '*'")
        }
        
        
        ## PREPARE DATA PRIOR POPULATING THE TABLE
        df_all_misspliced <- df_all_misspliced %>%
          relocate(ref_junID, novel_junID)
        
        
        
        #######################################
        ## CHECK INTEGRITY WITH PARENT TABLES
        #######################################
        
        print(paste0(Sys.time(), " --> checking integrity with parent tables ... "))
        
        master_novel <- df_novel %>%
          dplyr::filter(novel_junID %in% 
                          (df_all_misspliced %>%
                             pull(novel_junID))) %>% 
          dplyr::select(novel_coordinates) 
        
        if ( !(identical(df_all_misspliced$novel_coordinates %>% sort(), 
                         master_novel$novel_coordinates %>% sort())) ) {
          print("ERROR! Tables not identical")
          if ( all(intersect(setdiff(df_all_misspliced$novel_coordinates, 
                                     master_novel$novel_coordinates),
                             ambiguous_novel_junc$novel_junID) == setdiff(df_all_misspliced$novel_coordinates, 
                                                                          master_novel$novel_coordinates)) == T ) {
            df_all_misspliced <- df_all_misspliced %>%
              as.data.table() %>%
              dplyr::filter(!(novel_coordinates %in% ambiguous_novel_junc$novel_junID))
          }
        }
        
        if ( identical(df_all_misspliced$novel_coordinates %>% sort(), 
                       master_novel$novel_coordinates %>% sort()) ) {
          
          df_all_misspliced <- df_all_misspliced %>%
            dplyr::select(-novel_coordinates)
          
          
          ## CHECK PARENT INTEGRITY
          
          df <- df_novel %>%
            dplyr::select(novel_junID, ref_junID) %>%
            arrange(novel_junID) %>% 
            inner_join(df_all_misspliced %>%
                         dplyr::select(novel_junID, ref_junID) %>%
                         arrange(novel_junID),
                       by = "novel_junID")
          
          
          
          diff <- df %>% 
            dplyr::filter(ref_junID.x != ref_junID.y)
          
          
          if (diff %>% nrow() > 0) {
            df_all_misspliced <- df_all_misspliced %>%
              dplyr::filter(!(ref_junID %in% diff$ref_junID.y)) 
            print(paste0("ERROR!: ", diff, " --> mismatch junctions."))
            break;
          } else {
            print("good")
          }
          
          
          #####################################
          ## CALCULATE MSR MEASURES
          #####################################
          
          print(paste0(Sys.time(), " --> calculating MSR measures ... "))
          
          # df_all_misspliced %>% dplyr::filter(ref_junID == "17") %>% as.data.frame()
          
          db_introns <- df_all_misspliced %>%
            group_by(ref_junID, novel_type) %>%
            mutate(MSR = sum(novel_sum_counts)/(sum(novel_sum_counts) + ref_sum_counts))
          
          # db_introns %>% dplyr::filter(ref_junID == "17") %>% as.data.frame()
          
          db_introns <- db_introns %>% 
            spread(key = novel_type, value = MSR, fill = 0) %>%
            group_by(ref_junID) %>% 
            dplyr::mutate(MSR_Donor = max(novel_donor, na.rm = T)) %>%
            dplyr::mutate(MSR_Acceptor = max(novel_acceptor, na.rm = T))  %>%
            ungroup()
          
          
          #####################################
          ## GET THE GENE TPM
          #####################################
          
          if ( file.exists( paste0(results.folder, "/", project_id, "/tpm/", project_id, "_", cluster_id, "_tpm.rds")) ) {
            
            tpm <- readRDS(file = paste0(results.folder,  "/", project_id, "/tpm/", project_id, "_", cluster_id, "_tpm.rds")) %>% 
              dplyr::select(gene_id = gene, all_of(samples))
            
            
            tpm <- tpm  %>%
              dplyr::mutate(tpm_median = matrixStats::rowMedians(x = as.matrix(.[2:(ncol(tpm))]))) %>%
              dplyr::select(gene_id, tpm_median) 
            
            
            ## In case there are any duplicates, take the genes with the maximum tpms
            tpm <- tpm %>% 
              distinct(gene_id, .keep_all = T) %>%
              group_by(gene_id) %>% 
              summarize_all(max)
            
            tpm_tidy <- tpm %>%
              inner_join(y = df_gene %>%  as_tibble(),
                         by = c("gene_id" = "gene_id")) %>%
              inner_join(y = df_transcript %>% 
                           as_tibble(),
                         by = c("id" = "gene_id"),
                         multiple = "all") %>%
              dplyr::select(transcript_id = id.y, 
                            tpm_median)
            
            db_introns <- db_introns %>%
              left_join(y = tpm_tidy,
                        by = c("transcript_id" = "transcript_id")) %>% 
              dplyr::rename(gene_tpm = tpm_median)
          }
          
          
          #####################################
          ## ADD THE REFERENCE TYPE
          #####################################
          
          
          db_introns[, "ref_type"] <- ""
          db_introns %>% head()
          
          ## TYPE 'BOTH'
          indx <- which(db_introns$MSR_Donor > 0 & db_introns$MSR_Acceptor > 0)
          db_introns[indx, "ref_type"] <- "both"
          print(paste0("Introns type 'both': ", indx %>% length()))
          
          
          
          ## TYPE 'DONOR'
          indx <- which(db_introns$MSR_Donor > 0 & db_introns$MSR_Acceptor == 0)
          db_introns[indx, "ref_type"] <- "donor"
          print(paste0("Junctions type 'donor': ", indx %>% length()))
          
          
          ## TYPE 'ACCEPTOR'
          indx <- which(db_introns$MSR_Donor == 0 & db_introns$MSR_Acceptor > 0)
          db_introns[indx, "ref_type"] <- "acceptor"
          print(paste0("Junctions type 'acceptor': ", indx %>% length()))
          
          
          
          
          
          db_introns %>% nrow()
          db_introns %>% distinct(ref_junID) %>% nrow()
          
          
          ###############################
          ## CREATE INTRON TABLE
          ###############################
          
          print(paste0(Sys.time(), " --> creating 'intron' table ... "))
          
          # dbRemoveTable(conn = con, paste0(cluster_id, "_", project_id))
          query <- paste0("CREATE TABLE IF NOT EXISTS '", paste0(cluster_id, "_", project_id, "_misspliced"), "'", 
                          "(ref_junID INTEGER NOT NULL,
                          novel_junID INTEGER NOT NULL,
                          ref_n_individuals INTEGER NOT NULL,
                          ref_sum_counts INTEGER NOT NULL,
                          ref_type TEXT NOT NULL, 
                          novel_n_individuals INTEGER NOT NULL, 
                          novel_sum_counts INTEGER NOT NULL, 
                          MSR_D DOUBLE NOT NULL, 
                          MSR_A DOUBLE NOT NULL, 
                          gene_tpm DOUBLE, 
                          transcript_id INTEGER NOT NULL, 
                          FOREIGN KEY (ref_junID, novel_junID) REFERENCES novel (ref_junID, novel_junID),
                          FOREIGN KEY (transcript_id) REFERENCES 'transcript'(id))")
          
          res <- DBI::dbSendQuery(conn = con, statement = query)
          DBI::dbClearResult(res)
          
          
          #####################################
          ## POPULATE THE TABLE
          #####################################
          
          
          db_introns_final <- db_introns %>%
            dplyr::select(-novel_acceptor,-novel_donor)%>%
            dplyr::rename(MSR_D = MSR_Donor, MSR_A = MSR_Acceptor)
          
          DBI::dbAppendTable(conn = con,
                             name = paste0(cluster_id, "_", project_id, "_misspliced"), 
                             value = db_introns_final)
          
          ## CREATE INDEX
          query <- paste0("CREATE UNIQUE INDEX 'index_", paste0(cluster_id, "_", project_id, "_misspliced"), "' ON '",
                          paste0(cluster_id, "_", project_id, "_misspliced"), "'(ref_junID,novel_junID)");
          res <- DBI::dbSendQuery(conn = con, statement = query)
          DBI::dbClearResult(res)
          
          
          print(paste0(Sys.time(), ". '", paste0(cluster_id, "_", project_id, "_misspliced"), 
                       "' table populated!"))
          
          ####################################
          ## NEVER MISSPLICED
          ####################################
          
          print(paste0(Sys.time(), " --> getting never mis-spliced introns ... "))
          
          ## TYPE 'NONE'
          introns_never <- readRDS(file = paste0(results_folder_local, "/junction_pairing/", cluster_id, "/not-misspliced/", 
                                                 cluster_id, "_all_notmisspliced.rds"))
          ## The introns not misspliced in this tissue, should have not been detected as spliced.
          ## Thus, this should be zero
          if (intersect(introns_never, df_intron %>%
                        dplyr::filter(ref_junID %in% db_introns_final$ref_junID) %>%
                        pull(ref_coordinates) %>% 
                        unique()) %>% length() > 0 ){
            print("ERROR! The introns not misspliced in this tissue, should have not been detected as spliced.")
            break;
          }
          df_introns_never <- data.frame(ref_junID = introns_never)
          
          
          
          if (any(str_detect(string = df_introns_never$ref_junID, pattern = "\\*"))) {
            print("ERROR! * in the IDs")
            break;
          }
          
          split_read_counts <- split_read_counts %>%
            left_join(y = all_split_reads_details_105 %>% 
                        dplyr::select(junID, seqnames, start, end, strand),
                      by = "junID") %>%
            dplyr::select(-junID) %>%
            mutate(junID = paste0("chr", seqnames, ":", start, "-", end, ":", strand)) %>%
            dplyr::relocate(junID)
          
          split_read_counts_intron_never <- get_mean_coverage(split_read_counts = split_read_counts,
                                                              samples = samples,
                                                              junID = df_introns_never$ref_junID) %>%
            dplyr::rename(ref_n_individuals = n_individuals,
                          ref_sum_counts = sum_counts) 
          
          if (any(str_detect(split_read_counts_intron_never$junID,pattern = "\\*"))) {
            print("ERROR! some never mis-spliced junctions without the number of individuals")
            break;
          }
          
          
          df_never_merged <- df_introns_never %>%
            inner_join(y = split_read_counts_intron_never,
                       by = c("ref_junID" = "junID")) %>%
            mutate(MSR_D = 0, MSR_A = 0) %>%
            mutate(ref_type = "never")
          df_never_merged <- df_never_merged %>% as_tibble()
          
          
          
          if (any(df_never_merged$ref_n_individuals %>% is.na)) {
            print("ERROR! some never mis-spliced junctions without the number of individuals")
            break;
          }
          
          
          ## QC
          if (any(str_detect(df_never_merged$ref_junID, pattern = "\\*"))) {
            print("ERROR! * in the IDs")
            break;
          }
          
          
          ## TYPE 'NONE'
          print(paste0("Junctions type 'never': ", df_never_merged %>% distinct(ref_junID) %>% nrow()))
          
          ##################################################
          ## ADD THE INTRON REFERENCE KEY
          ##################################################
          
          print(paste0(Sys.time(), " --> adding the intron reference key to the never mis-spliced introns ... "))
          
          df_never_merged <- df_never_merged %>%
            inner_join(df_intron %>% 
                         dplyr::select(ref_junID, ref_coordinates, transcript_id),
                       by = c("ref_junID" = "ref_coordinates")) %>%
            dplyr::filter(!is.na(ref_junID)) %>%
            dplyr::select(-ref_junID) %>% 
            dplyr::rename(ref_junID = ref_junID.y) %>%
            relocate(ref_junID)
          
          
          if (df_never_merged %>% dplyr::filter(is.na(ref_junID)) %>% nrow > 0) {
            print("ERROR! IDs are NA")
          }
          if (any(df_never_merged$ref_type != "never")) {
            print("Error! Some never mis-spliced junctions have been stored as mis-spliced.")
          }
          if ((intersect(df_never_merged$ref_junID, db_introns$ref_junID) %>% length()) > 0) {
            print("Error! Some never mis-spliced junctions have been stored as mis-spliced.")
          }
          if (any(duplicated(df_never_merged$ref_junID))) {
            print("Error! Some never mis-spliced junctions are duplicated.")
          }
          
          ## Remove introns and novel junctions with not TPM data
          # df_never <- df_never %>%
          #   dplyr::filter(!(gene_tpm %>% is.na()))    
          
          #####################################
          ## GET THE GENE TPM
          #####################################
          
          
          if (exists("tpm_tidy")) {
            df_never_merged <- df_never_merged %>%
              left_join(y = tpm_tidy,
                        by = "transcript_id") %>% 
              dplyr::rename(gene_tpm = tpm_median) 
          }
          
          ## QC --------------------------------------------------------
          if (intersect(df_never_merged %>%
                        dplyr::filter(ref_type == "never") %>%
                        distinct(ref_junID) %>%
                        pull(), db_introns %>%
                        dplyr::filter(ref_type == "donor") %>%
                        distinct(ref_junID) %>%
                        pull()) %>% length() > 0) {
            print("Error!")
          }
          
          
          if (intersect(df_never_merged %>%
                        dplyr::filter(ref_type == "never") %>%
                        distinct(ref_junID) %>%
                        pull(), db_introns %>%
                        dplyr::filter(ref_type == "aceptor") %>%
                        distinct(ref_junID) %>%
                        pull()) %>% length() > 0) {
            print("Error!")
          }
          
          if (intersect(df_never_merged %>%
                        dplyr::filter(ref_type == "never") %>%
                        distinct(ref_junID) %>%
                        pull(), db_introns %>%
                        dplyr::filter(ref_type == "both") %>%
                        distinct(ref_junID) %>%
                        pull()) %>% length() > 0) {
            print("Error!")
          }
          
          if (intersect(db_introns %>%
                        dplyr::filter(ref_type == "acceptor") %>%
                        distinct(ref_junID) %>%
                        pull(), db_introns %>%
                        dplyr::filter(ref_type == "donor") %>%
                        distinct(ref_junID) %>%
                        pull()) %>% length() > 0) {
            print("Error!")
          }
          
          
          ###############################
          ## CREATE NEVER MIS-SPLICED TABLE
          ###############################
          
          query <- paste0("CREATE TABLE IF NOT EXISTS '", paste0(cluster_id, "_", project_id, "_nevermisspliced"), "'", 
                          "(ref_junID INTEGER NOT NULL,
                          ref_n_individuals INTEGER NOT NULL, 
                          ref_sum_counts INTEGER NOT NULL,
                          MSR_D DOUBLE NOT NULL, 
                          MSR_A DOUBLE NOT NULL, 
                          ref_type TEXT NOT NULL, 
                          gene_tpm DOUBLE,
                          transcript_id INTEGER NOT NULL,
                          FOREIGN KEY (ref_junID) REFERENCES intron (ref_junID),
                          FOREIGN KEY (transcript_id) REFERENCES 'transcript'(id))")
          
          res <- DBI::dbSendQuery(conn = con, statement = query)
          DBI::dbClearResult(res)
          
          #gene_tpm DOUBLE,
          ###############################
          ## POPULATE TABLE
          ###############################
          
          DBI::dbAppendTable(conn = con,
                             name = paste0(cluster_id, "_", project_id, "_nevermisspliced"), 
                             value = df_never_merged)
          
          
          ## CREATE INDEX
          query <- paste0("CREATE UNIQUE INDEX 'index_", 
                          paste0(cluster_id, "_", project_id, "_nevermisspliced"), "' ON '",
                          paste0(cluster_id, "_", project_id, "_nevermisspliced"),"'(ref_junID)");
          res <- DBI::dbSendQuery(conn = con, statement = query)
          DBI::dbClearResult(res)
          
          print(paste0(Sys.time(), ". '", 
                       paste0(cluster_id, "_", project_id, "_nevermisspliced"), 
                       "' table populated!"))
          
        } else {
          print("Error: novel junctions are distinct!")
          break;
        }
        
      }
      #}
      
    }
    gc()
  }
  
  
  dbListTables(con)
  DBI::dbDisconnect(conn = con)
  
}

