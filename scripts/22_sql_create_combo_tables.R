
#' Title
#' Creates the child tables for each of the recount3 ID projects indicated
#' @param database.path Local path to the .sqlite database
#' @param recount3.project.IDs Vector with the ID of the recount3 projects to work with
#' @param database.folder Local path to the folder that contains the database
#' @param results.folder Local path to the folder that contains the result files
#'
#' @return
#' @export
#'
#' @examples
SqlCreateComboTables <- function(database.path,
                                 database.folder,
                                 results.folder,
                                 dependencies.folder,
                                 recount3.project.IDs = NULL) {
  
  
  con <- dbConnect(RSQLite::SQLite(), database.path)
  tables <- DBI::dbListTables(conn = con)
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  
  logger::log_info( " --> SQL connection stablished!") 
  logger::log_info("Querying master tables ...")
  
  ## GET FROM MASTER TABLE
  query = paste0("SELECT * FROM 'metadata'")
  df_metadata <- dbGetQuery(con, query) 
  
  ## GET FROM GENE TABLE
  query = paste0("SELECT * FROM 'gene'")
  master_gene <- dbGetQuery(con, query) %>% as_tibble()
  master_gene %>% nrow()
  
  ## GET FROM TRANSCRIPT TABLE
  query = paste0("SELECT * FROM 'transcript'")
  master_transcript <- dbGetQuery(con, query) %>% as_tibble()
  master_transcript %>% nrow()
  
  DBI::dbDisconnect(conn = con) 
  
  if ( is.null(recount3.project.IDs) ){
    recount3.project.IDs <- (df_metadata$SRA_project %>% unique())
  }

  
  for ( project_id in recount3.project.IDs ) { 
    
    # project_id <- recount3.project.IDs[1]
    # project_id <- recount3.project.IDs[2]
    
    logger::log_info(" --> Working with '", project_id, "' ...")
    results_folder_local <- paste0(results.folder, "/", project_id, "/")
    
    clusters <- df_metadata %>%
      dplyr::filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    for (cluster_id in clusters) { 
      
      # cluster_id <- clusters[1]
      # cluster_id <- clusters[2]
      
      logger::log_info(project_id, " --> ", cluster_id)
      
      ###############################
      ## PREPARE DATA
      ###############################
      
      
      
      if ( file.exists(paste0(results_folder_local, "/base_data/", project_id, "_", cluster_id, "_all_split_reads_combos.rds")) && 
           file.exists(paste0(results_folder_local, "/base_data/", project_id, "_", cluster_id, "_split_read_counts_combos.rds")) ) {
        
      
        ## LOAD BASE DATA ONLY FOR THE CURRENT CLUSTER ID -------------------------------------------------
        
        logger::log_info(" --> ", cluster_id, " load base data ... ")
        
        
        ## Load all split reads
        all_split_reads <- readRDS(file = paste0(results_folder_local, "/base_data/", 
                                                 project_id, "_", cluster_id, "_all_split_reads_combos.rds")) %>%
          dplyr::select(-c("in_ref","n_projects","annotated")) %>%
          unnest(gene_id)
        
        ## Load split read counts
        split_read_counts <- readRDS(file = paste0(results_folder_local, "/base_data/", 
                                                   project_id, "_", cluster_id, "_split_read_counts_combos.rds")) 
        
        logger::log_info(" --> ", split_read_counts %>% nrow(), " split read counts loaded from '", cluster_id, "' cluster!")
        
      
        
        ## Load samples
        samples <- readRDS(file = paste0(results_folder_local, "/base_data/", 
                                         project_id, "_", cluster_id, "_samples_used.rds"))
        
        
        if ( !identical(names(split_read_counts)[-1] %>% sort(), samples %>% sort()) ) {
          logger::log_info("The number of samples used does not correspond to the number of columns in the 'split_read_counts' object!")
          break;
        }
        
        
        ## Add coverage detected for the introns in the current tissue
        split_read_counts_intron <- GenerateCoverage(split_read_counts = split_read_counts,
                                                      samples = samples,
                                                      junIDs = split_read_counts$junID) %>%
          dplyr::rename(ref_n_individuals = n_individuals,
                        ref_sum_counts = sum_counts) %>% 
          rowwise() %>%
          mutate(junID = ifelse(str_detect(string = junID, pattern = "\\*"), 
                                    str_replace(string = junID, pattern = "\\*", strand ),
                                    junID)) 
        
        split_read_counts_intron %>% head()
        
        
        
        
        ## Merge counts and split reads from novel combos
        db_introns <- split_read_counts_intron %>%
          inner_join(y = all_split_reads ,
                     by = "junID")
        
        
        ############################################
        ## ADD THE TRANSCRIPT FOREING KEY
        ############################################
        
        logger::log_info(" --> adding TRANSCRIPT foreing key to the novel combos...")
     
        ## Add the GENE ID for the foreign key
        db_introns_tidy <- db_introns %>% 
          dplyr::select(-gene_id) %>%
          unnest(tx_id_junction) %>% 
          left_join(master_transcript ,
                    by = c("tx_id_junction" = "transcript_id")) %>%
          dplyr::select(-gene_id, -tx_id_junction, -TSL, -MANE) %>%
          dplyr::rename(transcript_id = id)  %>%
          drop_na(transcript_id) %>%
          group_by(junID) %>%
          mutate(transcript_id_list = paste(transcript_id, collapse = ",")) %>%
          ungroup() %>%
          distinct(junID, .keep_all = T)
        
        rm(db_introns)
        
        logger::log_info(db_introns_tidy %>% distinct(junID) %>% nrow(), " combos to store....")
        
        summary(db_introns_tidy)
        
        #####################################
        ## GET THE GENE TPM
        #####################################
        
        
        if ( file.exists( paste0(results.folder, "/", project_id, "/tpm/", project_id, "_", cluster_id, "_tpm.rds")) ) {
          
          logger::log_info(" --> calculating TPM values ... ")
          
          tpm <- readRDS(file = paste0(results.folder,  "/", project_id, "/tpm/", project_id, "_", cluster_id, "_tpm.rds")) %>% 
            dplyr::select(gene_id, all_of(samples))
          
          tpm <- tpm  %>%
            dplyr::mutate(tpm_median = matrixStats::rowMedians(x = as.matrix(.[2:(ncol(tpm))]))) %>%
            dplyr::select(gene_id, tpm_median) 
          
          ## In case there are any duplicates, take the genes with the maximum tpms
          tpm <- tpm %>% 
            distinct(gene_id, .keep_all = T) %>%
            group_by(gene_id) %>% 
            summarize_all(max) %>%
            ungroup()
          
          tpm %>% as_tibble()

          tpm_tidy <- tpm %>%
            inner_join(y = master_gene %>% as_tibble(),
                       by = c("gene_id" = "gene_id")) %>%
            inner_join(y = master_transcript %>% as_tibble(),
                       by = c("id" = "gene_id"),
                       multiple = "all") %>%
            dplyr::select(transcript_id = id.y, 
                          tpm_median)
          
          db_introns_tidy <- db_introns_tidy %>%
            left_join(y = tpm_tidy,
                      by = c("transcript_id" = "transcript_id")) %>% 
            dplyr::rename(gene_tpm = tpm_median)
          
        }
        
        
        ######################################
        ## INTRONS - ADD MAXENTSCAN INFO 
        ######################################
        
        logger::log_info(" --> adding the MaxEntScan info ...")
        
        wd <- getwd()
        if ( !file.exists(paste0(dependencies.folder, "/Homo_sapiens.GRCh38.dna.primary_assembly.fa")) ) {
          stop("ERROR! File dependency 'Homo_sapiens.GRCh38.dna.primary_assembly.fa' does not exist within the specified dependencies folder.")
        }
        ## Add MaxEntScan score to the split reads
        db_introns_tidy <- GenerateMaxEntScore(junc_tidy = db_introns_tidy ,
                                               max_ent_tool_path = paste0(dependencies.folder, "/fordownload/"),
                                               homo_sapiens_fasta_path = paste0(dependencies.folder, "/Homo_sapiens.GRCh38.dna.primary_assembly.fa") ) %>% as_tibble()
        
        db_introns_tidy <- db_introns_tidy %>% 
          dplyr::select(-c(donorSeqStart, donorSeqStop, AcceptorSeqStart, AcceptorSeqStop)) %>%
          dplyr::rename(ref_donor_sequence = donor_sequence, ref_acceptor_sequence = acceptor_sequence)
        
        setwd(wd)
   
        db_introns_tidy <- db_introns_tidy %>%
          dplyr::rename(ref_mes5ss = ss5score, ref_mes3ss = ss3score) %>%
          dplyr::relocate(c(ref_mes5ss, ref_mes3ss), .before = transcript_id ) %>% 
          as_tibble()
        
        
        ################################################
        ## INTRONS - ADD THE CONSERVATION AND CDTS INFO
        ################################################
        
        db_introns_tidy %>% as_tibble()
        
        logger::log_info("adding CDTS and Conservation scores...")
        
        df_all_introns <- GenerateCdtsPhastconsScores(db_introns = db_introns_tidy, intron_size = 100, phastcons_type = 17) %>% as_tibble()
        df_all_introns <- df_all_introns %>% mutate_if(is.numeric, ~replace_na(., 0))
        
        df_all_introns %>% as_tibble()
        rm(db_introns_tidy)
        
        
        #########################################
        ## INTRONS - ADD THE TRANSCRIPT BIOTYPE
        #########################################
        
        df_biotype_junID <- readRDS(file = file.path(database.folder, "all_split_reads_qc_level1_PC_biotype.rds")) %>% as_tibble()
        
        any(str_detect(df_biotype_junID$junID, pattern = "\\*"))
        
        df_all_introns <- df_all_introns %>% inner_join(y = df_biotype_junID %>% dplyr::select(junID, protein_coding), by = c("junID"))
        
        logger::log_info(df_all_introns$junID %>% unique %>% length(), " introns to be stored!")
        

        #########################################################
        ## CREATE AND POPULATE CHILD 'NOVEL COMBO' TABLE
        #########################################################
        
        logger::log_info( " --> creating 'mis-spliced' table ... ")
        
        # dbRemoveTable(conn = con, paste0(cluster_id, "_", project_id))
        query <- paste0("CREATE TABLE IF NOT EXISTS '", paste0(cluster_id, "_", project_id, "_combo"), "'", 
                  "(junID NUMERIC PRIMARY KEY NOT NULL,
                  coordinates TEXT NOT NULL, 

                  seqnames TEXT NOT NULL,
                  start NUMERIC NOT NULL,
                  end NUMERIC NOT NULL,
                  strand TEXT NOT NULL, 

                  ref_n_individuals INTEGER NOT NULL,
                  ref_sum_counts INTEGER NOT NULL,
                  type TEXT NOT NULL, 

                  left_motif TEXT NOT NULL, 
                  right_motif TEXT NOT NULL, 

                  ref_length INTEGER NOT NULL, 
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
                  transcript_id_list TEXT NOT NULL,

                  FOREIGN KEY (transcript_id) REFERENCES 'transcript'(id))")
        
        ## Connect the database
        con <- dbConnect(RSQLite::SQLite(), database.path)
        DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
        ## Create the child table
        res <- DBI::dbSendQuery(conn = con, statement = query)
        DBI::dbClearResult(res)
        
        
        ## POPULATE THE TABLE
        df_all_introns_tidy_final <- df_all_introns %>%
          dplyr::rename(coordinates = junID,
                        ref_length = width) %>%
          tibble::rowid_to_column("junID") %>% 
          mutate(transcript_id = transcript_id %>% as.integer())
        
        summary(df_all_introns_tidy_final)
        
        DBI::dbAppendTable(conn = con,
                           name = paste0(cluster_id, "_", project_id, "_combo"), 
                           value = df_all_introns_tidy_final )
          
        
        ## Disconnect the database
        DBI::dbDisconnect(conn = con)
        
        
        logger::log_info("'",cluster_id, "_", project_id, "_combo' table populated!")
          
          
          
      } else {
        DBI::dbDisconnect(conn = con) 
        logger::log_info("Dependency files do not exist! '", cluster_id, "_", project_id, "'!")
      }
    }
    gc()
  }
}
