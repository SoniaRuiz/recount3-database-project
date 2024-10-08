
#' Title
#' Creates two 'child' SQLITE tables per sample cluster to store i) the annotated introns that do not have evidence for mis-splicing activity 
#' (never-misspliced table) and the ii) annotated introns with evidence of mis-splicing at their donor or acceptor splice sites.
#' These are referred to as 'child' tables because they inherit information from the 'intron' and 'novel' master tables
#' @param database.path Local path to the .sqlite database
#' @param recount3.project.IDs Vector with the ID of the recount3 projects to work with
#' @param database.folder Local path to the folder that contains the database
#' @param results.folder Local path to the folder that contains the result files
#'
#' @return
#' @export
#'
#' @examples
SqlCreateChildTables <- function(database.path,
                                 database.folder,
                                 results.folder,
                                 recount3.project.IDs = NULL) {
  
  
  ###########################
  ## LOAD DATA DEPENDENCIES
  ###########################
  
  if ( !file.exists(paste0(database.folder, "/all_split_reads_qc_level2.rds")) ) {
    logger::log_info("ERROR! Second-filter level of split reads file not found!")
    break;
  } else {
    logger::log_info("Loading second-filter level of split reads file ...")
    all_split_reads_details <- readRDS(file = paste0(database.folder, "/all_split_reads_qc_level2.rds")) %>% dplyr::select(junID, strand)
  }
  
  ###########################
  ## QUERY MASTER TABLES
  ###########################
  
  con <- dbConnect(RSQLite::SQLite(), database.path)
  tables <- DBI::dbListTables(conn = con)
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  
  logger::log_info("SQL connection stablished!") 
  logger::log_info("Querying master tables ...")
  
  ## GET FROM MASTER TABLE
  query = paste0("SELECT * FROM 'metadata'")
  master_metadata <- dbGetQuery(con, query) 
  
  ## GET FROM INTRON TABLE
  query = paste0("SELECT ref_junID, ref_coordinates, transcript_id, misspliced FROM 'intron'")
  master_intron <- dbGetQuery(con, query) %>% as_tibble()
  master_intron %>% nrow()
  master_intron %>% 
    dplyr::count(misspliced)
  
  ## GET FROM NOVEL JUNCTION TABLE
  query = paste0("SELECT ref_junID, novel_junID, novel_coordinates, novel_type FROM 'novel'")
  master_novel <- dbGetQuery(con, query) %>% as_tibble()
  master_novel %>% nrow()
  master_novel %>% 
    dplyr::count(novel_type) %>%
    print()
  
  setdiff(master_novel$ref_junID, master_intron$ref_junID)
  
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
    recount3.project.IDs <- (master_metadata$SRA_project %>% unique())
  }
  
  ## Loop through the projects parallely
  # doParallel::registerDoParallel(10)
  # foreach(i = seq(length(recount3.project.IDs))) %dopar%{
  #   
  #   project_id <- recount3.project.IDs[i]
  
  
  ###########################
  ## CREATE CHILD TABLES
  ###########################
  
  
  for ( project_id in recount3.project.IDs ) { 
    
    # project_id <- recount3.project.IDs[1]
    # project_id <- recount3.project.IDs[2]
    # project_id <- "LAML"
    # project_id <- "COAD" 
    
    logger::log_info(" --> Initiating data process of '", project_id, "' ...")
    results_folder_local <- paste0(results.folder, "/", project_id, "/")
    
    clusters <- master_metadata %>%
      dplyr::filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    for (cluster_id in clusters) { 
      
      # cluster_id <- clusters[1]
      # cluster_id <- clusters[3]
      
      logger::log_info(project_id, " --> ", cluster_id)
      
      
      
      ###############################
      ## PREPARE DATA
      ###############################
      con <- dbConnect(RSQLite::SQLite(), database.path)
      
      if ( !DBI::dbExistsTable(con, name = paste0(cluster_id, "_", project_id, "_misspliced")) &&
           !DBI::dbExistsTable(con, name = paste0(cluster_id, "_", project_id, "_nevermisspliced")) && 
           file.exists(paste0(results_folder_local, "/junction_pairing/", 
                              cluster_id, "/", cluster_id, "_raw_distances_tidy.rds")) ) {
        
        DBI::dbDisconnect(conn = con) 
        
        #########################################################
        ## LOAD BASE DATA ONLY FOR THE CURRENT CLUSTER ID 
        #########################################################
        
        logger::log_info(cluster_id, " loading base data ... ")
        
        ## Load split read counts
        split_read_counts <- readRDS(file = paste0(results_folder_local, "/base_data/", 
                                                   project_id, "_", cluster_id, "_split_read_counts.rds")) 
        
        ## Load samples
        samples <- readRDS(file = paste0(results_folder_local, "/base_data/", 
                                         project_id, "_", cluster_id, "_samples_used.rds"))
        
        ## LOAD INTRONS AND NOVEL JUNCTIONS paired for the current tissue 
        df_cluster_distances <- readRDS(file = paste0(results_folder_local, "/junction_pairing/", cluster_id, "/", 
                                                      cluster_id, "_raw_distances_tidy.rds")) %>% as_tibble()
        
        
        if ( is.null(names(split_read_counts)) ) {
          split_read_counts <- split_read_counts %>%
            as_tibble(rownames = "junID")
        }
        if ( !identical(names(split_read_counts)[-1] %>% sort(), samples %>% sort()) ) {
          logger::log_info("The number of samples used does not correspond to the number of columns in the 'split_read_counts' object!")
          break;
        }
        
        
       
        #########################################################
        ## ADD EXPRESSION DATA TO THE PAIRED JUNCTIONS
        #########################################################
       
        
        ## Add total number of supporting split reads across samples to introns and novel junctions
        df_local_intron_pairings_w_counts <- add_coverage_to_introns(df_cluster_distances, split_read_counts, samples)
        df_local_novel_pairings_w_counts <- add_coverage_to_novel_jxn(df_cluster_distances, split_read_counts, samples)

        
        
        ## QC - remove star from IDs
        df_local_intron_pairings_w_counts <- qc_replace_star_ID(db.introns = df_local_intron_pairings_w_counts)
        df_local_novel_pairings_w_counts <- qc_replace_star_ID(db.introns = df_local_novel_pairings_w_counts)
       
        
        
        ## Merge introns and novel junctions with the coverage data added
        logger::log_info("merging introns and novel junctions ... ")
        df_local_intron_w_novel_pairings <- df_local_novel_pairings_w_counts %>%
          dplyr::select(-c(seqnames, start, end, strand)) %>%
          inner_join(y = df_local_intron_pairings_w_counts %>% dplyr::select(-c(seqnames, start, end, strand)),
                     by = "ref_junID")
        
        
        
        ## QC
        if (!identical(df_local_intron_w_novel_pairings %>% nrow(),
                       df_cluster_distances %>% nrow())) {
          message("ERROR: some junctions have been lost in the process of adding the coverage.")
          break;
        }
        if (any(str_detect(df_local_intron_w_novel_pairings$ref_junID,pattern = "//*"))) {
          logger::log_info("ERROR: some IDs contain '*'")
          break;
        }
        if (any(str_detect(df_local_intron_w_novel_pairings$novel_junID,pattern = "//*"))) {
          logger::log_info("ERROR: some IDs contain '*'")
          break;
        }
        if (any(df_local_intron_w_novel_pairings$ref_junID %>% is.na())) {
          logger::log_info("ERROR: There are missing reference introns!!")
          break;
        }

        
        
        ################################################
        ## ADD FOREIGN KEY FROM THE MASTER TABLES
        ################################################
        
        logger::log_info("adding MASTER INTRON foreign key... ")
        
        ## JOIN data with MASTER INTRON table
        df_local_intron_w_novel_pairings_w_master <- df_local_intron_w_novel_pairings %>% 
          left_join(y = master_intron %>% 
                      dplyr::filter(misspliced == "Yes") %>%
                      dplyr::select(ref_junID, ref_coordinates, transcript_id),
                    by = c("ref_junID" = "ref_coordinates")) %>%
          dplyr::rename(ref_coordinates = ref_junID) %>%
          dplyr::rename(ref_junID = ref_junID.y) %>%
          dplyr::relocate(ref_junID)
        
        logger::log_info("adding MASTER NOVEL foreign key... ")
        
        df_local_intron_w_novel_pairings_w_master <- df_local_intron_w_novel_pairings_w_master %>% 
          inner_join(y = master_novel,
                     by = c("novel_junID" = "novel_coordinates",
                            "ref_junID" = "ref_junID")) %>%
          dplyr::rename(novel_coordinates = novel_junID) %>%
          dplyr::rename(novel_junID = novel_junID.y) %>%
          dplyr::relocate(ref_junID, novel_junID)
        
        
        
        
        # logger::log_info(df_local_intron_w_novel_pairings_w_master %>% nrow(), " non-ambiguous junction pairings to be stored... ")
       
        
        if (setdiff(df_local_intron_w_novel_pairings_w_master$ref_coordinates, 
                    master_intron$ref_coordinates) %>% length() > 0) {
          logger::log_info("ERROR! Some introns detected in this tissue are not stored in the master intron table.")
          break;
        }
        if (which(str_detect(df_local_intron_w_novel_pairings_w_master$novel_coordinates,pattern = "//*")) %>% length() > 0) {
          logger::log_info("ERROR: some nove_junIDs contain '*'")
          break;
        }
        
        
        ####################################
        ## LOAD INTRONS PARENTING AMBIGUOUS 
        ## NOVEL JXN
        ####################################
        
      
        ## Novel junctions are ambiguous and we do not include them in the database.
        ## However, we do store annotated introns parenting them that have passed the QC (i.e. stored on the 'intron' master database)
        
        introns_parent_ambigous <- intersect(master_intron$ref_coordinates, 
                                             setdiff(df_local_intron_w_novel_pairings$ref_junID,
                                                     df_local_intron_w_novel_pairings_w_master$ref_coordinates) %>% unique)
        
        
        ## TYPE 'MAYBE'
        introns_parent_ambigous <- data.frame(ref_junID = introns_parent_ambigous) %>% 
          as_tibble() %>%
          mutate(ref_type = "maybe")
        
        # missing_intron <- "chr1:155240078-155244505:-"
        # df_local_intron_pairings %>% filter(ref_junID == missing_intron)
        # split_read_counts_intron %>% filter(junID == missing_intron)
        # df_local_intron_pairings_w_counts %>% filter(ref_junID == missing_intron)
        # df_local_novel_pairings_w_counts %>% filter(ref_junID == missing_intron)
        # introns_parent_ambigous %>% filter(ref_junID == missing_intron)
        
        logger::log_info(introns_parent_ambigous %>% distinct(ref_junID) %>% nrow(), " introns parenting ambiguous novel junctions that also passed the QC...")
        
        
        #######################################
        ## CHECK INTEGRITY WITH PARENT TABLES
        #######################################
        
        
        logger::log_info("checking integrity with parent tables ... ")
        
        df_local_novel <- master_novel %>%
          dplyr::filter(novel_junID %in% (df_local_intron_w_novel_pairings_w_master %>% pull(novel_junID))) %>% 
          dplyr::select(novel_coordinates) 
        
        if ( !(identical(df_local_intron_w_novel_pairings_w_master$novel_coordinates %>% sort(), 
                         df_local_novel$novel_coordinates %>% sort())) ) {
          
          logger::log_info("ERROR! Local mis-splicing tables and master novel table are not identical")
          
          if ( all(intersect(setdiff(df_local_intron_w_novel_pairings_w_master$novel_coordinates, 
                                     master_novel$novel_coordinates),
                             ambiguous_novel_junc$novel_junID) == setdiff(df_local_intron_w_novel_pairings_w_master$novel_coordinates, 
                                                                          master_novel$novel_coordinates)) == T ) {
            df_local_intron_w_novel_pairings_w_master <- df_local_intron_w_novel_pairings_w_master %>%
              as.data.table() %>%
              dplyr::filter(!(novel_coordinates %in% ambiguous_novel_junc$novel_junID))
          }
          
          break;
          
        } else {
          
          
       
          
          
          df_local_intron_w_novel_pairings_w_master <- df_local_intron_w_novel_pairings_w_master %>%
            dplyr::select(-novel_coordinates)
          
          #####################################
          ## QC - CHECK INTEGRITY WITH PARENT TABLE
          #####################################
         
          df <- master_novel %>%
            dplyr::select(novel_junID, ref_junID) %>%
            arrange(novel_junID) %>% 
            inner_join(df_local_intron_w_novel_pairings_w_master %>%
                         dplyr::select(novel_junID, ref_junID) %>%
                         arrange(novel_junID),
                       by = "novel_junID")
          
          diff <- df %>% 
            dplyr::filter(ref_junID.x != ref_junID.y)
          
          
          if (diff %>% nrow() > 0) {
            
            df_local_intron_w_novel_pairings_w_master <- df_local_intron_w_novel_pairings_w_master %>%
              dplyr::filter(!(ref_junID %in% diff$ref_junID.y)) 
            logger::log_info("ERROR!: ", diff, " --> mismatch junctions.")
            break;
            
          } 
     
          
          #####################################
          ## CALCULATE MSR MEASURES
          #####################################
          
          
          df_local_pairings_w_master_w_MSR <- add_MSR_measures(db.introns = df_local_intron_w_novel_pairings_w_master)
          
          
          #####################################
          ## GET THE GENE TPM
          #####################################
          
          
          df_local_pairings_w_master_w_MSR_w_TPM <- add_median_TPM_values(results.folder = results.folder, 
                                                                          cluster.samples = samples, 
                                                                          master.gene = master_gene,
                                                                          master.transcript = master_transcript,
                                                                          project.id = project_id, 
                                                                          cluster.id = cluster_id, 
                                                                          db.introns = df_local_pairings_w_master_w_MSR)
          
          
          #####################################
          ## ADD THE TYPE OF INTRON
          #####################################
          
          
          df_local_pairings_w_master_w_MSR_w_TPM <- add_intron_category(db.introns = df_local_pairings_w_master_w_MSR_w_TPM)
          
          
          
          #########################################################
          ## CREATE AND POPULATE CHILD 'MIS-SPLICED' INTRON TABLE
          #########################################################
          
          
          db_introns_final <- df_local_pairings_w_master_w_MSR_w_TPM %>%
            dplyr::select(-novel_acceptor,-novel_donor,-ref_coordinates)%>%
            dplyr::rename(MSR_D = MSR_Donor, MSR_A = MSR_Acceptor)
          
          create_and_populate_misspliced_child_table(database.path, 
                                                     cluster.id = cluster_id, 
                                                     project.id = project_id, 
                                                     db.introns.final = db_introns_final)
          
          
  
          
          ####################################
          ## LOAD NON PAIRED INTRONS
          ####################################
          
          # These are the annotated introns that could not be paired
          ## These introns are not ambiguous and should be included with A MAYBE evidence of mis-splicing
          ## LOAD NON PAIRED INTRONS
          df_non_paired_introns <- readRDS(file = paste0(results_folder_local, "/junction_pairing/", cluster_id, "/not-misspliced/", 
                                                         cluster_id, "_all_misspliced_not_paired.rds"))
          never_additional_introns <- intersect(master_intron$ref_coordinates, df_non_paired_introns)
          logger::log_info("Additional never introns to be added: ", 
                           never_additional_introns %>% unique %>% length())
          
          
          ## TYPE 'MAYBE'
          never_additional_introns <- data.frame(ref_junID = never_additional_introns) %>% 
            as_tibble() %>%
            mutate(ref_type = "maybe")
          
          logger::log_info(never_additional_introns %>% distinct(ref_junID) %>% nrow(), " non-paried introns passing QC - type 'maybe'.")
          
          
          
     
          
          ####################################
          ## NEVER MISSPLICED INTRONS
          ####################################
          
          logger::log_info("getting never mis-spliced introns ... ")
          
          
          ## LOAD NEVERMIS-SPLICED INTRONS - INTRON TYPE 'NONE'
          introns_never <- data.frame(ref_junID = readRDS(file = paste0(results_folder_local, "/junction_pairing/", cluster_id, "/not-misspliced/", 
                                                                        cluster_id, "_all_notmisspliced.rds"))) %>% as_tibble() %>% mutate(ref_type = "never")
          
          logger::log_info(introns_never %>% distinct(ref_junID) %>% nrow(), " 'never' type introns...")
          
          
          ## JOIN NEVER MISSPLICED AND MAYBE NOT MISSPLICED
          
          df_introns_never <- rbind(introns_never, never_additional_introns, introns_parent_ambigous)  %>%
            inner_join(y = all_split_reads_details,
                       by = c("ref_junID" = "junID")) %>%
            mutate(strand = strand %>% as.character()) %>%
            rowwise() %>%
            mutate(ref_junID = ifelse(str_detect(string = ref_junID, pattern = "\\*"), 
                                      str_replace(string = ref_junID, pattern = "\\*", strand ),
                                      ref_junID)) %>%
            dplyr::select(-strand)
          
          
         
          
          ## The introns not misspliced in this tissue, should have not been detected as spliced.
          ## Thus, this should be zero
          if (intersect(df_introns_never,
                        master_intron %>% dplyr::filter(ref_junID %in% db_introns_final$ref_junID) %>% pull(ref_coordinates) %>% unique()) %>% 
              length() > 0 ){
            logger::log_info("ERROR! Some never-misspliced introns have also been classified as mis-spliced.")
            break;
          }
          
          if (any(str_detect(string = df_introns_never$ref_junID, pattern = "\\*"))) {
            logger::log_info("ERROR! * in the IDs")
            break;
          }
          
          
          
          split_read_counts_intron_never <- generate_coverage(split_read_counts = split_read_counts,
                                                              samples = samples,
                                                              junID = df_introns_never$ref_junID) %>%
            dplyr::rename(ref_n_individuals = n_individuals, ref_sum_counts = sum_counts) 
          
          
          if (any(str_detect(split_read_counts_intron_never$junID,pattern = "\\*"))) {
            logger::log_info("ERROR! some never mis-spliced junctions without the number of individuals")
            break;
          }
          
          
          #####################################
          ## CALCULATE MSR MEASURES
          #####################################
          
          
          db_never_introns_final <- df_introns_never %>%
            inner_join(y = split_read_counts_intron_never, by = c("ref_junID" = "junID")) %>%
            mutate(MSR_D = 0, MSR_A = 0) %>% as_tibble()
          
          
          ## QC
          if ( any(db_never_introns_final$ref_n_individuals %>% is.na()) ) {
            logger::log_info("ERROR! some never mis-spliced junctions without the number of individuals")
            break;
          }
          if (any(str_detect(db_never_introns_final$ref_junID, pattern = "\\*"))) {
            logger::log_info("ERROR! * in the IDs")
            break;
          }
        
          
          ## TYPE 'MAYBE'
          logger::log_info("Junctions type 'maybe': ", db_never_introns_final %>% filter(ref_type == "maybe") %>% distinct(ref_junID) %>% nrow())
          
          ## TYPE 'NONE'
          logger::log_info("Junctions type 'never': ", db_never_introns_final %>% filter(ref_type == "never") %>% distinct(ref_junID) %>% nrow())
          
          
          ##################################################
          ## ADD REFERENCE KEY TO THE MASTER INTRON TABLE 
          ##################################################
          
          logger::log_info( "adding the intron reference key to the never mis-spliced introns ... ")
          
          db_never_introns_final <- db_never_introns_final %>%
            inner_join(master_intron %>% dplyr::select(ref_junID, ref_coordinates, transcript_id),
                       by = c("ref_junID" = "ref_coordinates")) %>%
            dplyr::filter(!is.na(ref_junID)) %>%
            dplyr::select(-ref_junID) %>% 
            dplyr::rename(ref_junID = ref_junID.y) %>%
            relocate(ref_junID)
          
          
          if (db_never_introns_final %>% dplyr::filter(is.na(ref_junID)) %>% nrow > 0) {
            logger::log_info("ERROR! IDs are NA")
            break;
          }
          if ((intersect(db_never_introns_final$ref_junID, db_introns_final$ref_junID) %>% length()) > 0) {
            logger::log_info("Error! Some never mis-spliced junctions have been stored as mis-spliced.")
            break;
          }
          if (any(duplicated(db_never_introns_final$ref_junID))) {
            logger::log_info("Error! Some never mis-spliced junctions are duplicated.")
            break;
          }
          
          
          #####################################
          ## ADD GENE TPM
          #####################################

          db_never_introns_final <- add_median_TPM_values(results.folder = results.folder, 
                                                          cluster.samples = samples,
                                                          master.gene = master_gene,
                                                          master.transcript = master_transcript,
                                                          project.id = project_id, 
                                                          cluster.id = cluster_id, 
                                                          db.introns = db_never_introns_final)
          
          #####################################
          ## QC
          #####################################
          
          if (intersect(db_never_introns_final %>% dplyr::filter(ref_type == "never") %>% distinct(ref_junID) %>% pull(), 
                        db_introns_final %>% dplyr::filter(ref_type == "donor") %>% distinct(ref_junID) %>% pull()) %>% length() > 0) {
            logger::log_info("Error! Some never-misspliced introns are classified as mis-spliced as the donor splice site!")
            break;
          }
          
          
          if (intersect(db_never_introns_final %>% dplyr::filter(ref_type == "never") %>% distinct(ref_junID) %>% pull(), 
                        db_introns_final %>% dplyr::filter(ref_type == "aceptor") %>% distinct(ref_junID) %>% pull()) %>% length() > 0) {
            logger::log_info("Error! Some never-misspliced introns are classified as mis-spliced as the acceptor splice site!")
            break;
          }
          
          if (intersect(db_never_introns_final %>% dplyr::filter(ref_type == "never") %>% distinct(ref_junID) %>% pull(), 
                        db_introns_final %>% dplyr::filter(ref_type == "both") %>% distinct(ref_junID) %>% pull()) %>% length() > 0) {
            logger::log_info("Error! Some never-misspliced introns are classified as mis-spliced as both splice sites!")
            break;
          }
          
          if (intersect(db_never_introns_final %>% dplyr::filter(ref_type == "maybe") %>% distinct(ref_junID) %>% pull(), 
                        db_never_introns_final %>% dplyr::filter(ref_type == "never") %>% distinct(ref_junID) %>% pull()) %>% length() > 0) {
            logger::log_info("Error! Some never-misspliced introns are classified as maybe mis-spliced!")
            break;
          }
          
          if (intersect(db_never_introns_final %>% dplyr::filter(ref_type == "maybe") %>% distinct(ref_junID) %>% pull(), 
                        db_introns_final %>% distinct(ref_junID) %>% pull()) %>% length() > 0) {
            logger::log_info("Error! Some maybe misspliced introns are classified as mis-spliced!")
            break;
          }
          
          if (intersect(db_introns_final %>% dplyr::filter(ref_type == "acceptor") %>% distinct(ref_junID) %>% pull(), 
                        db_introns_final %>% dplyr::filter(ref_type == "donor") %>% distinct(ref_junID) %>% pull()) %>% length() > 0) {
            logger::log_info("Error! Some misspliced introns only at the acceptor are also classified as mis-spliced as the donor splice site!")
            break;
          }
          
          
          
          
          
          #############################################################
          ## CREATE AND POPULATE CHILD 'NEVER MIS-SPLICED' INTRON TABLE
          #############################################################
          
          create_and_populate_never_misspliced_child_table(database.path, 
                                                           cluster.id = cluster_id, 
                                                           project.id = project_id, 
                                                           db.introns.final = db_never_introns_final)

          
        } 
        
        
        rm(df_local_intron_pairings_w_counts)
        rm(df_local_novel_pairings_w_counts)
        rm(df_local_intron_w_novel_pairings)
        rm(df_local_intron_w_novel_pairings_w_master)
        rm(introns_parent_ambigous)
        rm(split_read_counts)
        rm(df_cluster_distances)
        rm(samples)
        rm(df_local_novel)
        rm(df_local_pairings_w_master_w_MSR)
        rm(df_local_pairings_w_master_w_MSR_w_TPM)
        rm(db_introns_final)
        rm(df_non_paired_introns)
        rm(never_additional_introns)
        rm(df_introns_never)
        rm(split_read_counts_intron_never)
        rm(db_never_introns_final)
        gc()
        
      } else {
        DBI::dbDisconnect(conn = con) 
        logger::log_info("Tables '", cluster_id, "_", project_id, "_nevermisspliced' and '", cluster_id, "_", project_id, "_misspliced' exist!")
      }
    }
    
  }
}


#' Title
#' Creates one child SQLITE table per sample cluster to store the 'novel combination' splicing events.
#' This table is referred to as 'child' table because it inherit information from the 'combo' master table
#' @return
#' @export
#'
#' @examples
SqlCreateChildTableCombo <-  function(con = NULL) {
  
  ## CONNECT THE DATABASE
  if (is.null(con)) {
    con <- dbConnect(RSQLite::SQLite(), database.path)
  }
  
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  
  tables <- DBI::dbListTables(conn = con)
  
  ## GET INFO FROM COMBO MASTER TABLE
  query <- paste0("SELECT ref_junID, ref_coordinates, transcript_id FROM 'combo'")
  master_combo <- dbGetQuery(con, query) %>% as_tibble()
  master_combo %>% nrow()
  
  
  for ( project_id in recount3.project.IDs ) { 
    
    # project_id <- recount3.project.IDs[1]
    
    logger::log_info(" --> Working with '", project_id, "' ...")
    results_folder_local <- file.path(results.folder, project_id)
    
    clusters <- df_metadata %>% dplyr::filter(SRA_project == project_id) %>% distinct(cluster) %>% pull()
    
    for (cluster_id in clusters) { 
      
      # cluster_id <- clusters[1]
      
      logger::log_info(project_id, " --> ", cluster_id)
      
      ###############################
      ## PREPARE DATA
      ###############################
      
      if ( file.exists(paste0(results_folder_local, "/base_data/", project_id, "_", cluster_id, "_all_split_reads_combos.rds")) && 
           file.exists(paste0(results_folder_local, "/base_data/", project_id, "_", cluster_id, "_split_read_counts_combos.rds")) ) {  
        
        ## Load split read counts
        split_read_counts <- readRDS(file = paste0(results_folder_local, "/base_data/", project_id, "_", cluster_id, "_split_read_counts_combos.rds"))
        
        logger::log_info(" --> ", split_read_counts %>% nrow(), " split read counts loaded from '", cluster_id, "' cluster!")
        
        ## Load samples
        samples <- readRDS(file = paste0(results_folder_local, "/base_data/", project_id, "_", cluster_id, "_samples_used.rds"))
        
        if ( !identical(names(split_read_counts)[-1] %>% sort(), samples %>% sort()) ) {
          stop("The number of samples used does not correspond to the number of columns in the 'split_read_counts' object!")
        }
        
        ## Add coverage detected for the introns in the current tissue
        split_read_counts_w_coverage <- GenerateCoverage(split.read.counts = split_read_counts,
                                                         samples = samples,
                                                         junIDs = split_read_counts$junID) %>%
          dplyr::rename(ref_n_individuals = n_individuals, ref_sum_counts = sum_counts) %>% 
          rowwise() %>%
          mutate(junID = ifelse(str_detect(string = junID, pattern = "\\*"), str_replace(string = junID, pattern = "\\*", strand ),
                                junID)) 
        
        split_read_counts_w_coverage %>% head()
        
        ## Merge counts and split reads from novel combos
        split_read_counts_all_details <- split_read_counts_w_coverage %>%
          inner_join(y = master_combo, by = c("junID" = "ref_coordinates"))
        
        #####################################
        ## AD GENE TPM
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
          
          split_read_counts_all_details_tidy <- split_read_counts_all_details_tidy %>%
            left_join(y = tpm_tidy,
                      by = c("transcript_id" = "transcript_id")) %>% 
            dplyr::rename(gene_tpm = tpm_median)
          
        }
        
        #########################################################
        ## CREATE AND POPULATE CHILD 'NOVEL COMBO' TABLE
        #########################################################
        
        logger::log_info( " --> creating 'novel combo' table ... ")
        
        # dbRemoveTable(conn = con, paste0(cluster_id, "_", project_id))
        query <- paste0("CREATE TABLE IF NOT EXISTS '", paste0(cluster_id, "_", project_id, "_combo"), "'", 
                        "(ref_junID INTEGER NOT NULL,
                          ref_n_individuals INTEGER NOT NULL,
                          ref_sum_counts INTEGER NOT NULL,
                          gene_tpm DOUBLE, 
                          transcript_id INTEGER NOT NULL, 
                          FOREIGN KEY (ref_junID) REFERENCES combo (ref_junID),
                          FOREIGN KEY (transcript_id) REFERENCES 'transcript'(id))")
        
        ## Connect the database
        con <- dbConnect(RSQLite::SQLite(), database.path)
        DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
        
        ## Create the NOVEL COMBO table
        res <- DBI::dbSendQuery(conn = con, statement = query)
        DBI::dbClearResult(res)
        
        
        ## POPULATE THE TABLE
        split_read_counts_all_details_final <- split_read_counts_all_details %>%
          dplyr::select(-junID) %>%
          mutate(transcript_id = transcript_id %>% as.integer()) %>%
          dplyr::relocate(ref_junID)
        
        summary(split_read_counts_all_details_final)
        
        DBI::dbAppendTable(conn = con, name = paste0(cluster_id, "_", project_id, "_combo"), 
                           value = split_read_counts_all_details_final )
        
        
        ## Disconnect the database
        DBI::dbDisconnect(conn = con)
        
        
        logger::log_info("'", cluster_id, "_", project_id, "_combo' table populated!")
        
        
        
      } else {
        DBI::dbDisconnect(conn = con) 
        logger::log_info("Dependency files do not exist! '", cluster_id, "_", project_id, "'!")
      }
    }
    gc()
  }
}



## HELPER FUNCTIONS ------------------------------------------------------------------------------------

add_coverage_to_introns <- function(df_cluster_distances, split_read_counts, samples) {
  
  ## 1. Obtain unique INTRONS from the pairings
  df_local_intron_pairings <- df_cluster_distances %>%
    distinct(ref_junID, .keep_all = T) %>%
    dplyr::select(ref_junID, 
                  seqnames = ref_seq, 
                  start = ref_start,
                  end = ref_end, 
                  strand = ref_strand)
  
  
  ## Add coverage detected for the introns in the current tissue
  split_read_counts_intron <- generate_coverage(split_read_counts = split_read_counts,
                                                samples = samples,
                                                junIDs = df_local_intron_pairings$ref_junID) %>%
    dplyr::rename(ref_n_individuals = n_individuals,
                  ref_sum_counts = sum_counts) %>% as_tibble()
  split_read_counts_intron %>% head()
  df_local_intron_pairings_w_counts <- df_local_intron_pairings %>%
    left_join(y = split_read_counts_intron,
              by = c("ref_junID" = "junID"))
  
  df_local_intron_pairings_w_counts %>% 
    return()
}

add_coverage_to_novel_jxn <- function(df_cluster_distances, split_read_counts, samples) {

  ## Obtain NOVEL JUNCTIONS from the pairings
  df_local_novel_pairings <- df_cluster_distances %>%
    distinct(ref_junID, novel_junID, .keep_all = T) %>%
    dplyr::select(novel_junID, 
                  ref_junID, 
                  seqnames = novel_seq, 
                  start = novel_start, 
                  end = novel_end, 
                  strand = novel_strand)
  
  ## Add coverage detected for the introns in the current tissue
  split_read_counts_novel <- generate_coverage(split_read_counts = split_read_counts,
                                               samples = samples,
                                               junID = df_local_novel_pairings$novel_junID) %>%
    dplyr::rename(novel_n_individuals = n_individuals,
                  novel_sum_counts = sum_counts) %>% as_tibble()
  df_local_novel_pairings_w_counts <- df_local_novel_pairings %>%
    left_join(y = split_read_counts_novel,
              by = c("novel_junID" = "junID")) %>% as_tibble()
  
  
  df_local_novel_pairings_w_counts %>%
    return()

}

qc_replace_star_ID <- function(db.introns) {
  
  
  if ( any(names(db.introns) == "ref_junID") ) {
    db.introns <- db.introns %>% 
      mutate(strand = strand %>% as.character()) %>%
      rowwise() %>%
      mutate(ref_junID = ifelse(str_detect(string = ref_junID, pattern = "\\*"), 
                                str_replace(string = ref_junID, pattern = "\\*", strand  ),
                                ref_junID)) 
  }
  if (any(names(db.introns) == "novel_junID")) {
    db.introns <- db.introns %>% 
      mutate(strand = strand %>% as.character()) %>%
      rowwise() %>%
      mutate(novel_junID = ifelse(str_detect(string = novel_junID, pattern = "\\*"), 
                                str_replace(string = novel_junID, pattern = "\\*", strand ),
                                novel_junID)) 
  }
  
  return(db.introns)
}

add_MSR_measures <- function(db.introns) {
  
  logger::log_info("calculating MSR measures ... ")
  
  db.introns <- db.introns %>%
    group_by(ref_junID, novel_type) %>%
    mutate(MSR = sum(novel_sum_counts)/(sum(novel_sum_counts) + ref_sum_counts)) %>%
    ungroup()
  
  db.introns <- db.introns %>%
    spread(key = novel_type, value = MSR, fill = 0) %>%
    group_by(ref_junID) %>% 
    ## If the annotated intron is mis-spliced at both splice sites, it will have some rows MSR=0 and some rows MSR>0 at the 5' and at the 3'ss. 
    ## We select the highest value
    dplyr::mutate(MSR_Donor = max(novel_donor, na.rm = T)) %>%
    dplyr::mutate(MSR_Acceptor = max(novel_acceptor, na.rm = T)) %>%
    ungroup()
  
  
  db.introns %>% 
    return()
}

add_median_TPM_values <- function(results.folder, cluster.samples, master.gene, master.transcript, project.id, cluster.id, db.introns) {
  
  if ( file.exists( paste0(results.folder, "/", project.id, "/tpm/", project.id, "_", cluster.id, "_tpm.rds")) ) {
    
    logger::log_info("calculating TPM values ... ")
    
    tpm <- readRDS(file = paste0(results.folder,  "/", project.id, "/tpm/", project.id, "_", cluster.id, "_tpm.rds")) %>% 
      dplyr::select(gene_id, all_of(cluster.samples))
    
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
      inner_join(y = master.gene %>% as_tibble(),
                 by = c("gene_id" = "gene_id")) %>%
      inner_join(y = master.transcript %>% as_tibble(),
                 by = c("id" = "gene_id"),
                 multiple = "all") %>%
      dplyr::select(transcript_id = id.y, 
                    tpm_median)
    
    db.introns <- db.introns %>%
      left_join(y = tpm_tidy,
                by = c("transcript_id" = "transcript_id")) %>% 
      dplyr::rename(gene_tpm = tpm_median)
    
    
    
    
  } 
  
  return(db.introns)
}

add_intron_category <- function(db.introns) {
  
  db.introns[, "ref_type"] <- ""
  db.introns %>% head()
  
  ## TYPE 'BOTH'
  indx <- which(db.introns$MSR_Donor > 0 & db.introns$MSR_Acceptor > 0)
  db.introns[indx, "ref_type"] <- "both"
  logger::log_info("Introns type 'both': ",  db.introns[indx, ] %>% distinct(ref_junID) %>% nrow())
  
  
  ## TYPE 'DONOR'
  indx <- which(db.introns$MSR_Donor > 0 & db.introns$MSR_Acceptor == 0)
  db.introns[indx, "ref_type"] <- "donor"
  logger::log_info("Junctions type 'donor': ", db.introns[indx, ] %>% distinct(ref_junID) %>% nrow())
  
  
  ## TYPE 'ACCEPTOR'
  indx <- which(db.introns$MSR_Donor == 0 & db.introns$MSR_Acceptor > 0)
  db.introns[indx, "ref_type"] <- "acceptor"
  logger::log_info("Junctions type 'acceptor': ", db.introns[indx, ] %>% distinct(ref_junID) %>% nrow())
  
  
  # db.introns %>% nrow()
  # db.introns %>% distinct(ref_junID) %>% nrow()
  
  db.introns %>% 
    return()
}

create_and_populate_misspliced_child_table <- function(database.path, cluster.id, project.id, db.introns.final) {
  
  logger::log_info( "creating 'mis-spliced' table ... ")
  
  # dbRemoveTable(conn = con, paste0(cluster.id, "_", project.id))
  query <- paste0("CREATE TABLE IF NOT EXISTS '", paste0(cluster.id, "_", project.id, "_misspliced"), "'", 
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
  
  ## Connect the database
  con <- dbConnect(RSQLite::SQLite(), database.path)
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  ## Create the child table
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  

  DBI::dbAppendTable(conn = con,
                     name = paste0(cluster.id, "_", project.id, "_misspliced"), 
                     value = db.introns.final)
  
  ## CREATE INDEX
  query <- paste0("CREATE UNIQUE INDEX 'index_", paste0(cluster.id, "_", project.id, "_misspliced"), "' ON '",
                  paste0(cluster.id, "_", project.id, "_misspliced"), "'(ref_junID,novel_junID)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  ## Disconnect the database
  DBI::dbDisconnect(conn = con)
  
  
  logger::log_info("'",cluster.id, "_", project.id, "_misspliced' table populated!")
  
}

create_and_populate_never_misspliced_child_table <- function(database.path, cluster.id, project.id, db.introns.final){
  
  query <- paste0("CREATE TABLE IF NOT EXISTS '", paste0(cluster.id, "_", project.id, "_nevermisspliced"), "'", 
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
  
  
  ## Connect the database
  con <- dbConnect(RSQLite::SQLite(), database.path)
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  
  ## Create the child table 
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  ## POPULATE TABLE
  DBI::dbAppendTable(conn = con,
                     name = paste0(cluster.id, "_", project.id, "_nevermisspliced"), 
                     value = db.introns.final)
  
  ## CREATE INDEX
  query <- paste0("CREATE UNIQUE INDEX 'index_", 
                  paste0(cluster.id, "_", project.id, "_nevermisspliced"), "' ON '",
                  paste0(cluster.id, "_", project.id, "_nevermisspliced"),"'(ref_junID)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  ## Disconnect the database
  DBI::dbDisconnect(conn = con)
  
  logger::log_info("'", cluster.id, "_", project.id, "_nevermisspliced' table populated!")
  
}
