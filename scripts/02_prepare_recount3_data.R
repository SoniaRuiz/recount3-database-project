
#' Title
#' Separates the quality-controlled split-read data by sample cluster for a given recount3 project
#' (e.g. by case/control cluster, by tissue, etc. )
#' @param recount3.project.IDs array with the recount3 projects to download
#' e.g. "SRP009615"
#' @param data.source source of the project in recount3. It can be:
#' - "data_sources/sra"  
#' - "data_sources/gtex" 
#' - "data_sources/tcga"
#'
#' @return
#' @export
#'
#' @examples
PrepareRecount3Data <- function(recount3.project.IDs, 
                                data.source,
                                results.folder,
                                levelqc1.folder,
                                supporting.reads,
                                num.cores,
                                subsampling = F) {
  
  
  
  
  if ( !file.exists(paste0(levelqc1.folder, "/all_split_reads_qc_level1.rds")) ) {
    
    logger::log_info("Error! The file with the split reads passing the LEVEL 1 filtering criteria does not exist!")
    break;
    
  } else {
    
    logger::log_info("loading the file with the split reads passing the LEVEL 1 filtering criteria...")
    all_split_reads_qc_level1 <- readRDS(file = paste0(levelqc1.folder, "/all_split_reads_qc_level1.rds"))
    
  }
  
  logger::log_info(all_split_reads_qc_level1 %>% nrow(), " initial number of split reads...")
  
  doParallel::registerDoParallel(num.cores)
  foreach(i = seq(length(recount3.project.IDs))) %dopar%{
    
    # i <- 1
    project_id <- recount3.project.IDs[i]
    
    
    
    ## If this is GTEx, the recount3.project.IDs object can contain up to 54 different values
    ## e.g. KIDNEY, BRAIN, BLOOD, etc
    # project_id <- recount3.project.IDs[1]
    # project_id <- "LUNG"
    # project_id <- "OV"
    
    local_folder_results <- paste0(results.folder, "/", project_id, "/base_data/")
    dir.create(file.path(local_folder_results), recursive = TRUE, showWarnings = F)
    
    

    logger::log_info(project_id, " - downloading junction data from recount3")
    
    
    
    rse_jxn <- recount3::create_rse_manual(
      project = project_id,
      project_home = data.source,
      organism = "human",
      annotation = "gencode_v29",
      type = "jxn"
    )
   
    metadata.info <- colData(rse_jxn)
    saveRDS(object = metadata.info, 
            file = paste0(local_folder_results, "/", project_id, "_samples_raw_metadata.rds"))
    
    
    #################################
    ## GET SAMPLE CLUSTERS
    #################################
    
    # metadata.info <- readRDS(file = paste0(local_folder_results, "/", project_id, "_samples_raw_metadata.rds"))
    metadata_tidy <- separate_clusters(recount3.project.IDs,
                                       project.metadata = metadata.info, 
                                       data.source)

    if ( metadata_tidy %>% nrow() > 0 ) {
      
      if (subsampling) {
        
        set.seed(12)
        
        if ( !all(is.na(metadata_tidy$all_mapped_reads)) ) {
          ## From the samples obtained, only keep the ones matching similar RIN numbers
          m.out <- MatchIt::matchit(cluster ~ all_mapped_reads,
                                    data = metadata_tidy,
                                    distance = metadata_tidy$all_mapped_reads,
                                    method = "nearest",
                                    caliper = c(all_mapped_reads = 500000),
                                    std.caliper = F)
        } else {
          ## From the samples obtained, only keep the ones matching similar age numbers
          m.out <- MatchIt::matchit(cluster ~ age,
                                    data = metadata_tidy,
                                    distance = metadata_tidy$age,
                                    method = "nearest",
                                    caliper = c(age = 1),
                                    std.caliper = F)
        }
        
        metadata_tidy_filter <- MatchIt::match.data(m.out) %>%
          dplyr::select(-c("distance", "weights", "subclass"))
        metadata_tidy_filter %>% 
          distinct(external_id, .keep_all = T) %>%
          dplyr::count(cluster)
        metadata_tidy_filter %>% 
          distinct(external_id, .keep_all = T) %>%
          as_tibble()
        
      } else {
        
        metadata_tidy_filter <- metadata_tidy
      }
      
      if ( !all(is.na(metadata_tidy$rin)) ) {
        stopifnot(
          "There are samples with RIN lower than 6!" =
            all(metadata_tidy_filter$rin >= 6)
        )
      }
      
      saveRDS(object = metadata_tidy_filter, 
              file = paste0(local_folder_results, "/", project_id, "_samples_metadata.rds"))
      
      # metadata_tidy_filter$all_mapped_reads %>% min
      # metadata_tidy_filter %>% as_tibble()
      # metadata_tidy_filter$external_id %>% sort()
      
      
      #############################################################
      ## FILTERS NEEDED TO CALL A SPLICE SITE 
      #############################################################
      
      project_samples <- metadata_tidy_filter %>% 
        distinct(sample_id) %>% 
        pull()
      
      logger::log_info("removing split reads from samples not passing the filtering criteria ...")
      
      ## Only samples passing the filtering criteria
      all_counts <- assay(rse_jxn, "counts")[,(colnames(assay(rse_jxn, "counts")) %in% project_samples), drop=FALSE]
      all_counts %>% nrow()
      
      rm(rse_jxn)
     
      logger::log_info("removing split reads that did not pass the LEVEL1 filtering criteria ...")
      
      ## Only consider novel and annotated junctions passing the LEVEL 1
      all_counts <- all_counts[drop=FALSE, (rownames(all_counts) %in% all_split_reads_qc_level1$junID),]
      all_counts %>% nrow()
      
      logger::log_info("removing split reads that do not present at least 1 split read ...")
      
      ##  Only consider novel and annotated junctions with at least N number of supporting split reads accross clusters
      all_counts <- all_counts[drop=FALSE, rowSums(all_counts) >= supporting.reads, ]
      all_counts %>% nrow()
      
      
      ## Convert to matrix for easy data manament across clusters
      all_counts <- all_counts %>% as.matrix()
      
      
      gc()
      
      ###################################
      ## GET SPLIT READ DATA PER CLUSTER
      ###################################
      
      ## Separate data per cluster
      clusters_ID <- metadata_tidy$cluster %>% unique() %>% as.character()
      logger::log_info(paste(clusters_ID, collapse = " & "))
      
      clusters_used <- NULL
      
      for ( cluster_id in clusters_ID ) {
        
        # cluster_id <- clusters_ID[1]
        # cluster_id <- clusters_ID[2]

        ################
        ## Get clusters
        ################
        
        ## Only consider samples that passed the filtering process performed above
        ## (e.g. RIN score and other filters, in case of the 'gtex' project) 
        cluster_samples <- metadata_tidy_filter %>% 
          filter(cluster == cluster_id) %>%
          distinct(sample_id) %>% 
          pull()
        
        logger::log_info(cluster_id, " samples: ", cluster_samples %>% length())
        
        ## When working with 'gtex' data, only tissues with at least 70 samples were considered
        if ( ((data.source == "data_sources/gtex" && cluster_samples %>% length() >= 70) ||
              (data.source != "data_sources/gtex" && cluster_samples %>% length() >= 1)) && 
             ((colnames(all_counts) %in% cluster_samples) %>% length() > 0) ) {
          
          saveRDS(object = cluster_samples, 
                  file = paste0(local_folder_results, "/", project_id, "_", cluster_id, "_samples_used.rds"))
          
          clusters_used <- c(clusters_used, cluster_id)
          
          ###################################
          ## Get counts for current cluster
          ###################################

          # 'drop=FALSE' to maintain matrix structure
          local_counts <- all_counts[,(colnames(all_counts) %in% cluster_samples), drop=FALSE]
          
          ##  Only consider novel and annotated junctions with at least N supporting split read across the samples of the current project
          local_counts <- local_counts[drop=FALSE,rowSums(local_counts) >= supporting.reads, ]
          local_counts %>% nrow()
          local_counts %>% ncol()
          
          ## Convert to tibble
          local_counts <- local_counts %>% as_tibble(rownames = "junID")
          local_counts %>% nrow()
          
          logger::log_info(cluster_id, " --> ", local_counts %>% nrow(), " final split reads.")
          
          ## QC check and save data
          logger::log_info("saving data....")
          
          stopifnot(
            "Still there are split reads with a lower number of supporting reads indicated by parameter." =
              local_counts %>%
              mutate(sumCounts = rowSums(dplyr::select(., !contains("junID")))) %>%
              filter(sumCounts < 1) %>%
              nrow() == 0
          )
          logger::log_info(local_counts %>% nrow(), " unique split read counts.")
          
          #######################
          ## Get split reads ID
          #######################

          all_split_reads <- local_counts %>%
            dplyr::select(junID) %>%
            data.table::as.data.table() %>%
            inner_join(y = all_split_reads_qc_level1,
                       by = "junID")
          
          ## Separate from novel combos
          all_split_reads_combos <- all_split_reads %>%
            filter(type == "novel_combo")
          
          all_split_reads <- all_split_reads %>%
            filter(!(junID %in% all_split_reads_combos$junID))
          
          #######################
          ## QC and save results
          #######################
          
          if (any(all_split_reads$width < 25) |
              any(str_detect(string = all_split_reads$chr, pattern = "random")) |
              any(str_detect(string = str_to_lower(all_split_reads$chr), pattern = "u")) |
              any(!(all_split_reads$type %in% c("annotated", "novel_acceptor", "novel_donor")))) {
            logger::log_info("ERROR! The split reads do not meet the minimum LEVEL1 filtering criteria.")
            break;
          }
          
          ## Check how much memory this object uses
          #logger::log_info(object.size(all_split_reads %>% data.table::as.data.table()), units = "GB")
          
          logger::log_info(all_split_reads %>% nrow(), " unique split reads.")
          
          ## Save split reads objects
          saveRDS(object = all_split_reads %>% as_tibble(),
                  file = paste0(local_folder_results, "/", project_id, "_", cluster_id, "_all_split_reads.rds"))
          saveRDS(object = all_split_reads_combos %>% as_tibble(),
                  file = paste0(local_folder_results, "/", project_id, "_", cluster_id, "_all_split_reads_combos.rds"))
          
          ## Save counts
          saveRDS(object = local_counts %>% filter(junID %in% all_split_reads$junID),
                  file = paste0(local_folder_results, "/", project_id, "_", 
                                cluster_id, "_split_read_counts.rds"))
          saveRDS(object = local_counts %>% filter(junID %in% all_split_reads_combos$junID),
                  file = paste0(local_folder_results, "/", project_id, "_", 
                                cluster_id, "_split_read_counts_combos.rds"))
          logger::log_info(local_counts %>% filter(junID %in% all_split_reads_combos$junID) %>% nrow(), " unique split read counts from combos.")
          gc()
          
          ## Free up some local memory
          rm(all_split_reads_combos)
          rm(all_split_reads)
          rm(local_counts)
          gc()
        }
        
        rm(cluster_samples)
        gc()
      }
      
      if ( !is.null(clusters_used) ) {
        saveRDS(object = clusters_used, 
                file = paste0(local_folder_results, "/", project_id, "_clusters_used.rds"))
      }
      
      rm(clusters_used)
      rm(all_counts)
      gc()
      
      
    } else {
      
      logger::log_info(paste0(project_id, " does not have any sample that meet the minimum criteria for inclusion."))
      
    }
    
    rm(metadata_tidy_filter)
    rm(metadata.info)
    rm(metadata_tidy)
    
    gc()
  }

}