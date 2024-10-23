#' Title
#' Obtains all junction pairings across all projectsID
#' @param recount3.project.IDs List of recount3 projects to analyse
#' @param all.clusters Clusters of samples. In GTEx projects, samples were clustered by tissue (eg. 'Puituitary', 'Thyroid', etc)
#' @param database.folder Path to the local folder that stores the database to be produced and the files needed to produce it
#' @param results.folder Local path to the folder that contains the results of the analyses performed
#'
#' @return
#' @export
#'
#' @examples
GetAllRawJxnPairings <- function(recount3.project.IDs,
                                 all.clusters = NULL,
                                 database.folder,
                                 num.cores,
                                 results.folder) {
  
  
  doParallel::registerDoParallel(num.cores)
  
  df_all_distances_pairings_raw <- foreach(i = seq(length(recount3.project.IDs)), .combine = "rbind") %dopar%{
    project_id <- recount3.project.IDs[i]
    
    # project_id <- recount3.project.IDs[1]
    # project_id <- "KIDNEY"
    
    logger::log_info(paste0(Sys.time(), " --> Working with '", project_id, "' DataBase..."))
    folder_results_root <- paste0(results.folder, "/", project_id, "/")
    
    if (is.null(all.clusters) && file.exists(paste0(folder_results_root, "/base_data/", project_id, "_clusters_used.rds"))) {
      all.clusters <-  readRDS(file = paste0(folder_results_root, "/base_data/", project_id, "_clusters_used.rds"))
    }
    
    if (!is.null(all.clusters)) {
      
      map_df(all.clusters, function(cluster) {
        
        # cluster <- all.clusters[1]
        
        logger::log_info(paste0(Sys.time(), " - ", project_id, " loading '", cluster, "'  data ..."))
        
        ## Load samples
        if (file.exists(paste0(folder_results_root, "/base_data/", project_id, "_", cluster, "_samples_used.rds"))) {
          
          samples <- readRDS(file = paste0(folder_results_root, "/base_data/", project_id, "_", cluster,  "_samples_used.rds"))
          
          if (samples %>% length() > 0) {
            
            folder_cluster_pairings <- paste0(folder_results_root, "/junction_pairing/", cluster, "/")
            
            if (!file.exists(paste0(folder_cluster_pairings, "/", cluster, "_raw_distances_tidy.rds"))) {
              
              ## Obtain the distances across all samples
              df_all <- map_df(samples, function(sample) { 
                # sample <- samples[1]
                file_name <- paste0(folder_cluster_pairings, "/", cluster, "_", sample, "_distances.rds")
                
                if (file.exists(file_name)) {
                  logger::log_info(paste0(cluster, " - ", sample))
                  readRDS(file = file_name) %>% return()
                } else {
                  return(NULL)
                }
              })
              
              if (nrow(df_all) > 0) {
                saveRDS(object = df_all %>% distinct(novel_junID, ref_junID, .keep_all = T) %>% mutate(tissue = cluster),
                        file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds"))
              }
              
            } else {
              logger::log_info(paste0("File '", cluster, "_raw_distances_tidy.rds' already exists!"))
              df_all <- readRDS( file = paste0(folder_cluster_pairings, "/", cluster, "_raw_distances_tidy.rds") )
            }
            
            
            if (nrow(df_all) > 0) {
              
              df_all %>% distinct(novel_junID, ref_junID, .keep_all = T) %>% mutate(project = project_id) %>%
                return()
              
            } else {
              return(NULL)
            }          
            
          } else {
            logger::log_info(paste0("ERROR: no samples available for the tissue: ", project_id))
            return(NULL)
          }
          
        } else {
          return(NULL)
        }
      })  
    }
  }
  
  logger::log_info("\t\t Saving 'df_all_distances_pairings' file for the database!")
  
  
  saveRDS(object = df_all_distances_pairings_raw %>% distinct(project) %>% pull(),
          file = file.path(results.folder, "/all_final_projects_used.rds"))
  
  saveRDS(object = df_all_distances_pairings_raw %>% distinct(novel_junID, ref_junID, .keep_all = T),
          file = file.path(database.folder,"/all_raw_jxn_pairings.rds"))
  
}


GetAllRawNovelCombos <- function(recount3.project.IDs,
                                 database.folder,
                                 results.folder) {
  
  ## Get all novel combos across all tables
  all_split_reads_combos <- map_df(recount3.project.IDs, function(project_id) { 
    
    # project_id <- recount3.project.IDs[1]
    
    logger::log_info("Working with '", project_id, "' ...")
    results_folder_local <- paste0(results.folder, "/", project_id, "/")
    
    clusters <- readRDS(file = paste0(results_folder_local, "/base_data/", project_id, "_clusters_used.rds"))
    
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
  
  all_split_read_combos_RBPs <- all_split_reads_combos %>% distinct(junID, .keep_all = T) %>% dplyr::rename(ref_junID = junID)
  
  saveRDS(object = all_split_read_combos_RBPs, file = file.path(database.folder, "all_novel_combos.rds"))
}
