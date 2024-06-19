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
get_all_raw_jxn_pairings <- function(recount3.project.IDs,
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
    
    if ( is.null(all.clusters) && 
         file.exists(paste0(folder_results_root, "/base_data/", 
                            project_id, "_clusters_used.rds")) ) {
      all.clusters <-  readRDS(file = paste0(folder_results_root, "/base_data/", 
                                             project_id, "_clusters_used.rds"))
    }
    
    if ( ! is.null(all.clusters) ) {
      
      map_df(all.clusters, function(cluster) {
        
        # cluster <- all.clusters[1]
        
        logger::log_info(paste0(Sys.time(), " - ", project_id, " loading '", cluster, "'  data ..."))
        
        ## Load samples
        if ( file.exists(paste0(folder_results_root, "/base_data/", project_id, "_", cluster, "_samples_used.rds")) ) {
          
          samples <- readRDS(file = paste0(folder_results_root, "/base_data/", project_id, "_", cluster,  "_samples_used.rds"))
          
          if ( samples %>% length() > 0 ) {
            
            folder_cluster_pairings <- paste0(folder_results_root, "/junction_pairing/", cluster, "/")
            
            if ( !file.exists(paste0(folder_cluster_pairings, "/", cluster, "_raw_distances_tidy.rds")) ) {
              
              ## Obtain the distances across all samples
              df_all <- map_df(samples, function(sample) { 
                
                # sample <- samples[1]
                
                file_name <- paste0(folder_cluster_pairings, "/", cluster, "_", sample, "_distances.rds")
                
                
                if ( file.exists(file_name) ) {
                  logger::log_info(paste0(cluster, " - ", sample))
                  df <- readRDS(file = file_name)
                  
                  return(df)
                } else {
                  return(NULL)
                }
                
              })
              
              if ( nrow(df_all) > 0 ) {
                saveRDS(object = df_all %>%
                          distinct(novel_junID, ref_junID, .keep_all = T) %>%
                          mutate(tissue = cluster),
                        file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds"))
              }
              
            } else {
              logger::log_info(paste0("File '", cluster, "_raw_distances_tidy.rds' already exists!"))
              df_all <- readRDS( file = paste0(folder_cluster_pairings, "/", cluster, "_raw_distances_tidy.rds") )
            }
            
            
            if ( nrow(df_all) > 0 ) {
              
              df_all %>%
                distinct(novel_junID, ref_junID, .keep_all = T) %>%
                mutate(project = project_id) %>%
                return()
              
            } else {
              return(NULL)
            }          
            # df_all2 <- readRDS(file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds"))
            
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
  
  logger::log_info(paste0(Sys.time(), " - saving 'df_all_distances_pairings' for the database!"))
  
  
  saveRDS(object = df_all_distances_pairings_raw %>% 
            distinct(project) %>% pull(),
          file = paste0(results.folder, "/all_final_projects_used.rds"))
  
  saveRDS(object = df_all_distances_pairings_raw %>%
            distinct(novel_junID, ref_junID, .keep_all = T),
          file = paste0(database.folder,"/all_raw_jxn_pairings.rds"))
  
}