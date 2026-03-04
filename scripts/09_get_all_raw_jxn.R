#' Title
#' Loops through the projects from the current recount3 project to obtain all original unique split reads found across their samples
#' @param recount3.project.IDs List of recount3 projects to analyse
#' @param all.clusters Clusters of samples. In GTEx projects, samples were clustered by tissue (eg. 'Puituitary', 'Thyroid', etc)
#' @param database.folder Local path to the folder that will contain the database and results files related with the database
#' @param results.folder Local path to the folder that contains the results of the analyses performed
#'
#' @return
#' @export
#'
#' @examples
GetAllAnnotatedSplitReads <- function(recount3.project.IDs,
                                      database.folder,
                                      results.folder,
                                      num.cores,
                                      replace,
                                      all.clusters = NULL) {
  
  
  
  #############################################
  ## These are all split reads from all tissues
  ## obtained from 'USE ME' SAMPLES and samples with more than 6 RIN
  #############################################
  
  if (replace || !file.exists(file.path(database.folder, "all_split_reads_qc_level2.rds"))) {
    
    doParallel::registerDoParallel(num.cores)
    all_split_reads_details_all_sample_clusters <- foreach(i = seq(length(recount3.project.IDs)), .combine = "rbind") %dopar%{
                                                             
      project_id <- recount3.project.IDs[i]
      
      # project_id <- recount3.project.IDs[1]
      # project_id <- recount3.project.IDs[2]
      
      folder_results_root <- file.path(results.folder, project_id)
            
      if (is.null(all.clusters) && 
          file.exists(file.path(folder_results_root, "base_data", paste0(project_id, "_clusters_used.rds")))) {        
        all.clusters <-  readRDS(file = file.path(folder_results_root, "base_data", paste0(project_id, "_clusters_used.rds")))
      } 
      
      message(project_id, " ", all.clusters)
      
      all_jxn_qc <- if (!is.null(all.clusters)) {
        
        all_jxn_qc <- map_df(all.clusters, function(cluster) {
          
          # cluster <- all.clusters[1]
          logger::log_info(project_id, " loading '", cluster, "'  data ...")
          
          if (file.exists(file.path(folder_results_root, "base_data", paste0(project_id, "_", cluster, "_all_split_reads.rds")))) {
            
            all_split_reads_details <- readRDS(file = file.path(folder_results_root, "base_data", paste0(project_id, "_", cluster, "_all_split_reads.rds")))
            all_split_reads_details |>
              distinct(junID, .keep_all = T) |>
              as_tibble()
            
          } 
          
        })
        
        if (all_jxn_qc %>% nrow() > 0 ) {
          all_jxn_qc |> distinct(junID, .keep_all = T)
        } 
      }
    }
        
    logger::log_info("Saving 'all_split_reads_qc_level2.rds' for the database!")
        
    ## SAVE DATA
    dir.create(path = database.folder, recursive = T, showWarnings = F)

    # Ensure junID doesn't have the 'chr' string 
    all_split_reads_details_level2 <- all_split_reads_details_all_sample_clusters |>
      distinct(junID, .keep_all = T)
    
    #all_split_reads_details_level2 <- all_split_reads_details_all_sample_clusters %>% distinct(junID, .keep_all = T) %>% GRanges()
    #seqlevelsStyle(all_split_reads_details_level2) <- "Ensembl"
    #all_split_reads_details_level2$junID <- as.character(all_split_reads_details_level2)

    saveRDS(object = all_split_reads_details_level2 |> 
              as_tibble() |>
              distinct(junID, .keep_all = T),
            file = file.path(database.folder, "all_split_reads_qc_level2.rds") )
  }
}

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
                                 database.folder,
                                 results.folder,
                                 num.cores,
                                 replace,
                                 all.clusters = NULL) {
  
  if (replace || !file.exists(file.path(database.folder,"all_raw_jxn_pairings.rds"))) {
    
    doParallel::registerDoParallel(num.cores)
    df_all_distances_pairings_raw <- foreach(i = seq(length(recount3.project.IDs)), .combine = "rbind") %dopar%{
      
      project_id <- recount3.project.IDs[i]
      
      # project_id <- recount3.project.IDs[3]
      
      folder_results_root <- file.path(results.folder, project_id)
      folder_results_base <- file.path(folder_results_root, "base_data")
      folder_results_jxn <- file.path(folder_results_root, "junction_pairing")
      
      if (is.null(all.clusters) && file.exists(file.path(folder_results_base, paste0(project_id, "_clusters_used.rds")))) {
        all.clusters <-  readRDS(file = file.path(folder_results_base, paste0(project_id, "_clusters_used.rds")))
      }
      
      if (!is.null(all.clusters)) {
        
        map_df(all.clusters, function(cluster) {
          
          # cluster <- all.clusters[2]
                    
          ## Load samples
          if (file.exists(file.path(folder_results_base, paste0(project_id, "_", cluster, "_samples_used.rds")))) {
            
            samples <- readRDS(file = file.path(folder_results_base, paste0(project_id, "_", cluster,  "_samples_used.rds")))
            
            if (samples %>% length() > 0) {
              
              folder_cluster_pairings <- file.path(folder_results_jxn, cluster)
              
              df_all <- if (!file.exists(file.path(folder_cluster_pairings, paste0(cluster, "_raw_distances_tidy.rds"))) || replace) {
                
                ## Obtain the distances across all samples
                df_all <- map_df(samples, function(sample) { 
                  # sample <- samples[1]
                  file_name <- file.path(folder_cluster_pairings, paste0(cluster, "_", sample, "_distances.rds"))
                  
                  if (file.exists(file_name)) {
                    
                    logger::log_info(project_id, " - ", cluster, " - ", sample)
                    readRDS(file = file_name)
                  } 
                })
                
                if (nrow(df_all) > 0) {
                  saveRDS(object = df_all %>% distinct(novel_junID, ref_junID, .keep_all = T) %>% mutate("cluster" = cluster),
                          file = file.path(folder_cluster_pairings, paste0(cluster, "_raw_distances_tidy.rds")))
                }
                
                df_all
                
              } else {
                logger::log_info("File '", cluster, "_raw_distances_tidy.rds' already exists!")
                readRDS(file = file.path(folder_cluster_pairings, paste0(cluster, "_raw_distances_tidy.rds")))
              }
              
              
              if (nrow(df_all) > 0) {
                print(nrow(df_all))
                df_all %>% 
                  distinct(novel_junID, ref_junID, .keep_all = T) %>% 
                  mutate(project = project_id)
                
              }           
              
            } else {
              stop("ERROR: no samples available for the project: '", project_id, "'")
            }
            
          }
        })  
      }
    }
    
    logger::log_info("Saving 'all_raw_jxn_pairings.rds' file for the database...")
    print(head(df_all_distances_pairings_raw))    
    saveRDS(object = df_all_distances_pairings_raw |>
              distinct(project) |>
              pull(),
            file = file.path(results.folder, "all_final_projects_used.rds"))
    
    saveRDS(object = df_all_distances_pairings_raw %>% distinct(novel_junID, ref_junID, .keep_all = T),
            file = file.path(database.folder,"all_raw_jxn_pairings.rds"))
  } else {
    message("File 'all_raw_jxn_pairings.rds' already exists!")
  }
  
}


GetAllRawNovelCombos <- function(recount3.project.IDs,
                                 database.folder,
                                 results.folder,
                                 replace) {
  
  if (replace || !file.exists(file.path(database.folder, "all_raw_novel_combos.rds"))) {
    
    logger::log_info("Function 'GetAllRawNovelCombos' ...")
    
    ## Get all novel combos across all tables
    all_split_reads_combos <- map_df(recount3.project.IDs, function(project_id) { 
      
      # project_id <- recount3.project.IDs[6]
      
      # logger::log_info("Working with '", project_id, "' ...")
      results_folder_local <- file.path(results.folder, project_id, "base_data")
      
      if (file.exists(file.path(results_folder_local, paste0(project_id, "_clusters_used.rds")))) {
        clusters <- readRDS(file = file.path(results_folder_local, paste0(project_id, "_clusters_used.rds")))  
    
        map_df(clusters, function(cluster_id) { 
          
          # cluster_id <- clusters[1]
          # logger::log_info(project_id, " --> ", cluster_id)
          
          if (file.exists(file.path(results_folder_local, paste0(project_id, "_", cluster_id, "_all_split_reads_combos.rds")))) {
            
            ## Load all split reads
            readRDS(file = file.path(results_folder_local, paste0(project_id, "_", cluster_id, "_all_split_reads_combos.rds"))) %>%
              dplyr::select(-any_of(c("in_ref", "n_projects", "annotated"))) %>%
              unnest(gene_id)
          }
        })
      }
    })

    logger::log_info("Saving 'all_raw_novel_combos.rds' file for the database...")
    
    saveRDS(object = all_split_reads_combos %>% distinct(junID, .keep_all = T),
            file = file.path(database.folder, "all_raw_novel_combos.rds"))
  } else {
    message("File 'all_raw_novel_combos.rds' already exists!")
  }
  
}


GetAllRawAmbiguous <- function(recount3.project.IDs,
                               database.folder,
                               results.folder,
                               replace) {
  
  if (replace || !file.exists(file.path(database.folder, "all_raw_ambiguous.rds"))) {
    
    logger::log_info("Function 'GetAllRawAmbiguous' ...")
    
    ## Get all ambiguous junctions across all tables
    all_split_reads_ambig <- map_df(recount3.project.IDs, function(project_id) { 
      
      # project_id <- recount3.project.IDs[1]
 
      results_folder_local <- file.path(results.folder, project_id, "base_data")

      if (file.exists(file.path(results_folder_local, paste0(project_id, "_clusters_used.rds")))) {
        clusters <- readRDS(file = file.path(results_folder_local, paste0(project_id, "_clusters_used.rds")))
        
        map_df(clusters, function(cluster_id) { 
          
          # cluster_id <- clusters[1]
          # logger::log_info(project_id, " --> ", cluster_id)
          
          if (file.exists(file.path(results_folder_local, paste0(project_id, "_", cluster_id, "_all_split_reads_ambig.rds")))) {
            
            ## Load all split reads
            readRDS(file = file.path(results_folder_local, paste0(project_id, "_", cluster_id, "_all_split_reads_ambig.rds"))) %>%
              dplyr::select(-any_of(c("in_ref","n_projects","annotated"))) %>%
              unnest(gene_id) %>% 
              distinct(junID, .keep_all = T)
          }
        })
      }
    })
    
    if (nrow(all_split_reads_ambig) == 0) {
      stop("No ambiguous junctions found!")
    }
    logger::log_info("Saving 'all_raw_ambiguous.rds' file for the database...")
    
    saveRDS(object = all_split_reads_ambig %>% distinct(junID, .keep_all = T), 
            file = file.path(database.folder, "all_raw_ambiguous.rds"))
  } else {
    message("File 'all_raw_ambiguous.rds' already exists!")
  } 
  
}

GetAllRawUnannotated <- function(recount3.project.IDs,
                                 database.folder,
                                 results.folder,
                                 replace) { 
  
  if (replace || !file.exists(file.path(database.folder, "all_raw_novel_unannotated.rds"))) {
    
    logger::log_info("Function 'GetAllRawUnannotated' ...")
    
    ## Get all novel unannotated across all tables
    all_split_reads_unannot <- map_df(recount3.project.IDs, function(project_id) { 
      
      # project_id <- recount3.project.IDs[1]
      
      # logger::log_info("Working with '", project_id, "' ...")
      results_folder_local <- file.path(results.folder, project_id, "base_data")
      
      if (file.exists(file.path(results_folder_local, paste0(project_id, "_clusters_used.rds")))) {
        clusters <- readRDS(file = file.path(results_folder_local, paste0(project_id, "_clusters_used.rds")))
        
        map_df(clusters, function(cluster_id) { 
          
          # cluster_id <- clusters[1]
          # logger::log_info(project_id, " --> ", cluster_id)
          
          if (file.exists(file.path(results_folder_local, paste0(project_id, "_", cluster_id, "_all_split_reads_unannotated.rds")))) {
            
            ## Load all split reads
            readRDS(file = file.path(results_folder_local, paste0(project_id, "_", cluster_id, "_all_split_reads_unannotated.rds"))) %>%
              dplyr::select(-any_of(c("in_ref","n_projects","annotated")))
          }
        })
      }
    })
        
    logger::log_info("Saving 'all_raw_novel_unannotated.rds' file for the database...")
    
    saveRDS(object = all_split_reads_unannot %>% distinct(junID, .keep_all = T), 
            file = file.path(database.folder, "all_raw_novel_unannotated.rds"))
  } else {
    message("File 'all_raw_novel_unannotated.rds' already exists!")
  } 
}