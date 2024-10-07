
#' Title
#' Function to pair novel junctions with annotated introns across the samples of each sample cluster 
#' (i.e. case/control, tissue, etc)
#' @param recount3.project.IDs array with the recount3 projects to download
#' e.g. "SRP009615" 
#'
#' @return
#' @export
#'
#' @examples
JunctionPairing <- function(recount3.project.IDs, 
                            results.folder,
                            replace,
                            num.cores) {
  
  
  doParallel::registerDoParallel(num.cores)
  foreach(i = seq(length(recount3.project.IDs)), .combine = "rbind") %dopar%{
    
    # i<- 1
    project_id <- recount3.project.IDs[i]
    
    # project_id <- recount3.project.IDs[3]
    # project_id <- recount3.project.IDs[2]
    
    folder_root <- file.path(results.folder, project_id)
    folder_base_data <- file.path(folder_root, "base_data")
    
    logger::log_info("\t\t Getting data from '", project_id, "' project...")
    
    ## Load clusters
    
    if (file.exists(paste0(folder_base_data, "/", project_id, "_clusters_used.rds"))) {
      
      clusters_ID <- readRDS(file = paste0(folder_base_data, "/", project_id, "_clusters_used.rds"))
      
      for (cluster_id in clusters_ID) {
        
        # cluster_id <- clusters_ID[1]
        
        logger::log_info(paste0(Sys.time(), " - loading '", cluster_id, "' source data ..."))
        
        ############################################
        ## LOAD DATA FOR THE CURRENT PROJECT
        ############################################
        
        if ( file.exists(paste0(folder_base_data, "/", project_id, "_", cluster_id, "_samples_used.rds")) && 
             file.exists(paste0(folder_base_data, "/", project_id, "_", cluster_id, "_all_split_reads.rds")) && 
             file.exists(paste0(folder_base_data, "/", project_id, "_", cluster_id, "_split_read_counts.rds")) ) {
          
          ## Load samples
          samples_used <- readRDS(file = paste0(folder_base_data, "/", project_id, "_", cluster_id, "_samples_used.rds")) %>% unique()
        
          
          ## Load split read data
          all_split_reads_details <- readRDS(file = paste0(folder_base_data, "/", project_id, "_", cluster_id, "_all_split_reads.rds")) %>% as_tibble()
          
          
          ## Load split read counts
          split_read_counts <- readRDS(file = paste0(folder_base_data, "/", project_id, "_", cluster_id, "_split_read_counts.rds")) %>% as_tibble()
  
          
          if (!identical((split_read_counts %>% names())[-1] %>% sort(), samples_used %>% sort())) {
            logger::log_info("ERROR! different number of samples used!")
            stop("ERROR! different number of samples used!");
          }
          
          if (!identical(all_split_reads_details$junID, split_read_counts$junID)) {
            logger::log_info("ERROR! The number of junctions considered is not correct.")
            stop("ERROR! The number of junctions considered is not correct.");
          }
          
          ############################################
          ## DISTANCES SUITE OF FUNCTIONS
          ############################################
  
          folder_pairing_results <- file.path(folder_root, "junction_pairing", cluster_id)
          dir.create(file.path(folder_pairing_results), recursive = TRUE, showWarnings = T)
  
          GetDistances(project.id = project_id,
                       cluster = cluster_id,
                       samples = samples_used,
                       split.read.counts = split_read_counts,
                       all.split.reads.details = all_split_reads_details,
                       folder.name = folder_pairing_results,
                       replace = replace)
          gc()
  
  
          ExtractDistances(cluster = cluster_id,
                           samples = samples_used,
                           folder.name = folder_pairing_results,
                           replace = replace)
          gc()
  
  
          GetNeverMisspliced(cluster = cluster_id,
                             samples = samples_used,
                             split.read.counts = split_read_counts,
                             all.split.reads.details = all_split_reads_details,
                             folder.name = folder_pairing_results,
                             replace = replace)
  
          
          rm(all_split_reads_details)
          rm(split_read_counts)
          gc()
          
        }
      }
    }
  }
}
