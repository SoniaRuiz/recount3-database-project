#' Title
#' Load all annotated introns with no evidence of mis-splicing activity across all recount3 projects indicated
#' @param recount3.project.IDs List of recount3 projects to analyse
#' @param all.clusters Clusters of samples. In GTEx projects, samples were clustered by tissue (eg. 'Puituitary', 'Thyroid', etc)
#' @param database.folder Path to the local folder that stores the database to be produced and the files needed to produce it
#' @param results.folder Local path to the folder that contains the results of the analyses performed
#'
#' @return
#' @export
#'
#' @examples
get_all_intron_never_misspliced <- function (recount3.project.IDs,
                                             all.clusters = NULL,
                                             database.folder,
                                             results.folder) {
  
  
  df_never <- map_df(recount3.project.IDs, function(project_id) {
    
    # project_id <- recount3.project.IDs[10]
    
    # logger::log_info(paste0(Sys.time(), " --> Working with '", project_id, "' DataBase..."))
    base_folder <- paste0(results.folder, "/", project_id)
    
    if ( is.null(all.clusters) && 
         file.exists(paste0(base_folder, "/base_data/", project_id, "_clusters_used.rds"))) {
      
      all.clusters <-  readRDS(file = paste0(base_folder, "/base_data/", project_id, "_clusters_used.rds"))
      
    } 
    
    if ( !is.null(all.clusters) ) {
      
      map_df(all.clusters, function(cluster) { 
        
        # cluster <- all.clusters[1]
        
        if ( file.exists(paste0(base_folder, "/junction_pairing/", cluster, 
                                "/not-misspliced/", cluster, "_all_notmisspliced.rds")) ) {
          # logger::log_info(paste0(Sys.time(), " --> ", cluster))
          df_introns_never <- readRDS(file = paste0(base_folder, "/junction_pairing/", cluster, 
                                                    "/not-misspliced/",cluster, "_all_notmisspliced.rds")) %>% as_tibble()
          return(data.frame(ref_junID = df_introns_never$value))
        } else {
          return(NULL)
        }
        
      })
    }
  })
  
  df_never %>% distinct(ref_junID) %>% return()
}
