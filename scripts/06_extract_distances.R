
#' Title
#'
#' @param cluster 
#' @param samples 
#' @param folder_name 
#'
#' @return
#' @export
#'
#' @examples
extract_distances <- function(cluster,
                              samples,
                              folder_name,
                              replace) {
  
  
  if ( replace ) {
    
    ## Obtain the distances across all samples
    df_all <- map_df(samples, function(sample) { 
      
      # sample <- samples[1]
      # sample <- "GTEX-ZVT2-0426-SM-5E44S.1"
      
      file_name <- paste0(folder_name, "/", cluster, "_", sample, "_distances.rds")
      
      if (file.exists(file_name)) {
        
        
        df <- readRDS(file = file_name)
        
        logger::log_info(paste0("sample: ", sample, " - ", df %>% nrow, " split read pairings!"))
        
        
        return(df)
        
      } else {
        
        logger::log_info(paste0("sample: ", sample, " doesn't exist!"))
        break;
      }
      
    })
    
    saveRDS(object = df_all %>%
              distinct(novel_junID, ref_junID, .keep_all = T) %>%
              mutate(tissue = cluster),
            file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds") )
    
   } else {
     logger::log_info(paste0("Sample '", cluster, "_raw_distances_tidy.rds' exists!"))
    }
  
  ## RELEASE SOME MEMORY
  rm(df_all)
  gc()
  
}
