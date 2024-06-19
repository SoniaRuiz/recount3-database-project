#' Title
#' Cl
#' @param sample.metadata 
#' @param project.id 
#'
#' @return
#' @export
#'
#' @examples
age_stratification_set_metadata <- function (sample.metadata, 
                                             project.id) {
  
  logger::log_info(" --> ", project.id, " - performing the sample clustering.")
  
  ## 1. Get key metadata
  
  clusters <- sample.metadata$cluster %>% unique()
  age_groups <- sample.metadata$age %>% unique()
  samples <- sample.metadata$sample_id %>% unique()
  
  ## 2. Loop through the individuals
  
  gtex_age_clustering <- map_df(samples, function(individual) {
    
    # individual <- (samples)[1]
    
    person <- sample.metadata %>%
      filter(sample_id == individual)
    
    # print(person$smtsd %>% unique)
    
    if ( person$age %>% unique() %in% age_groups ) {
      
      if ( any(person$cluster %in% clusters) ) {
        
        sex <- ifelse(person$gender == 1, "male", "female")
        
        data.frame(individual = person$sample_id %>% unique(),
                   age = person$age %>% unique(),
                   rin = person$rin,
                   mapped_read_count = person$all_mapped_reads,
                   sex = sex,
                   #sample_recount_id = person$gtex.sampid,
                   cluster = person$cluster,
                   project = project.id) %>% return()
      }
      
    }
    
  }) 
  
  gtex_age_clustering %>% 
    return()
  
}