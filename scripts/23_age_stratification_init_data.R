#' Title
#' Clusters the samples by age supergroup
#' @param projects.id 
#' @param gtf.version 
#' @param project.name 
#' @param results.folder 
#'
#' @return
#' @export
#'
#' @examples
age_stratification_init_data <- function (projects.id,
                                          gtf.version,
                                          project.name,
                                          results.folder) {
  
  # project_id <- projects.id[1]
  map_df(projects.id, function(project_id) {
    
    logger::log_info(" --> ", project_id, " - getting sample's metadata.")
    
    ## Samples are loaded from the splicing project
    sample_metadata <- readRDS(file = file.path(results.folder, project_id, "/base_data/", 
                                                paste0(project_id, "_samples_metadata.rds")) )
    ## Cluster the samples by age
    age_samples_clusters_tidy <- age_stratification_set_metadata(sample.metadata = sample_metadata, project.id = project_id) %>%
      mutate(age_group = "")
    
    for ( age_cluster in (age_samples_clusters_tidy$age %>% unique()) ) {
      
      # age_cluster = (age_samples_clusters_tidy$age %>% unique())[1]

      
      switch(age_cluster, 
             '20-29' = {
               indx <- which(age_samples_clusters_tidy$age == age_cluster)
               age_samples_clusters_tidy[indx, "age_group"] <- "20-39"
             },
             '30-39' = {
               indx <- which(age_samples_clusters_tidy$age == age_cluster)
               age_samples_clusters_tidy[indx, "age_group"] <- "20-39"
             },
             '40-49' = {
               indx <- which(age_samples_clusters_tidy$age == age_cluster)
               age_samples_clusters_tidy[indx, "age_group"] <- "40-59"
             },
             '50-59' = {
               indx <- which(age_samples_clusters_tidy$age == age_cluster)
               age_samples_clusters_tidy[indx, "age_group"] <- "40-59"
             },
             '60-69' = {
               indx <- which(age_samples_clusters_tidy$age == age_cluster)
               age_samples_clusters_tidy[indx, "age_group"] <- "60-79"
             },
             '70-79' = {
               indx <- which(age_samples_clusters_tidy$age == age_cluster)
               age_samples_clusters_tidy[indx, "age_group"] <- "60-79"
             },
             {
               indx <- which(age_samples_clusters_tidy$age == age_cluster)
               age_samples_clusters_tidy[indx, "age_group"] <- NA
               logger::log_info('age cluster ', age_cluster, 'not found')
             }
      )
    }
    
    age_samples_clusters_tidy <- age_samples_clusters_tidy %>%
      mutate(rin_transformed = rin / 10) %>%
      mutate(age_group_trans = (age_group %>% as.factor %>% as.integer())/3)
    
    ## If the tissue has the three age groups with at least 25 samples
    if ( (age_samples_clusters_tidy %>%
          dplyr::count(age_group) %>%
          drop_na() %>%
          filter(n >= 25) %>%
          pull(n)) %>% length() == 3 ) {
      
      
      # "40-59" ----------------------------------------------
      
      data_combined <- age_samples_clusters_tidy %>%
        filter(age_group %in% c("20-39","40-59")) 
      
      m.out0 <- MatchIt::matchit(age_group_trans ~ rin_transformed,
                                 data = data_combined,
                                 method = "optimal",
                                 ratio = 1, replace = F)
      
      st_match <- MatchIt::match.data(m.out0) %>%
        dplyr::select(-c(distance,weights, subclass))
      
      # m.out <- MatchIt::matchit(formula = age_group ~ rin_transformed, 
      #                           data = age_samples_clusters_tidy, 
      #                           distance = age_samples_clusters_tidy$rin_transformed,
      #                           method = "nearest", 
      #                           caliper = c(rin_transformed = 0.005), 
      #                           std.caliper = FALSE, 
      #                           ratio = 1)
      # 
      # st_match <- MatchIt::match.data(m.out) %>%
      #   dplyr::select(-c(distance,weights, subclass))
      
      
      
      
      # "60-79" ----------------------------------------------
      
      df_youngest <- st_match %>%
        filter(age_group == "20-39") 
      
      df_eldest <- age_samples_clusters_tidy %>%
        filter(age_group == "60-79")
      
      data_combined <- rbind(df_youngest, df_eldest)
      
      
      m.out <- MatchIt::matchit(age_group_trans ~ rin_transformed,
                                data = data_combined,
                                method = "optimal",
                                ratio = 1, 
                                replace = F)
      
      nd_match <- MatchIt::match.data(m.out) %>%
        dplyr::select(-c( distance, weights, subclass))
      
      
      
      # Join both ----------------------------------------------
      
      
      age_samples_clusters_subsampled <- rbind(st_match, nd_match) %>%
        distinct(individual, .keep_all = T) 
      
      logger::log_info(project_id, " --> ", paste(age_samples_clusters_subsampled %>% dplyr::count(age_group),collapse = " - "))
      
      
      
      
      return(age_samples_clusters_subsampled)
      
    } else {
      return(NULL)
    }
  })
  
}
