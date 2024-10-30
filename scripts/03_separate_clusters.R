#' Title
#' Each recount3 project has a different data source, hence the metadata is structured differently
#' This function extracts the metadata (clusters, RIN, etc) from each recount3 project depending on its data source
#' @param project.metadata raw metadata object as provided by recount3 (https://rdrr.io/bioc/recount3/man/read_metadata.html)
#' @param data.source source of the project in recount3. It can be:
#' - "data_sources/sra"  
#' - "data_sources/gtex" 
#' - "data_sources/tcga"
#'
#' @return
#' an standardised metadata object 
#' @export
#'
#' @examples
SeparateClusters <- function(recount3.project.IDs,
                             project.metadata,
                             data.source) {
  
  
  if ( data.source == "data_sources/gtex" ) {
    
    project_metadata_tidy <- project.metadata %>%
      as_tibble()  %>%
      mutate(gtex.smrin = gtex.smrin %>% as.double()) %>%
      filter(gtex.smafrze != "EXCLUDE",
             #!(gtex.smtsd %in% c("Brain - Cortex", "Brain - Cerebellum")),
             gtex.smrin >= 6.0) ## Only fresh-frozen preserved tissues
    
    
    stopifnot(
      "Still there are samples with RIN < 6 and labelled as EXCLUDE" =
        !any(project_metadata_tidy$gtex.smrin < 6)
    )
    stopifnot(
      "Still there are samples labelled as EXCLUDE" =
        !any(project_metadata_tidy$gtex.smafrze == "EXCLUDE")
    )
    
    
    ## TODO same standard column names as with data_sources/sra projects
    project_metadata_tidy <- data.frame(sample_id =  project_metadata_tidy$external_id %>% as.character(),
                                        age = project_metadata_tidy$gtex.age %>% as.character(),
                                        rin = project_metadata_tidy$gtex.smrin %>% as.double(),
                                        gender = project_metadata_tidy$gtex.sex %>% as.character(),
                                        cluster = project_metadata_tidy$gtex.smtsd,
                                        smnabtcht = project_metadata_tidy$gtex.smnabtcht,
                                        smafrze = project_metadata_tidy$gtex.smafrze,
                                        avg_read_length = project_metadata_tidy$recount_seq_qc.avg_len,
                                        all_mapped_reads = project_metadata_tidy$recount_qc.star.all_mapped_reads,
                                        SRA_project = project_metadata_tidy$recount_project.project) 
    
    
    
    
  } else if (data.source == "data_sources/sra") {
    
    project_metadata_tidy <- project.metadata %>%
      as_tibble() %>%
      dplyr::select(external_id, 
                    sra.experiment_title, 
                    sra.sample_attributes, 
                    all_mapped_reads = recount_qc.star.all_mapped_reads) %>%
      mutate(rn = row_number()) %>%
      separate_rows(sra.sample_attributes, sep = "\\|\\s*") %>%
      separate(sra.sample_attributes, into = c('col1', 'col2'), sep = ";;") %>% 
      pivot_wider(names_from = col1, values_from = col2) %>% 
      dplyr::select(-rn) %>%
      as.data.frame() %>%
      dplyr::select(-any_of("sample_id")) %>%
      dplyr::rename(sample_id = external_id)
    
    if ( recount3.project.IDs == "SRP151040" ) {
      
      project_metadata_tidy <- project_metadata_tidy %>%
        rowwise() %>%
        dplyr::mutate(cluster = paste0(`cell type`, " - " ,`disease state`) %>% str_remove_all(pattern = "\\'")) %>%
        ungroup()
      
    } else if ( recount3.project.IDs == "SRP100948" ) {
      
      project_metadata_tidy <- project_metadata_tidy %>%
        rowwise() %>%
        dplyr::mutate(cluster = diagnosis %>% str_remove_all(pattern = "\\'")) %>%
        ungroup() %>%
        mutate(rin_score = if (exists('rin_score', where = project_metadata_tidy)) rin_score else NA) %>%
        mutate_at(vars(one_of('rin_score')), as.double) %>%
        dplyr::rename(age = "age at death",
                      rin = "rin_score",
                      gender = if (exists('Sex', where = project_metadata_tidy)) "Sex" else "gender") %>%
        mutate(age = age %>% as.integer()) %>%
        as_tibble()
      
    } else if ( recount3.project.IDs == "SRP181886" ) {
      
      project_metadata_tidy <- project_metadata_tidy %>%
        rowwise() %>%
        dplyr::mutate(cluster = diagnosis %>% str_remove_all(pattern = "\\'")) %>%
        ungroup() %>%
        mutate(rin_score = if (exists('rin_score', where = project_metadata_tidy)) rin_score else NA) %>%
        mutate_at(vars(one_of('rin_score')), as.double) %>%
        dplyr::rename(rin = "rin_score",
                      gender = if (exists('Sex', where = project_metadata_tidy)) "Sex" else "gender") %>%
        mutate(age = age %>% as.integer()) %>%
        as_tibble()
      
    } else if ( recount3.project.IDs == "SRP058181" ) {
      
      project_metadata_tidy <- project_metadata_tidy %>%
        mutate(cluster = ifelse( test = str_detect(sra.experiment_title, pattern="P"),
                                 yes = "Parkinsons Disease",
                                 no = "Control")) %>%
        mutate(rin_score = if (exists('rin_score', where = project_metadata_tidy)) rin_score else NA) %>%
        mutate_at(vars(one_of('rin_score')), as.double) %>%
        dplyr::rename(age = "age at death",
                      rin = "rin_score",
                      gender = if (exists('Sex', where = project_metadata_tidy)) "Sex" else "gender") %>%
        mutate(age = age %>% as.integer()) %>%
        as_tibble()
    }
    
    
    project_metadata_tidy <- project_metadata_tidy %>%
      mutate(cluster = cluster %>% as.factor(),
             #rin_score = rin_score %>% as.double(),
             avg_read_length = project.metadata %>% as_tibble %>% filter(external_id %in% project_metadata_tidy$sample_id) %>% pull(recount_seq_qc.avg_len),
             SRA_project = project.metadata %>% as_tibble %>% filter(external_id %in% project_metadata_tidy$sample_id) %>% pull(recount_project.project) )
    
    #if ( any(project_metadata_tidy %>% names == "rin") && !all(is.na(project_metadata_tidy$rin)) ) {
    #  project_metadata_tidy <- project_metadata_tidy %>%
    #    filter(rin >= 6.0) ## Only fresh-frozen preserved tissues
    #}
  
    
  } else if (data.source == "data_sources/tcga") {
    
    
      project_metadata_tidy <- project.metadata %>%
        as_tibble()
      
      # project_metadata_tidy %>% dplyr::count(tcga.cgc_sample_sample_type)
      # project_metadata_tidy[1:2,] %>% gather(key = "rail_id") %>% view
      # project_metadata_tidy$tcga.gdc_cases.samples.sample_type_id
      
      if ( any(is.na(project_metadata_tidy$tcga.cgc_sample_sample_type)) ) {
        project_metadata_tidy <- project_metadata_tidy[-which(is.na(project_metadata_tidy$tcga.cgc_sample_sample_type)),]
      }
      if ( any(project_metadata_tidy$tcga.cgc_sample_sample_type == "") ) {
        
        message(project_metadata_tidy$recount_project.project %>% unique," has empty tcga.cgc_sample_sample_type! Replacing...")
         
        idx <- which(project_metadata_tidy$tcga.cgc_sample_sample_type == "")
        project_metadata_tidy[idx,"tcga.cgc_sample_sample_type"] <- project_metadata_tidy[idx,"tcga.gdc_cases.samples.sample_type"][[1]]
      
        if (any(project_metadata_tidy$tcga.cgc_sample_sample_type == "") ) {
          message("ERROR...still some clusters are detected as string empty ''")
          break;
        }
         
      }
      
      # project_metadata_tidy %>% dplyr::count(tcga.cgc_sample_sample_type)
      
      project_metadata_tidy <- data.frame(sample_id =  project_metadata_tidy$external_id %>% as.character(),
                                          age = (project_metadata_tidy$tcga.gdc_cases.demographic.year_of_death - 
                                                  project_metadata_tidy$tcga.gdc_cases.demographic.year_of_birth) %>% as.character(),
                                          #rin = project_metadata_tidy$gtex.smrin %>% as.double(),
                                          gender = project_metadata_tidy$tcga.gdc_cases.demographic.gender %>% str_to_sentence() %>% as.character(),
                                          ethnicity =  project_metadata_tidy$tcga.gdc_cases.demographic.ethnicity %>% str_to_sentence() %>% as.character(),
                                          race = project_metadata_tidy$tcga.gdc_cases.demographic.race %>% str_to_sentence() %>% as.character(),
                                          
                                          individual_id = project_metadata_tidy$tcga.gdc_cases.submitter_id %>% as.character(),
                                          
                                          cluster = project_metadata_tidy$tcga.cgc_sample_sample_type %>% str_to_sentence() %>% as.character(),
                                          tumor = project_metadata_tidy$tcga.gdc_cases.project.name %>% str_to_sentence() %>% as.character(),
                                          site = project_metadata_tidy$tcga.cgc_case_primary_site %>% str_to_sentence() %>% as.character(),
                                          sample_type = project_metadata_tidy$tcga.gdc_cases.samples.sample_type_id ,
                                          
                                          smnabtcht = project_metadata_tidy$tcga.cgc_case_batch_number,
                                          #smafrze = project_metadata_tidy$gtex.smafrze,
                                          avg_read_length = project_metadata_tidy$recount_qc.star.average_input_read_length,
                                          all_mapped_reads = project_metadata_tidy$recount_qc.star.all_mapped_reads,
                                          SRA_project = project_metadata_tidy$recount_project.project) 
      
  }
  
  
  return(project_metadata_tidy)
  
}
