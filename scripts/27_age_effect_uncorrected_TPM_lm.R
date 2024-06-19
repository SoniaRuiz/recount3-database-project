#' Title
#'
#' @param project.id 
#' @param tpm.uncorrected 
#' @param sample.metadata 
#' @param gene.list 
#' @param results.folder 
#'
#' @return
#' @export
#'
#' @examples
age_effect_uncorrected_TPM_lm <- function(project.id,
                                          tpm.uncorrected,
                                          sample.metadata,
                                          gene.list,
                                          results.folder) {
  
  ################################
  ## LOAD GENES OF INTEREST
  ################################
  
  ## Check the samples used are the same
  if ( !( identical(x = colnames(tpm.uncorrected)[-1] %>% unique %>% sort(),
                    y = colnames(sample.metadata)[-1] %>% unique %>% sort())) ) {
    print("ERROR! The columns of the two objects provided do not match!")
    break;
  }
  
  
  genes_of_interest <- (tpm.uncorrected$gene) %>% unique
  
  
  ## Per each gene, obtain TPM and metadata values and run lm
  df_lm_output <- map_df(genes_of_interest, function(gene_id) {
    
    # gene_id <- genes_of_interest[1]
    # gene_id <- "ENSG00000167281"
    
    # print(gene_id)
    
    tpm <- tpm.uncorrected %>%
      filter(gene == gene_id) %>%
      mutate(gene = "tpm") %>%
      as_tibble()
    
    ## If the gene is expressed across the samples of the current tissue
    if ( rowSums(tpm[, c(2:ncol(tpm))], na.rm = T) != 0 ) {
      
      
      ## Prepare the data prior modelling
      df_gene <- rbind(tpm %>% dplyr::rename(covariates = "gene"), sample.metadata) %>%
        gather(sample, value, -covariates)  %>%
        drop_na() %>%
        spread(key = covariates, value = value) %>%
        drop_na() %>%
        dplyr::mutate(gtex.age = gtex.age %>% as.double(),
                      gtex.dthhrdy = gtex.dthhrdy %>% as.double(),
                      gtex.sex = gtex.sex %>% as.double(),
                      gtex.smcenter = gtex.smcenter %>% as.double(),
                      gtex.smgebtch = gtex.smgebtch %>% as.double(),
                      gtex.smgebtchd = gtex.smgebtchd %>% as.double(),
                      gtex.smnabtch = gtex.smnabtch %>% as.double(),
                      gtex.smnabtchd = gtex.smnabtchd %>% as.double(),
                      gtex.smnabtcht = gtex.smnabtcht %>% as.double(),
                      gtex.smrin = gtex.smrin %>% as.double(),
                      ## Transform TPM to log10 so resilduals after modelling in the linear model are normally distributed
                      ## https://stats.stackexchange.com/questions/40907/in-regression-analysis-what-does-taking-the-log-of-a-variable-do
                      tpm = (tpm %>% as.double()) %>% log10) %>%
        as_tibble() %>%
        filter_all(all_vars(!is.infinite(.))) %>%
        drop_na()
      
      
      ## Linear models to predict the log10 TPM value of each RBP
      lm_output <- lm(tpm ~ 
                        gtex.smrin + 
                        gtex.smcenter +
                        gtex.smgebtch +
                        gtex.smgebtchd +
                        gtex.smnabtch +
                        gtex.smnabtchd +
                        gtex.smnabtcht +
                        gtex.dthhrdy +
                        gtex.sex +
                        gtex.smtsd +
                        gtex.age,
                      data = df_gene)
      
      lm_output <- lm_output %>% summary()
      df_lm_output <- lm_output$coefficients %>% as_tibble(rownames = "covariate")
      
      df_lm_output <- df_lm_output %>%
        mutate(gene_id = gene_id) %>%
        dplyr::select(covariate, Estimate, pval = `Pr(>|t|)`, gene_id)
  
      return(df_lm_output)
      
    } else {
      return(NULL)
    }
    
  }) 
  
  ## Filter by age covariate and order by estimate size
  df_lm_output_tidy <- df_lm_output %>%
    mutate(q = p.adjust(p = pval, method = "fdr"))  %>%
    arrange(Estimate)
  
  ## Add gene SYMBOL info
  df_lm_output_tidy <- left_join(x = df_lm_output_tidy, 
                                 y = gene.list %>% drop_na(), 
                                 by = c("gene_id" = "id")) %>%
    group_by(gene_id) 
  
  ## Save results
  saveRDS(object = df_lm_output_tidy,
          file = file.path(results.folder, paste0(project.id, "_tpm_lm_ensembl105.rds")))
  
  ###############################################
  ## Get GENE that decrease expression with AGE
  ###############################################
  
  ## Genes with expression levels affected by age
  genes_age_lm <- df_lm_output_tidy %>%
    filter(covariate == "gtex.age") %>%
    arrange(name)
  
  return(genes_age_lm %>%
           mutate(project = project.id))
  
  
}