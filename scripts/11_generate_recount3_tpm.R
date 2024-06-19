library(tidyverse)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(DBI)
library(dplyr)
library(here)
library(doParallel)
library(tidyverse)
library(protr)

#' Title
#' Obtains the raw counts data from recount3 for a given projectID. 
#' It transforms the raw counts and calculates the TPM of the genes across all samples from the projectID.
#' As the raw counts represent the total coverage for a given gene within a given sample, raw counts need to be scaled by library size 
#' as longer genes may obtain higher number of counts, which does not mean higher expression.
#'
#' Then, it separates the TPM values across the samples of each sample cluster.
#' @param recount3.project.IDs Vector with all the recount3 projectIDs to process
#' @param data.source source of the data within recount3. It can be:
#' - "data_sources/sra"  
#' - "data_sources/gtex" 
#' - "data_sources/tcga"
#' @param results.folder Path to the local folder storing the results produced
#'
#' @return
#' @export
#'
#' @examples
generate_recount3_tpm <- function(recount3.project.IDs,
                                  data.source,
                                  tpm.folder,
                                  results.folder) {
  
  
  ## Loop through the recount3 projects received by parameter
  # doParallel::registerDoParallel(5)
  
  for(i in seq(length(recount3.project.IDs)) ) {
  # Register cluster
  #cluster <- makeCluster(5)
  #registerDoParallel(cluster)
  #foreach(i = seq(length(recount3.project.IDs))) %dopar% {
    # i <- 5
    # project_id <- "COAD"
    
    project_id <- recount3.project.IDs[i]
    results_folder_local <- paste0(results.folder, "/", project_id, "/")
    results_folder_local_tpm <- paste0(results_folder_local, "/tpm/")
    dir.create(file.path(results_folder_local_tpm), recursive = TRUE, showWarnings = F)
    dir.create(file.path(tpm.folder), recursive = TRUE, showWarnings = F)
    recount_tpm <- NULL
    
    if ( !file.exists(paste0(tpm.folder, "/", project_id, "_tpm.rds")) ) {
    
      logger::log_info("Downloading raw counts from ", project_id, "...") 
      
      ## 1. Get expression data from recount3, transform raw counts and calculate TPM
      rse <- NULL
      tryCatch(
        eval(
          rse <- recount3::create_rse_manual(
            project = project_id,
            project_home = data.source,
            organism = "human",
            annotation = "gencode_v29",
            type = "gene"
          )
          
         ),
        error = function(e) {
          rse = NULL
        }
      )
      
      ## The counts download has been successful, so we calculate the TPM
      if ( !is.null(rse) ) {
        message(project_id)
        SummarizedExperiment::assays(rse)$counts <- recount3::transform_counts(rse)
        # SummarizedExperiment::assayNames(rse) %>% logger::log_info()
        
        logger::log_info("Computing TPM for genes found in ", project_id)
        recount_tpm <- recount::getTPM(rse)

        ## Save tpm values for all genes across all samples
        saveRDS(object = recount_tpm,
                file = paste0(tpm.folder, "/", project_id, "_tpm.rds"))
        
        ## Release some memory
        rm(rse)
        gc()
      }
      
    } else {

      logger::log_info("Loading TPM data for ", project_id, "...")
      recount_tpm <- readRDS(file = paste0(tpm.folder, "/", project_id, "_tpm.rds"))

    }
    
    
    ## 2. For each tissue within the current project, filter the RSE by its samples

    if ( !is.null(recount_tpm) && 
         file.exists(paste0(results_folder_local, "/base_data/", project_id, "_samples_metadata.rds")) &&
         file.exists(paste0(results_folder_local, "/base_data/", project_id, "_clusters_used.rds")) ) {

      metadata.info <- readRDS(file = paste0(results_folder_local, "/base_data/", project_id, "_samples_metadata.rds"))
      clusters_ID <- readRDS(file = paste0(results_folder_local, "/base_data/", project_id, "_clusters_used.rds")) %>% unique()

      for ( cluster_id in clusters_ID ) {

        # cluster_id <- clusters_ID[1]
        logger::log_info("Getting results for ", cluster_id)

        cluster_samples <- readRDS(file = paste0(results_folder_local, "/base_data/",
                                                 project_id, "_", cluster_id, "_samples_used.rds"))

        ## Filter the object for the set of samples corresponding to the current cluster
        recount_tpm_local <- recount_tpm %>% 
          as_tibble(rownames = "gene_id") %>%
          #rename_with(~stringr::str_replace_all(., '\\.1', '')) %>%
          mutate(gene_id = gsub(pattern = "\\..*",replacement = "",x = gene_id)) %>% 
          dplyr::select(c("gene_id", all_of(cluster_samples)))

        ## Save results
        saveRDS(object = recount_tpm_local,
                file = paste0(results_folder_local_tpm, project_id, "_", cluster_id, "_tpm.rds"))

        rm(recount_tpm_local)
        rm(cluster_samples)
        gc()

      }
    }


    ## 3. We do the same for the ageing analysis

    if ( !is.null(recount_tpm) && 
         file.exists(paste0(results_folder_local, "/base_data/", project_id, "_age_samples_metadata.rds")) &&
         file.exists(paste0(results_folder_local, "/base_data/", project_id, "_age_clusters_used.rds")) ) {

      metadata.info <- readRDS(file = paste0(results_folder_local, "/base_data/", project_id, "_age_samples_metadata.rds"))
      clusters_ID <- readRDS(file = paste0(results_folder_local, "/base_data/", project_id, "_age_clusters_used.rds")) %>% unique()

      for ( cluster_id in clusters_ID ) {

        # cluster_id <- clusters_ID[1]
        logger::log_info("Getting results for ", cluster_id)

        cluster_samples <- readRDS(file = paste0(results_folder_local, "/base_data/",
                                                 project_id, "_", cluster_id, "_samples_used.rds"))

        ## Filter the object for the set of samples corresponding to the current cluster
        recount_tpm_local <- recount_tpm %>% #head %>%
          as_tibble(rownames = "gene_id") %>%
          #rename_with(~stringr::str_replace_all(., '\\.1', '')) %>%
          mutate(gene_id = gsub(pattern = "\\..*",replacement = "",x = gene_id)) %>%
          dplyr::select(c("gene_id", all_of(cluster_samples)))

        recount_tpm_local %>% head %>% logger::log_info()

        ## Save results
        saveRDS(object = recount_tpm_local,
                file = paste0(results_folder_local_tpm, project_id, "_", cluster_id, "_tpm.rds"))

        rm(recount_tpm_local)
        rm(cluster_samples)
        gc()

      }
    }

    logger::log_info(project_id, " finished!") 
    

    rm(recount_tpm)
    gc()
    
  }
  #stopCluster(cl = cluster)
}

# TCGA -------------------------------------------------------------------------------------------------------------

# recount3.project.IDs <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
#                           "KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD",
#                           "LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM") %>%  sort()
# data.source <- "data_sources/tcga"
# tpm.folder <- "/mnt/PROJECTS/recount3-database-project/results/TCGA_1read_subsampleFALSE/tpm/"
# results.folder <- "/mnt/PROJECTS/recount3-database-project/results/TCGA_1read_subsampleFALSE/111/"

# PD/control --------------------------------------------------------------------------------------------------------

# recount3.project.IDs <- "SRP058181"
# data.source <- "data_sources/sra"
# tpm.folder <- "/mnt/PROJECTS/recount3-database-project/results/tpm/"
# results.folder <- "/mnt/PROJECTS/recount3-database-project/results/SRP058181_1read/110/"

# AD/control --------------------------------------------------------------------------------------------------------

# recount3.project.IDs <- "SRP100948"
# data.source <- "data_sources/sra"
# tpm.folder <- "/mnt/PROJECTS/recount3-database-project/results/tpm/"
# results.folder <- "/mnt/PROJECTS/recount3-database-project/results/SRP100948_1read/110/"

# GTEx --------------------------------------------------------------------------------------------------------------
# source("/mnt/PROJECTS/splicing-accuracy-manuscript/scripts/11_generate_recount3_tpm.R")

#
#
# recount3.project.IDs <- c( "ADIPOSE_TISSUE",  "ADRENAL_GLAND",   "BLADDER",         "BLOOD",           "BLOOD_VESSEL",    "BONE_MARROW",
#                            "BRAIN",           "BREAST",          "CERVIX_UTERI",    "COLON",           "ESOPHAGUS",       "FALLOPIAN_TUBE",
#                            "HEART",           "KIDNEY",          "LIVER",           "LUNG",            "MUSCLE",          "NERVE",
#                            "OVARY",           "PANCREAS",        "PITUITARY",       "PROSTATE",        "SALIVARY_GLAND",  "SKIN",
#                            "SMALL_INTESTINE", "SPLEEN",          "STOMACH",         "TESTIS",          "THYROID",         "UTERUS",
#                            "VAGINA"  )
#
# ## This is the name of the project producing the database
# data.source <- "data_sources/gtex"
# tpm.folder <- "/mnt/PROJECTS/recount3-database-project/results/GTEX_1read_subsampleFALSE/tpm/"
# results.folder <- "/mnt/PROJECTS/recount3-database-project/results/GTEX_1read_subsampleFALSE/111/"
