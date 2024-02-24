#####################################
## init.R file
## Contains the functions to download, qc, pair
## and database the exon-exon split read data 
## from a RNA-sequencing data project from recount3
#####################################

library(tidyverse)
library(GenomicRanges)
library(DBI)

# source("/home/sruiz/PROJECTS/recount3-database-project/init.R")


base_folder <- here::here()

print(base_folder)

# base_folder <- "/mnt/PROJECTS/splicing-accuracy-manuscript/"
dependencies_folder <- paste0("~/PROJECTS/splicing-accuracy-manuscript/dependencies/")


## Array with the IDs of the project to database as provided by recount3
## recount3 projects ID be checked here: https://jhubiostatistics.shinyapps.io/recount3-study-explorer/
recount3_project_IDs <- c("SRP100948") # SRP058181 # SRP100948


## This distinction is made for recount3 project that, albeit are separated in multiple independent projects (such as GTEx) in recount3,
## they belong to the same main project
## For instance, GTEx is stored in recount3 in multiple independent projects ID (e.g. BRAIN, SKIN, BLOOD, etc), but all of them
## belong to GTEx
supportive_reads <- 1
project_name <- paste0(recount3_project_IDs, "_", supportive_reads, "read") #SRP100948_

data_source <- "data_sources/sra"

#####################################
## LOAD SOURCE SCRIPTS
#####################################

setwd(file.path("~/PROJECTS/splicing-accuracy-manuscript/", "/scripts"))
files.sources = list.files()
sapply(files.sources, source)
setwd(file.path(base_folder))


#####################################
## INIT - DATABASE RECOUNT3 PROJECT
#####################################

## This is the Ensembl gtf transcriptome version 
gtf_versions <- c(110)


for (gtf_version in gtf_versions) {
  
  # gtf_version <- gtf_versions[1]
  
  database_folder <- paste0(base_folder, "/database/", project_name, "/", gtf_version, "/")
  results_folder <- paste0(base_folder, "/results/", project_name, "/", gtf_version, "/")
  #levelqc1_folder <- paste0(base_folder, "/database/")
  tpm_folder <- paste0(base_folder, "/results/tpm/")
  

  # download_recount3_data(recount3.project.IDs = recount3_project_IDs,
  #                        project.name = project_name,
  #                        gtf.version = gtf_version,
  #                        data.source = data_source,
  #                        database.folder = database_folder,
  #                        results.folder = results_folder)

  prepare_recount3_data(recount3.project.IDs = recount3_project_IDs,
                        data.source = data_source,
                        results.folder = results_folder,
                        subsampling = T,
                        levelqc1.folder = database_folder,
                        supporting.reads = supportive_reads,
                        num.cores = 2)

  junction_pairing(recount3.project.IDs = recount3_project_IDs,
                   results.folder = results_folder,
                   replace = T,
                   num.cores = 2)

  get_all_annotated_split_reads(recount3.project.IDs = recount3_project_IDs,
                                database.folder = database_folder,
                                results.folder = results_folder)

  get_all_raw_jxn_pairings(recount3.project.IDs = recount3_project_IDs,
                           database.folder = database_folder,
                           results.folder = results_folder)


  recount3_project_IDs <- readRDS(file = paste0(results_folder, "/all_final_projects_used.rds"))


  tidy_data_pior_sql(recount3.project.IDs = recount3_project_IDs,
                     database.folder = database_folder,
                     levelqc1.folder = database_folder,
                     results.folder = results_folder)


  generate_transcript_biotype_percentage(gtf.version = gtf_version,
                                         database.folder = database_folder,
                                         results.folder = results_folder)


  generate_recount3_tpm(recount3.project.IDs = recount3_project_IDs,
                        data.source = data_source,
                        tpm.folder = tpm_folder,
                        results.folder = results_folder)


  database_path <- paste0(database_folder,  "/", project_name, ".sqlite")

  sql_database_generation(database.path = database_path,
                          recount3.project.IDs = recount3_project_IDs,
                          remove.all = T,
                          database.folder = database_folder,
                          results.folder = results_folder,
                          gtf.version = gtf_version)
}

