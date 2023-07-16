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


## Array with the IDs of the project to database as provided by recount3
## recount3 projects ID be checked here: https://jhubiostatistics.shinyapps.io/recount3-study-explorer/
recount3_project_IDs <- c("SRP058181")#,"SRP100948")
## This distinction is made for recount3 project that, albeit are separated in multiple independent projects (such as GTEx) in recount3,
## they belong to the same main project
## For instance, GTEx is stored in recount3 in multiple independent projects ID (e.g. BRAIN, SKIN, BLOOD, etc), but all of them
## belong to GTEx
project_name <- c("SRP058181")#,"SRP100948")
data_source <- "data_sources/sra"


setwd(normalizePath("."))

dependencies_folder <- paste0("/home/sruiz//PROJECTS/splicing-accuracy-manuscript/dependencies/")


source(paste0(getwd(), "/database_qc_data.R"))
source(paste0(getwd(), "/database_junction_pairing.R"))
source(paste0(getwd(), "/database_SQL_helper.R"))
source(paste0(getwd(), "/database_SQL_generation.R"))


#####################################
## INIT - DATABASE RECOUNT3 PROJECT
#####################################

## This is the Ensembl gtf transcriptome version 
gtf_versions <- c(105)


for (gtf_version in gtf_versions) {
  
  # gtf_version <- gtf_versions[1]
  
  database_folder <- paste0(getwd(), "/database/", project_name, "/", gtf_version, "/")
  results_folder <- paste0(getwd(), "/results/", project_name, "/", gtf_version, "/")
  
  init_recount3_data(recount3.project.IDs = recount3_project_IDs,
                    project.name = project_name,
                    gtf.version = gtf_version,
                    data.source = data_source,
                    database.folder = database_folder,
                    results.folder = results_folder)
  
  tidy_sample_cluster_data(recount3.project.IDs = recount3_project_IDs,
                           project.name = project_name,
                           data.source = data_source,
                           gtf.version = gtf_version,
                           database.folder = database_folder,
                           results.folder = results_folder)

  junction_pairing(recount3.project.IDs = recount3_project_IDs,
                   project.name = project_name,
                   gtf.version = gtf_version,
                   database.folder = database_folder,
                   results.folder = results_folder)


  get_all_annotated_split_reads(recount3.project.IDs = recount3_project_IDs,
                                project.name = project_name,
                                gtf.version = gtf_version,
                                database.folder = database_folder,
                                results.folder = results_folder)


  get_all_raw_distances_pairings(recount3.project.IDs = recount3_project_IDs,
                                 project.name = project_name,
                                 gtf.version = gtf_version,
                                 database.folder = database_folder,
                                 results.folder = results_folder)


  tidy_data_pior_sql(recount3.project.IDs = recount3_project_IDs,
                     gtf.version = gtf_version,
                     database.folder = database_folder,
                     results.folder = results_folder)


  generate_transcript_biotype_percentage(recount3.project.IDs = recount3_project_IDs,
                                         project.name = project_name,
                                         gtf.version = gtf_version,
                                         database.folder = database_folder,
                                         results.folder = results_folder)

  generate_recount3_tpm(recount3.project.IDs = recount3_project_IDs,
                        project.name = project_name,
                        gtf.version = gtf_version,
                        data.source = data_source,
                        database.folder = database_folder,
                        results.folder = results_folder)
  
  

  database_path <- paste0(database_folder,  "/", project_name, ".sqlite")
  
  sql_database_generation(database.path = database_path,
                          recount3.project.IDs = recount3_project_IDs,
                          project.name = project_name,
                          gtf.version = gtf_version,
                          remove.all = F,
                          database.folder = database_folder,
                          results.folder = results_folder)
  
  gc()
}

