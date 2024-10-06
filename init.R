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

# base_folder <- "/mnt/PROJECTS/splicing-accuracy-manuscript/"
base_folder <- here::here()

#print(base_folder)


# dependencies_folder <- paste0("/mnt/PROJECTS/splicing-accuracy-manuscript//dependencies/")
dependencies_folder <- file.path(base_folder, "dependencies/")




## Array with the IDs of the project to database as provided by recount3
## recount3 projects ID be checked here: https://jhubiostatistics.shinyapps.io/recount3-study-explorer/


recount3_project_IDs <- c( "ADIPOSE_TISSUE",  "ADRENAL_GLAND",   "BLADDER",         "BLOOD",           "BLOOD_VESSEL",    "BONE_MARROW",
                           "BRAIN",           "BREAST",          "CERVIX_UTERI",    "COLON",           "ESOPHAGUS",       "FALLOPIAN_TUBE",
                           "HEART",           "KIDNEY",          "LIVER",           "LUNG",            "MUSCLE",          "NERVE",
                           "OVARY",           "PANCREAS",        "PITUITARY",       "PROSTATE",        "SALIVARY_GLAND",  "SKIN",
                           "SMALL_INTESTINE", "SPLEEN",          "STOMACH",         "TESTIS",          "THYROID",         "UTERUS",
                           "VAGINA"  )

# recount3_project_IDs <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
#                          "KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD",
#                          "LUSC","MESO","OV","PAAD","PCPG","PRAD","READ",
#                          "SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM") %>%  sort()

#recount3_project_IDs <- recount3_project_IDs[11:17]




#recount3_project_IDs <-"SRP151040"
#recount3_project_IDs <-"SRP100948"
#recount3_project_IDs <-"SRP058181"
#recount3_project_IDs <-"SRP181886" # Extended AD/control

## This distinction is made for recount3 project that, albeit are separated in multiple independent projects (such as GTEx) in recount3,
## they belong to the same main project
## For instance, GTEx is stored in recount3 in multiple independent projects ID (e.g. BRAIN, SKIN, BLOOD, etc), but all of them
## belong to GTEx

main_project_identifier <- "GTEX"
#main_project_identifier <- "TCGA"
#main_project_identifier <- "SRP151040"
#main_project_identifier <- "SRP100948"
#main_project_identifier <- "SRP058181"
#main_project_identifier <- "SRP181886" 


#data_source <- "data_sources/tcga" 
data_source <- "data_sources/gtex"
#data_source <- "data_sources/sra"


supportive_reads <- 1
data_subsample <- F
project_name <- paste0(main_project_identifier, "_", supportive_reads, "read_subsample", data_subsample)




#####################################
## LOAD SOURCE SCRIPTS
#####################################

# scripts_folder <- paste0("/mnt/PROJECTS/splicing-accuracy-manuscript/scripts/")

scripts_folder <- file.path(base_folder, "scripts/")
setwd(file.path(scripts_folder))
files.sources = list.files()
sapply(files.sources, source)
setwd(file.path(base_folder))

#####################################
## INIT LOGS
#####################################

log_file <- here::here(paste0("logs/Splicing_Database_Generation_",project_name,".log"))
logger::log_appender(logger::appender_tee(log_file, append = T))
logger_layout <- logger::layout_glue_generator(format = '[{time}] [{level}] {msg}')
logger::log_layout(logger_layout)

#####################################
## INIT - DATABASE RECOUNT3 PROJECT
#####################################

## This is the Ensembl gtf transcriptome version 
gtf_versions <- c(76, 81, 90, 104)


for (gtf_version in gtf_versions) {
  
  # gtf_version <- gtf_versions[1]
  
  database_folder <- file.path(base_folder, "database", project_name, gtf_version)
  dir.create(path = database_folder, recursive = T, showWarnings = F)
  
  results_folder <- file.path(base_folder, "results", project_name, gtf_version)
  dir.create(path = results_folder, recursive = T, showWarnings = F)
  
  tpm_folder <- file.path(base_folder, "results", project_name, "tpm")
  dir.create(path = tpm_folder, recursive = T, showWarnings = F)
  

  download_recount3_data(recount3.project.IDs = recount3_project_IDs,
                         project.name = project_name,
                         gtf.version = gtf_version,
                         data.source = data_source,
                         database.folder = database_folder,
                         results.folder = results_folder)

  prepare_recount3_data(recount3.project.IDs = recount3_project_IDs,
                        data.source = data_source,
                        results.folder = results_folder,
                        subsampling = data_subsample,
                        levelqc1.folder = database_folder,
                        supporting.reads = supportive_reads,
                        num.cores = 4)

  junction_pairing(recount3.project.IDs = recount3_project_IDs,
                   results.folder = results_folder,
                   replace = T,
                   num.cores = 8)

  get_all_annotated_split_reads(recount3.project.IDs = recount3_project_IDs,
                                database.folder = database_folder,
                                results.folder = results_folder,
                                num.cores = 8)

  get_all_raw_jxn_pairings(recount3.project.IDs = recount3_project_IDs,
                           database.folder = database_folder,
                           results.folder = results_folder,
                           num.cores = 8)


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
                          remove.all = F,
                          database.folder = database_folder,
                          results.folder = results_folder,
                          gtf.version = gtf_version)
  
}

