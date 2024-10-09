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

#####################################
## SET MAIN VARIABLES
#####################################

## This is the Ensembl gtf transcriptome version 
gtf_version <- c(111)

## This is the name of the project producing the database
supportive_reads <- 1
data_subsample = F
main_project_identifier <- "SRP181886" 
project_name <- paste0(main_project_identifier, "_", supportive_reads, "read_subsample", data_subsample)


args <-
  list(
    base_folder = here::here(),
    dependencies_folder = file.path(here::here(), "dependencies"),
    #data_source = "data_sources/tcga" 
    #data_source = "data_sources/gtex"
    data_source = "data_sources/sra",
    recount3_project_IDs ="SRP181886",
    gtf_version = gtf_version,
    database_folder = file.path(here::here(), "database", project_name, gtf_version),
    results_folder = file.path(here::here(), "results", project_name, gtf_version),
    figures_folder = file.path(here::here(), "results", project_name, gtf_version),
    tpm_folder = file.path(here::here(), "results", project_name, "tpm")
  )


dir.create(path = args$database_folder, recursive = T, showWarnings = F)
dir.create(path = args$results_folder, recursive = T, showWarnings = F)
dir.create(path = args$figures_folder, recursive = T, showWarnings = F) 
dir.create(path = args$tpm_folder, recursive = T, showWarnings = F)


## Example other recount3 project identifiers

#main_project_identifier <- "GTEX"
#main_project_identifier <- "TCGA"
#main_project_identifier <- "SRP151040"
#main_project_identifier <- "SRP100948"
#main_project_identifier <- "SRP058181"

## Array with the IDs of the project to database as provided by recount3
## recount3 projects ID be checked here: https://jhubiostatistics.shinyapps.io/recount3-study-explorer/

# recount3_project_IDs <- c( "ADIPOSE_TISSUE",  "ADRENAL_GLAND",   "BLADDER",         "BLOOD",           "BLOOD_VESSEL",    "BONE_MARROW",
#                            "BRAIN",           "BREAST",          "CERVIX_UTERI",    "COLON",           "ESOPHAGUS",       "FALLOPIAN_TUBE",
#                            "HEART",           "KIDNEY",          "LIVER",           "LUNG",            "MUSCLE",          "NERVE",
#                            "OVARY",           "PANCREAS",        "PITUITARY",       "PROSTATE",        "SALIVARY_GLAND",  "SKIN",
#                            "SMALL_INTESTINE", "SPLEEN",          "STOMACH",         "TESTIS",          "THYROID",         "UTERUS",
#                            "VAGINA"  )

# recount3_project_IDs <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
#                          "KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD",
#                          "LUSC","MESO","OV","PAAD","PCPG","PRAD","READ",
#                          "SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM") %>%  sort()

#recount3_project_IDs <-"SRP151040"
#recount3_project_IDs <-"SRP100948"
#recount3_project_IDs <-"SRP058181"
#recount3_project_IDs <-"SRP181886" # Extended AD/control


#####################################
## LOAD SOURCE SCRIPTS
#####################################

# scripts_folder <- paste0("/mnt/PROJECTS/splicing-accuracy-manuscript/scripts/")

scripts_folder <- file.path(args$base_folder, "scripts/")
setwd(file.path(scripts_folder))
files.sources = list.files()
sapply(files.sources, source)
setwd(file.path(args$base_folder))



#####################################
## INIT LOGS
#####################################

log_file <- here::here(paste0("logs/Splicing_Database_Generation_",project_name,".log"))
logger::log_appender(logger::appender_tee(log_file, append = T))
logger_layout <- logger::layout_glue_generator(format = '[{time}] [{level}] {msg}')
logger::log_layout(logger_layout)

## attr(attr(loadEdb,"srcref"),"srcfile")

#####################################
## INIT - DATABASE RECOUNT3 PROJECT
#####################################

## This is the Ensembl gtf transcriptome version 
#gtf_versions <- c(111)


#for (gtf_version in gtf_versions) {
  
  # gtf_version <- gtf_versions[1]
  
  # database_folder <- file.path(args$base_folder, "database", project_name, gtf_version)
  # dir.create(path = database_folder, recursive = T, showWarnings = F)
  # 
  # results_folder <- file.path(args$base_folder, "results", project_name, gtf_version)
  # dir.create(path = results_folder, recursive = T, showWarnings = F)
    

  #################################################
  ## DOWNLOAD AND PREPARE JUNCTIONS FROM RECOUNT3 PROJECTS
  
  
  DownloadRecount3Data(recount3.project.IDs = args$recount3_project_IDs,
                       project.name = project_name,
                       gtf.version = args$gtf_version,
                       blacklist.path = file.path(args$dependencies_folder, "hg38-blacklist.v2.bed"),
                       gtf.path = file.path(args$dependencies_folder, paste0("/Homo_sapiens.GRCh38.", args$gtf_version, ".chr.gtf")),
                       data.source = args$data_source,
                       database.folder = args$database_folder,
                       results.folder = args$results_folder)

  PrepareRecount3Data(recount3.project.IDs = args$recount3_project_IDs,
                      data.source = args$data_source,
                      results.folder = args$results_folder,
                      subsampling = data_subsample,
                      levelqc1.folder = args$database_folder,
                      supporting.reads = supportive_reads,
                      num.cores = 4)

  #################################################
  ## JUNCTION PAIRING
  
  JunctionPairing(recount3.project.IDs = args$recount3_project_IDs,
                   results.folder = args$results_folder,
                   replace = T,
                   num.cores = 8)
  
  GetAllAnnotatedSplitReads(recount3.project.IDs = args$recount3_project_IDs,
                            database.folder = args$database_folder,
                            results.folder = args$results_folder,
                            num.cores = 8)
  
  GetAllRawJxnPairings(recount3.project.IDs = args$recount3_project_IDs,
                           database.folder = args$database_folder,
                           results.folder = args$results_folder,
                           num.cores = 8)
  
  #################################################
  ## DATA BASE
  
  recount3_project_IDs_used <- readRDS(file = file.path(args$results_folder, "all_final_projects_used.rds"))
  
  
  TidyDataPiorSQL(recount3.project.IDs = recount3_project_IDs_used,
                  database.folder = args$database_folder,
                  levelqc1.folder = args$database_folder,
                  results.folder = args$results_folder)
  
  
  GenerateTranscriptBiotypePercentage(gtf.version = args$gtf_version,
                                      database.folder = args$database_folder,
                                      results.folder = args$results_folder)
  
  
  GenerateRecount3TPM(recount3.project.IDs = recount3_project_IDs_used,
                      data.source = data_source,
                      tpm.folder = args$tpm_folder,
                      results.folder = args$results_folder)
  
  
  
  database_sqlite_path <- paste0(args$database_folder,  "/", project_name, ".sqlite")
  
  SqlDatabaseGeneration(database.path = database_sqlite_path,
                        recount3.project.IDs = recount3_project_IDs_used,
                        remove.all = F,
                        database.folder = args$database_folder,
                        results.folder = args$results_folder,
                        gtf.version = args$gtf_version,
                        dependencies.folder = args$dependencies_folder,
                        discard.minor.introns = T)
  
#}

