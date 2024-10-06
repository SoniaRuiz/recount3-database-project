library(tidyverse)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(DBI)
library(dplyr)
library(here)
library(logger)
library(foreach)
library(doParallel)
library(Biostrings)
library(tidyverse)
library(protr)
library(optparse)

# source("/home/sruiz/PROJECTS/recount3-database-project/init_ENCODE.R")

# base_folder <- "/mnt/PROJECTS/splicing-accuracy-manuscript/"
# base_folder <- "~/PROJECTS/splicing-accuracy-manuscript/"



#####################################
## LOAD SOURCE SCRIPTS
#####################################

setwd(file.path(base_folder,"scripts"))
files.sources = list.files()
sapply(files.sources, source)
setwd(file.path(base_folder))


#####################################
## SET MAIN VARIABLES
#####################################

## This is the Ensembl gtf transcriptome version 
gtf_version <- c(111)

## This is the name of the project producing the database
supportive_reads <- 1
analysis_type = "shRNA"
main_project <- paste0("ENCODE_SR_", supportive_reads, "read_", analysis_type)

args <-
  list(
    base_folder = here::here(),
    dependencies_folder = file.path(here::here(), "dependencies"),
    
    results_folder = file.path(here::here(), "results", main_project, gtf_version),
    figures_folder = file.path(here::here(), "results", main_project, gtf_version)
  )


dir.create(path = args$results_folder, recursive = T, showWarnings = F)
dir.create(path = args$figures_folder, recursive = T, showWarnings = F)


log_file <- here::here("logs/", paste0("Splicing_Analysis_ENCODE_",gtf_version,"_",analysis_type,".log"))
logger::log_appender(logger::appender_tee(log_file, append = T))
logger_layout <- logger::layout_glue_generator(format = '[{time}] [{level}] {msg}')
logger::log_layout(logger_layout)

#######################################################
## Load the ENCODE Metadata per RBP
## If the file does not exists, download the data
#######################################################


metadata_path <- paste0(args$results_folder, "/metadata_",analysis_type,"_samples.tsv")

if ( !file.exists(metadata_path) ) {
  
  ## 1. Download metadata per RBP from ENCODE
  metadata <- ENCODE_download_metadata(experiment_type = analysis_type,
                                       results_path = args$results_folder,
                                       dependencies_path = args$dependencies_folder)
  
  ## 2. Download .bam files
  ENCODE_download_bams(metadata_df = metadata, results_path = args$results_folder,
                       regtools_path = "/home/grocamora/tools/regtools/build/",
                       samtools_path = "/home/grocamora/tools/samtools/bin/")
  
  
} else {
  
  metadata <- readr::read_delim(metadata_path, show_col_types = F)
}

#######################################################
## Get the RBPs
## If it is an shRNA project, only focus on RBPs from
## Van Nostrand et al (2020)
#######################################################


if (analysis_type == "shRNA") {
  metadata <- metadata %>%
    dplyr::filter(if_any(c("Splicing regulation", Spliceosome, "Novel RBP", "Exon Junction Complex", NMD), ~ . != 0))
}

target_RBPs <- metadata %>% dplyr::pull(target_gene) %>% unique()

#target_RBPs <- target_RBPs[-38]
metadata_RBPs <- metadata %>% dplyr::filter(target_gene %in% target_RBPs)


#####################################
## INIT - DATABASE PROJECT
#####################################

for (gtf_version in gtf_versions) {
  
  # gtf_version <- gtf_versions[1]
  
  #################################################
  ## SET UP VARIABLES
  
  database_base_folder <- paste0(base_folder, "/database/", project_name, "/")
  dir.create(database_base_folder, recursive = T, showWarnings = F)
  database_folder <- paste0(base_folder, "/database/", project_name, "/", gtf_version)
  dir.create(database_folder, recursive = T, showWarnings = F)
  
  results_folder <- paste0(base_folder, "/results/", project_name, "/", gtf_version)
  dir.create(results_folder, recursive = T, showWarnings = F)
  
  levelqc1_folder <- database_folder
  dir.create(levelqc1_folder, recursive = T, showWarnings = F)
  
  tpm_folder <- paste0(base_folder, "/results/", project_name, "/tpm/")
  dir.create(tpm_folder, recursive = T, showWarnings = F)

  
  #################################################
  ## PREPARE JUNCTIONS FROM KNOCKDOWN EXPERIMENTS

  source(here::here("scripts/28_ENCODE_prepare_encode_data.R"))

  logger::log_info(paste0(Sys.time(), "\t\t starting 'prepare_encode_data' function..."))
  prepare_encode_data(metadata = metadata_RBPs,
                      RBP_source_path = RBPs_base_folder,
                      results_path = results_folder,
                      database_path = database_base_folder,
                      gtf_version,
                      gtf_path = file.path(dependencies_folder, paste0("Homo_sapiens.GRCh38.",gtf_version,".chr.gtf")),
                      ENCODE_silencing_series = analysis_type,
                      num_cores = 8)
  gc()

  #################################################
  ## JUNCTION PAIRING

  junction_pairing(recount3.project.IDs = target_RBPs,
                   results.folder = results_folder,
                   replace = T,
                   num.cores = 10)


  get_all_annotated_split_reads(recount3.project.IDs = target_RBPs,
                                database.folder = database_folder,
                                results.folder = results_folder,
                                num.cores = 10)


  get_all_raw_jxn_pairings(recount3.project.IDs = target_RBPs,
                           database.folder = database_folder,
                           results.folder = results_folder,
                           num.cores = 10)

  
  #################################################
  ## DATA BASE

  all_final_projects_used <- readRDS(file.path(RBPs_base_folder, "all_final_projects_used.rds"))


  tidy_data_pior_sql(recount3.project.IDs = all_final_projects_used,
                     database.folder = database_folder,
                     levelqc1.folder = levelqc1_folder,
                     results.folder = results_folder)


  generate_transcript_biotype_percentage(gtf.version = gtf_version,
                                         database.folder = database_folder,
                                         results.folder = results_folder)


  database_path <- paste0(database_folder,  "/", project_name, ".sqlite")


  sql_database_generation(database.path = database_path,
                          recount3.project.IDs = all_final_projects_used,
                          remove.all = T,
                          database.folder = database_folder,
                          results.folder = results_folder,
                          gtf.version = gtf_version,
                          discard.minor.introns = F)

  gc()
}
