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
    
    database_folder = file.path(base_folder, "database", main_project, gtf_version),
    results_folder = file.path(here::here(), "results", main_project, gtf_version),
    figures_folder = file.path(here::here(), "results", main_project, gtf_version)
  )


dir.create(path = args$database_folder, recursive = T, showWarnings = F)
dir.create(path = args$results_folder, recursive = T, showWarnings = F)
dir.create(path = args$figures_folder, recursive = T, showWarnings = F)


log_file <- here::here("logs/", paste0("Splicing_Analysis_ENCODE_",gtf_version,"_",analysis_type,".log"))
logger::log_appender(logger::appender_tee(log_file, append = T))
logger_layout <- logger::layout_glue_generator(format = '[{time}] [{level}] {msg}')
logger::log_layout(logger_layout)





#####################################
## LOAD SOURCE SCRIPTS
#####################################

setwd(file.path(args$base_folder,"scripts"))
files.sources = list.files()
sapply(files.sources, source)
setwd(file.path(args$base_folder))


#######################################################
## Load the ENCODE Metadata per RBP
## If the file does not exists, download the data
#######################################################


metadata_path <- paste0(args$results_folder, "/metadata_",analysis_type,"_samples.tsv")

if ( !file.exists(metadata_path) ) {
  
  ## 1. Download metadata per RBP from ENCODE
  metadata <- ENCODEDownloadMetadata(experiment_type = analysis_type,
                                     results_path = args$results_folder,
                                     dependencies_path = args$dependencies_folder)
  
} else {
  
  metadata <- readr::read_tsv(file = metadata_path)
}



#######################################################
## Get the RBPs
## If it is an shRNA project, only focus on RBPs from
## Van Nostrand et al (2020)
#######################################################


if (analysis_type == "shRNA") {
  metadata <- metadata %>%
    dplyr::filter(if_any(c("Splicing.regulation", Spliceosome, "Novel.RBP", "Exon.Junction.Complex", NMD), ~ . != 0))
}

target_RBPs <- (metadata %>% dplyr::pull(target_gene) %>% unique())

# target_RBPs <- target_RBPs[1]
metadata_RBPs <- metadata %>% dplyr::filter(target_gene %in% target_RBPs)


#######################################################
## Download BAM files
#######################################################

download_bams = T

if (download_bams) {
  
  ## 2. Download WesternBlotting PCR efficiency
  DownloadKnockdownEfficiencyWB(metadata = metadata_RBPs, results.path = args$results_folder)
  
  ## 3. Download .bam files
  ENCODEDownloadBams(metadata_df = metadata_RBPs, results_path = args$results_folder,
                     regtools_path = "/home/grocamora/tools/regtools/build/",
                     samtools_path = "/home/grocamora/tools/samtools/bin/")
}

#####################################
## INIT - DATABASE PROJECT
#####################################

for (gtf_version in gtf_versions) {
  
  # gtf_version <- gtf_versions[1]
  
  #################################################
  ## SET UP VARIABLES
  
  database_base_folder <- paste0(args$base_folder, "/database/", main_project, "/")
  dir.create(database_base_folder, recursive = T, showWarnings = F)
  
  levelqc1_folder <- args$database_folder
  dir.create(levelqc1_folder, recursive = T, showWarnings = F)
  
  tpm_folder <- paste0(args$base_folder, "/results/", main_project, "/tpm/")
  dir.create(tpm_folder, recursive = T, showWarnings = F)

  
  #################################################
  ## PREPARE JUNCTIONS FROM KNOCKDOWN EXPERIMENTS

  logger::log_info(paste0(Sys.time(), "\t\t starting 'PrepareEncodeData' function..."))
  PrepareEncodeData(metadata = metadata_RBPs,
                    RBP.source.path = args$results_folder,
                    results.path = args$results_folder,
                    database.path = database_base_folder,
                    gtf.version = gtf_version,
                    blacklist.path = file.path(args$dependencies_folder, "hg38-blacklist.v2.bed"),
                    gtf.path = file.path(args$dependencies_folder, paste0("Homo_sapiens.GRCh38.", gtf_version, ".chr.gtf")),
                    ENCODE.silencing.series = analysis_type,
                    num.cores = 8)
  gc()

  #################################################
  ## JUNCTION PAIRING

  JunctionPairing(recount3.project.IDs = target_RBPs,
                  results.folder = args$results_folder,
                  replace = T,
                  num.cores = 10)


  GetAllAnnotatedSplitReads(recount3.project.IDs = target_RBPs,
                            database.folder = args$database_folder,
                            results.folder = args$results_folder,
                            num.cores = 10)


  GetAllRawJxnPairings(recount3.project.IDs = target_RBPs,
                           database.folder = args$database_folder,
                           results.folder = args$results_folder,
                           num.cores = 10)

  
  #################################################
  ## DATA BASE

  all_final_projects_used <- readRDS(file.path(args$results_folder, "all_final_projects_used.rds"))


  TidyDataPiorSQL(recount3.project.IDs = all_final_projects_used,
                  database.folder = args$database_folder,
                  levelqc1.folder = levelqc1_folder,
                  results.folder = args$results_folder)


  GenerateTranscriptBiotypePercentage(gtf.version = gtf_version,
                                      dependencies.folder = args$dependencies_folder,
                                      database.folder = args$database_folder,
                                      results.folder = args$results_folder)


  database_sqlite_file <- paste0(args$database_folder,  "/", main_project, ".sqlite")
  SqlDatabaseGeneration(database.path = database_sqlite_file,
                        recount3.project.IDs = all_final_projects_used,
                        remove.all = T,
                        database.folder = args$database_folder,
                        results.folder = args$results_folder,
                        dependencies.folder = args$dependencies_folder,
                        gtf.version = gtf_version,
                        discard.minor.introns = F)

  gc()
}
