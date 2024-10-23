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
library(doSNOW)
library(reticulate)
library("org.Hs.eg.db")

## .libPaths("~/R/x86_64-pc-linux-gnu-library/4.3")
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
    dependencies_folder = file.path("~/rds/hpc-work/SR/recount3-database-project/dependencies"),
    tools_folder = file.path("/home/sg2173/tools"),
    
    database_folder = file.path("~/rds/hpc-work/SR/recount3-database-project/database", main_project, gtf_version),
    database_root_folder = file.path("~/rds/hpc-work/SR/recount3-database-project/database", main_project),
    
    results_folder = file.path("~/rds/hpc-work/SR/recount3-database-project/results", main_project, gtf_version),
    logs_folder = file.path("~/rds/hpc-work/SR/recount3-database-project/logs"),
    tpm_folder = file.path("~/rds/hpc-work/SR/recount3-database-project/results", main_project, "tpm")
  )

dir.create(path = args$database_folder, recursive = T, showWarnings = F)
dir.create(path = args$database_root_folder, recursive = T, showWarnings = F)
dir.create(path = args$results_folder, recursive = T, showWarnings = F)
dir.create(path = args$logs, recursive = T, showWarnings = F)

dir.create(path = args$tpm_folder, recursive = T, showWarnings = F)
dir.create(path = args$tools_folder, recursive = T, showWarnings = F)



#####################################
## CREATE LOG
#####################################

log_file <- file.path(args$logs_folder, paste0("Splicing_Analysis_ENCODE_", gtf_version, "_", analysis_type, ".log"))
logger::log_appender(logger::appender_tee(log_file, append = T))
logger_layout <- logger::layout_glue_generator(format = '[{time}] [{level}] {msg}')
logger::log_layout(logger_layout)



#####################################
## LOAD SOURCE SCRIPTS
#####################################

files.sources = list.files(path = file.path(args$base_folder,"scripts"))
sapply(file.path(args$base_folder,"scripts", files.sources), source)


# logger::log_info("Source files loaded!")


#######################################################
## Load the ENCODE Metadata
## If the metadata file does not exists, then download
## the .bam files corresponding to each knockdown 
## experiment
#######################################################

metadata_path <- paste0(args$results_folder, "/metadata_", analysis_type, "_samples.tsv")

metadata <- if (!file.exists(metadata_path)) {
  
  ## 1. Download metadata per RBP from ENCODE
  metadata <- ENCODEDownloadMetadata(experiment_type = analysis_type,
                                     results_path = args$results_folder,
                                     dependencies_path = args$dependencies_folder)
  
  if (analysis_type == "shRNA") {
    metadata <- metadata %>%
      dplyr::filter(if_any(c("Splicing regulation", Spliceosome, "Novel RBP", "Exon Junction Complex", NMD), ~ . != 0))
  }
  
  ## 2. Download WesternBlotting PCR efficiency per RBP
  DownloadKnockdownEfficiencyWB(metadata, results.path = args$results_folder)
  
  ## 3. Download TPM gene data per RBP experiment
  DownloadKnockdownEfficiencyTPM(metadata, results.path = args$results_folder)
  
} else { readr::read_tsv(file = metadata_path, show_col_types = F) }

target_RBPs <- (metadata %>% dplyr::pull(target_gene) %>% unique())

###############################################################
target_RBPs <- "ADAR"
metadata <- metadata %>% dplyr::filter(target_gene == "ADAR")
###############################################################

## If not all sample experiments have been downloaded and their junctions extracted, download the bam files
if (nrow(CheckDownloadedFiles(RBP.metadata = metadata %>% filter(target_gene == target_RBPs[1]), 
                              RBP.path = file.path(args$results_folder, target_RBPs[1], "/"))) != 
    nrow(metadata %>% filter(target_gene == target_RBPs[1]))) {
  
  ## 4. Download .bam files
  ENCODEDownloadBams(metadata = metadata, 
                     results.path = args$results_folder,
                     regtools.path = file.path(args$tools_folder, "regtools/build/"),
                     samtools.path = file.path(args$tools_folder, "samtools-1.21/bin/"))
}



#####################################
## INIT - DATABASE PROJECT
#####################################

#for (gtf_version in gtf_versions) {
  
  # gtf_version <- gtf_versions[1]
  
  #################################################
  ## SET UP VARIABLES
  
  
  
  #################################################
  ## PREPARE JUNCTIONS FROM KNOCKDOWN EXPERIMENTS

  # logger::log_info("Starting 'PrepareEncodeData' function ...")
  PrepareEncodeData(metadata = metadata,
                    RBP.source.path = args$results_folder,
                    results.path = args$results_folder,
                    database.path = args$database_root_folder,
                    gtf.version = gtf_version,
                    blacklist.path = file.path(args$dependencies_folder, "hg38-blacklist.v2.bed"),
                    gtf.path = file.path(args$dependencies_folder, paste0("Homo_sapiens.GRCh38.", gtf_version, ".chr.gtf")),
                    ENCODE.silencing.series = analysis_type,
                    num.cores = 4)
  gc()

  
  
  
  # #################################################
  # ## JUNCTION PAIRING
  # 
  # JunctionPairing(recount3.project.IDs = target_RBPs,
  #                 results.folder = args$results_folder,
  #                 replace = T,
  #                 num.cores = 10)
  # 
  # 
  # GetAllAnnotatedSplitReads(recount3.project.IDs = target_RBPs,
  #                           database.folder = args$database_folder,
  #                           results.folder = args$results_folder,
  #                           num.cores = 10)
  # 
  # 
  # GetAllRawJxnPairings(recount3.project.IDs = target_RBPs,
  #                          database.folder = args$database_folder,
  #                          results.folder = args$results_folder,
  #                          num.cores = 10)
  # 
  # 
  # 
  # 
  # #################################################
  # ## DATA BASE
  # 
  # all_final_projects_used <- readRDS(file.path(args$results_folder, "all_final_projects_used.rds"))
  # 
  # 
  # TidyDataPiorSQL(recount3.project.IDs = all_final_projects_used,
  #                 database.folder = args$database_folder,
  #                 levelqc1.folder = args$database_folder,
  #                 results.folder = args$results_folder)
  # 
  # 
  # GenerateTranscriptBiotypePercentage(gtf.version = gtf_version,
  #                                     dependencies.folder = args$dependencies_folder,
  #                                     database.folder = args$database_folder,
  #                                     results.folder = args$results_folder)
  # 
  # 
  # database_sqlite_file <- paste0(args$database_folder,  "/", main_project, ".sqlite")
  # SqlDatabaseGeneration(database.path = database_sqlite_file,
  #                       recount3.project.IDs = all_final_projects_used,
  #                       remove.all = T,
  #                       database.folder = args$database_folder,
  #                       results.folder = args$results_folder,
  #                       dependencies.folder = args$dependencies_folder,
  #                       gtf.version = gtf_version,
  #                       discard.minor.introns = F)
  # 
  # gc()
#}
