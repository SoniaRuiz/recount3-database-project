library(tidyverse)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(DBI)
library(here)
library(logger)
library(foreach)
library(doParallel)
library(Biostrings)
library(protr)
library(optparse)
library(doSNOW)
library(reticulate)
library(org.Hs.eg.db)

## .libPaths("~/R/x86_64-pc-linux-gnu-library/4.3")
#####################################
## SET MAIN VARIABLES
#####################################

## This is the Ensembl gtf transcriptome version 
gtf_version <- c(111)

## This is the name of the project producing the database
min_supporting_reads <- 1
analysis_type = "shRNA" ## This can also be CRISPR

main_project <- paste0("ENCODE_SR_", min_supporting_reads, "read_", analysis_type)

args <-
  list(
    replace = T,
    num_cores = 1,
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

## 1. Download metadata per RBP from ENCODE
metadata <- if (!file.exists(metadata_path)) {
  metadata <- ENCODEDownloadMetadata(experiment_type = analysis_type,
                                     results_path = args$results_folder,
                                     dependencies_path = args$dependencies_folder)
} else { readr::read_tsv(file = metadata_path, show_col_types = F) }


target_RBPs <- (metadata %>% dplyr::pull(target_gene) %>% unique())

###############################################################
if (analysis_type == "shRNA") {
  metadata <- metadata %>%
    dplyr::filter(if_any(c("Splicing regulation", Spliceosome, "Novel RBP", "Exon Junction Complex", NMD), ~ . != 0))
}
target_RBPs <- c("EFTUD2","FUBP1","HNRNPC","NCBP2","PCBP2","PUF60",
                 "RAVER1","RBM22","U2AF1","U2AF2","PRPF4","SART3","SF3B4",
                 "MAGOH","SAFB2") %>% sort()
  # c("ADAR","AQR","TARDBP","UPF1","UPF2")
metadata <- metadata %>% dplyr::filter(target_gene %in% target_RBPs)
metadata %>% distinct(target_gene) %>% arrange(target_gene)
###############################################################


## 2. Download WesternBlotting PCR efficiency per RBP
# DownloadKnockdownEfficiencyWB(metadata, 
#                               results.path = args$results_folder,
#                               replace = args$replace,
#                               num.cores = args$num_cores)



## 3. Download TPM gene data per RBP experiment
DownloadKnockdownEfficiencyTPM(metadata, 
                               results.path = args$results_folder,
                               replace = F,
                               num.cores = args$num_cores)


## 4. If not all sample experiments have been downloaded and their junctions extracted, download the .bam files
for (RBP in target_RBPs) {
  if (nrow(CheckDownloadedFiles(RBP.metadata = metadata %>% filter(target_gene == RBP), 
                                RBP.path = file.path(args$results_folder, RBP, "/"))) != 
      nrow(metadata %>% filter(target_gene == RBP))) {
    
    ## Download .bam files
    ENCODEDownloadBams(metadata = metadata %>% filter(target_gene == RBP), 
                       results.path = args$results_folder,
                       regtools.path = file.path(args$tools_folder, "regtools/build/"),
                       samtools.path = file.path(args$tools_folder, "samtools-1.21/bin/"))
  }
  
}



#####################################
## INIT - DATABASE PROJECT
#####################################

# Loop to do annotation across multiple Ensembl versions
#for (gtf_version in gtf_versions) {
  
  # gtf_version <- gtf_versions[1]

  
  #################################################
  ## PREPARE JUNCTIONS FROM KNOCKDOWN EXPERIMENTS


  logger::log_info("Starting 'PrepareEncodeData' function ...")
  PrepareEncodeData(metadata = metadata,
                    RBP.source.path = args$results_folder,
                    results.path = args$results_folder,
                    database.path = args$database_root_folder,
                    gtf.version = gtf_version,
                    blacklist.path = file.path(args$dependencies_folder, "hg38-blacklist.v2.bed"),
                    gtf.path = file.path(args$dependencies_folder, paste0("Homo_sapiens.GRCh38.", gtf_version, ".chr.gtf")),
                    ENCODE.silencing.series = analysis_type,
                    num.cores = args$num_cores,
                    replace = args$replace)
  
  
  
  
  #################################################
  ## JUNCTION PAIRING AND QC
  

  JunctionPairing(recount3.project.IDs = target_RBPs,
                  results.folder = args$results_folder,
                  num.cores = args$num_cores,
                  replace = args$replace)


  GetAllAnnotatedSplitReads(recount3.project.IDs = target_RBPs,
                            database.folder = args$database_folder,
                            results.folder = args$results_folder,
                            num.cores = args$num_cores,
                            replace = args$replace)


  GetAllRawJxnPairings(recount3.project.IDs = target_RBPs,
                       database.folder = args$database_folder,
                       results.folder = args$results_folder,
                       num.cores = args$num_cores,
                       replace = args$replace)
  
  
  
  GetAllRawNovelCombos(recount3.project.IDs = target_RBPs,
                       database.folder = args$database_folder,
                       results.folder = args$results_folder,
                       replace = args$replace)


  
  #################################################
  ## DATA BASE PREP
  
   
  all_final_projects_used <- readRDS(file.path(args$results_folder, "all_final_projects_used.rds"))


  TidyDataPiorSQL(recount3.project.IDs = all_final_projects_used,
                  database.folder = args$database_folder,
                  levelqc1.folder = args$database_folder,
                  results.folder = args$results_folder,
                  replace = args$replace)


  GenerateTranscriptBiotypePercentage(gtf.version = gtf_version,
                                      dependencies.folder = args$dependencies_folder,
                                      database.folder = args$database_folder,
                                      replace = args$replace)
   
  
  
  #################################################
  ## DATA BASE BUILD
   
  database_sqlite_file <- paste0(args$database_folder,  "/", main_project, ".sqlite")
  SqlDatabaseGeneration(database.sqlite = database_sqlite_file,
                        recount3.project.IDs = all_final_projects_used,
                        database.folder = args$database_folder,
                        results.folder = args$results_folder,
                        dependencies.folder = args$dependencies_folder,
                        gtf.version = gtf_version,
                        max.ent.tool.path = paste0(args$dependencies_folder, "/fordownload/"),
                        bedtools.path = paste0(args$dependencies_folder, "/bedtools2/"),
                        hs.fasta.path = paste0(args$dependencies_folder, "/Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
                        phastcons.bw.path = paste0(args$dependencies_folder, "/hg38.phastCons17way.bw"),
                        cdts.bw.path = file.path(args$dependencies_folder, "CDTS_percentile_N7794_unrelated_all_chrs.bw"),
                        remove.all = T,
                        discard.minor.introns = F)

  
  
#}
