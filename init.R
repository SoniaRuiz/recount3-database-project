#####################################
## init.R file
## Contains the functions to download, qc, pair
## and database the exon-exon split read data 
## from a RNA-sequencing data project from recount3
#####################################

library(tidyverse)
library(GenomicRanges)
library(DBI)

# source("/home/sg2173/PROYECTS/SR/recount3-database-project/init.R")

# base_folder <- "/mnt/PROJECTS/splicing-accuracy-manuscript/"

#####################################
## SET MAIN VARIABLES
#####################################

## This is the Ensembl gtf transcriptome version 
gtf_version <- c(111)

## This is the name of the project producing the database
supportive_reads <- 1
data_subsample = F
main_project_identifier <- "SRP058181" 
project_name <- paste0(main_project_identifier, "_", supportive_reads, "read_subsample", data_subsample)


args <-
  list(
    replace = F,
    num_cores = 4,
    base_folder = here::here(),
    dependencies_folder = file.path("~/rds/hpc-work/SR/recount3-database-project/dependencies"),
    #data_source = "data_sources/tcga" 
    #data_source = "data_sources/gtex"
    data_source = "data_sources/sra",
    recount3_project_IDs = main_project_identifier,
    gtf_version = gtf_version,
    database_folder = file.path("~/rds/hpc-work/SR/recount3-database-project/database", project_name, gtf_version),
    results_folder = file.path("~/rds/hpc-work/SR/recount3-database-project/results", project_name, gtf_version),
    tpm_folder = file.path("~/rds/hpc-work/SR/recount3-database-project/results", project_name, "tpm"),
    logs_folder = file.path("~/rds/hpc-work/SR/recount3-database-project/logs")
    
  )


dir.create(path = args$database_folder, recursive = T, showWarnings = F)
dir.create(path = args$results_folder, recursive = T, showWarnings = F)
dir.create(path = args$tpm_folder, recursive = T, showWarnings = F)
dir.create(path = args$logs_folder, recursive = T, showWarnings = F)


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

files.sources = list.files(path = file.path(args$base_folder,"scripts"))
sapply(file.path(args$base_folder,"scripts", files.sources), source)


#####################################
## INIT LOGS
#####################################

log_file <- here::here(file.path(args$logs_folder, paste0("Splicing_Database_Generation_",project_name,".log")))
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
                       results.folder = args$results_folder,
                       replace = args$replace)

  PrepareRecount3Data(recount3.project.IDs = args$recount3_project_IDs,
                      data.source = args$data_source,
                      results.folder = args$results_folder,
                      levelqc1.folder = args$database_folder,
                      supporting.reads = supportive_reads,
                      num.cores = args$num_cores,
                      replace = args$replace,
                      subsampling = data_subsample)

  
  
  
  #################################################
  ## JUNCTION PAIRING

  JunctionPairing(recount3.project.IDs = args$recount3_project_IDs,
                  results.folder = args$results_folder,
                  num.cores = args$num_cores,
                  replace = args$replace)

  
  GetAllAnnotatedSplitReads(recount3.project.IDs = args$recount3_project_IDs,
                            database.folder = args$database_folder,
                            results.folder = args$results_folder,
                            num.cores = args$num_cores,
                            replace = args$replace)

  GetAllRawJxnPairings(recount3.project.IDs = args$recount3_project_IDs,
                       database.folder = args$database_folder,
                       results.folder = args$results_folder,
                       num.cores = args$num_cores,
                       replace = args$replace)

  GetAllRawNovelCombos(recount3.project.IDs = args$recount3_project_IDs,
                       database.folder = args$database_folder,
                       results.folder = args$results_folder,
                       replace = args$replace)

  
  
  
  #################################################
  ## DATA BASE PREP

  recount3_project_IDs_used <- readRDS(file = file.path(args$results_folder, "all_final_projects_used.rds"))


  TidyDataPiorSQL(recount3.project.IDs = recount3_project_IDs_used,
                  database.folder = args$database_folder,
                  levelqc1.folder = args$database_folder,
                  results.folder = args$results_folder,
                  replace = args$replace)


  GenerateTranscriptBiotypePercentage(gtf.version = args$gtf_version,
                                      dependencies.folder = args$dependencies_folder,
                                      database.folder = args$database_folder,
                                      replace = args$replace)


  GenerateRecount3TPM(recount3.project.IDs = c(recount3_project_IDs_used),
                      data.source = args$data_source,
                      tpm.folder = args$tpm_folder,
                      results.folder = args$results_folder,
                      replace = args$replace)


  
  
  #################################################
  ## DATA BASE BUILD
  

  database_sqlite_path <- paste0(args$database_folder,  "/", project_name, ".sqlite")

  SqlDatabaseGeneration(database.sqlite = database_sqlite_path,
                        recount3.project.IDs = recount3_project_IDs_used,
                        database.folder = args$database_folder,
                        results.folder = args$results_folder,
                        dependencies.folder = args$dependencies_folder,
                        gtf.version = args$gtf_version,
                        max.ent.tool.path = paste0(args$dependencies_folder, "/fordownload/"),
                        bedtools.path = paste0(args$dependencies_folder, "/bedtools2/"),
                        hs.fasta.path = paste0(args$dependencies_folder, "/Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
                        phastcons.bw.path = paste0(args$dependencies_folder, "/hg38.phastCons17way.bw"),
                        cdts.bw.path = file.path(args$dependencies_folder, "CDTS_percentile_N7794_unrelated_all_chrs.bw"),
                        remove.all = F,
                        discard.minor.introns = T)
  
#}

