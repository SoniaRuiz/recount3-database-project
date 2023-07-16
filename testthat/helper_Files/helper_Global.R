library(DBI)
library(testthat)
library(GenomicRanges)
library(tidyverse, warn.conflicts = FALSE)

#################### Main testing parameters
gtf_version <- 105
main_project <- "SRP100948"
## Change this path accordingly. It should point to the 'splicing-accuracy-manuscript' root folder
main_path <- normalizePath(path = "./")

database_folder <-paste0(main_path, "/database/", main_project, "/", gtf_version, "/")
results_folder <- paste0(main_path, "/results/",  main_project, "/", gtf_version, "/")

database_path <- paste0(database_folder, "/", main_project, ".sqlite")

dependencies_folder <- paste0("/home/sruiz//PROJECTS/splicing-accuracy-manuscript/dependencies/")

#################### Control what is being tested
### Specific tests
test_DatabaseTables <- TRUE
test_DataTypes <- TRUE
test_IntronTable <- TRUE
test_NovelTable <- TRUE
test_GeneTable <- TRUE
test_TranscriptTable <- TRUE
test_ClusterTables <- TRUE

### Specific sub-tests
test_clinvar <- TRUE
test_conservation_CDTS <- TRUE
test_MANE <- TRUE
test_MaxEntScan <- TRUE
test_TSL <- TRUE
test_biotypes <- TRUE
test_cluster_data <- T

### Clusters to test
#test_clusters <- "all" # Input a list of cluster such as c("Adipose - Subcutaneous_ADIPOSE_TISSUE", "Brain - Frontal Cortex (BA9)_BRAIN")
test_clusters <- c(paste0("AD_",main_project), paste0("control_", main_project))
num_cores <- 8

#################### Additional files and paths
## Additional Files
reference_transcriptome_path <- path.expand(paste0(dependencies_folder, "/Homo_sapiens.GRCh38.105.chr.gtf"))
clinvar_path <- path.expand(paste0(dependencies_folder, "/clinvar_intronic_tidy.rda"))
CNC_CDTS_CONS_gr_path <- path.expand(paste0( dependencies_folder, "/CNC_CDTS_CONS_gr.rda"))
hg_mane_transcripts_path <- path.expand(paste0(dependencies_folder, "/MANE.GRCh38.v1.0.ensembl_genomic.gtf"))

# df_all_introns_introverse_tidy_path <- path.expand(paste0(main_path, "/database/v", gtf_version, "/", main_project, "/df_all_introns_database_tidy.rds"))
# all_annotated_SR_details_length_105_raw_biotype_path <-paste0(main_path, "/results/all_split_reads_all_tissues_PC_biotype.rds")

df_all_introns_introverse_tidy_path <- path.expand(paste0(database_folder, "/all_split_reads_qc_level2.rds"))
all_annotated_SR_details_length_105_raw_biotype_path <-paste0(results_folder, "/all_split_reads_qc_level2_PC_biotype.rds")

## MaxEntScan requirements
bedtools_path <- path.expand(path = "/tools/bedtools2/")
fasta_path <- path.expand(paste0(dependencies_folder, "/Homo_sapiens.GRCh38.dna.primary_assembly.fa"))
fordownload_path <- path.expand(paste0(dependencies_folder, "/fordownload/"))

#################### Variables to validate the data
### Database tables:
expected_main_tables = c("intron", "novel", "gene", "master", "transcript") %>% sort()
### Intron table:
valid_intron_strands = c("+", "-") %>% sort()
valid_intron_seqnames = c(paste0("chr", seq(1, 22)), "chrX", "chrY")
valid_intron_clinvar = c("acceptor", "donor", "") %>% sort()
valid_intron_TSL <- c(1, 2, 3, 4, 5, 10)
### Novel table:
valid_novel_types <- c("novel_donor", "novel_acceptor") %>% sort()
valid_novel_strands = c("+", "-") %>% sort()
valid_novel_seqnames <- c(paste0("chr", seq(1, 22)), "chrX", "chrY")

#################### Functions

#' Generates path to original split reads data per cluster
#'
#' @param table_name Name of the table from which to read the cluster and
#'   project.
#' @param projects_path Path to the projects' folder.
#' @param gtf_version Version of the reference transcriptome.
#' @param main_project Name of the main project to study. Only tested
#'   "introverse" and "splicing"
#'
#' @return Path to split read counts.
#' @export
generateClusterSplitReadsPath <- function(table_name, results_folder, gtf_version, main_project){
  full_name <- stringr::str_split_fixed(table_name, "_", n = 2)
  project <- full_name[2]
  cluster <- full_name[1]
  
  cluster_split_reads_path <- paste0(results_folder, "/", project,  "/base_data/", project, "_", cluster, "_split_read_counts.rds")
  return(cluster_split_reads_path)
}

#' Generates path to annotated intron details data per cluster
#'
#' @param table_name Name of the table from which to read the cluster and
#'   project.
#' @param projects_path Path to the projects' folder.
#' @param gtf_version Version of the reference transcriptome.
#' @param main_project Name of the main project to study. Only tested
#'   "introverse" and "splicing"
#'
#' @return Path to annotated intron details files.
#' @export
generateAnnotatedSRdetailsPath <- function(table_name, results_folder, gtf_version, main_project){
  full_name <- stringr::str_split_fixed(table_name, "_", n = 2)
  project <- full_name[2]
  cluster <- full_name[1]
  
  cluster_split_reads_path <- paste0(results_folder, "/", project,  "/base_data/", 
                                     project, "_", cluster, "_all_split_reads.rds")
  return(cluster_split_reads_path)
}

#' Generates path to gene TPM data per cluster
#'
#' @param table_name Name of the table from which to read the cluster and
#'   project.
#' @param projects_path Path to the projects' folder.
#' @param gtf_version Version of the reference transcriptome.
#' @param main_project Name of the main project to study. Only tested
#'   "introverse" and "splicing"
#'
#' @return Path to gene TPM files.
#' @export
generateTPMpath <- function(table_name, results_folder, gtf_version, main_project){
  full_name <- stringr::str_split_fixed(table_name, "_", n = 2)
  project <- full_name[2]
  cluster <- full_name[1]
  
  tpm_path <- paste0(results_folder, "/", project, "/tpm/", project, "_", cluster, "_tpm.rds")
  return(tpm_path)
}

#' Get the cluster tables from the database
#'
#' @param con DBIConnection object used to direct commands to the database
#'   engine.
#'
#' @return List of cluster tables in the database.
#' @export
getTissueTables <- function(con){
  tables <- DBI::dbListTables(con)
  tissue_tables <- tables[grepl("misspliced", tables)] %>% sort()
  return(tissue_tables)
}

#' Removes special characters from a text
#'
#' @param text Text to process. 
#'
#' @return Modified text.
#' @export
clearSpecialCharacters <- function(text){
  return(gsub(" ", "", gsub("-", "_", text)))
}

#' Remove all variables from the global environment
#'
#' @param con DBIConnection object used to direct commands to the database
#'   engine.
#'
#' @return
#' @export
clearAllVariables <- function(con = NULL){
  if(!is.null(con)) DBI::dbDisconnect(con)
  rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)
  invisible(gc())
}

#' Get the cluster names from the database
#'
#' @param con DBIConnection object used to direct commands to the database
#'   engine.
#'
#' @return List of clusters in the database.
#' @export
getClusters <- function(con){
  tissue_tables <- getTissueTables(con)
  tissue_tables <- gsub("_misspliced", "", tissue_tables)
  tissue_tables <- gsub("_nevermisspliced", "", tissue_tables)
  clusters <- tissue_tables %>% unique
  
  return(clusters)
}

#' Split junction coordinates to dataframe
#'
#' Given a junction identifier and coordinates in the form of locus (i.e.
#' seqname:start-end:strand), it splits each element and generates a dataframe
#' with the different fields.
#'
#' @param junID Junction identifier.
#' @param jun_coordinates Junction coordinates in the form of locus (i.e.
#'   seqname:start-end:strand).
#' @param fix_data_types whether to set the start and end positions as
#'   numerical.
#' @param remove_chr whether to remove the characters "chr" from the seqnames.
#'
#' @return dataframe containing the extracted information from the locus.
#' @export
getGRdata <- function(junID,
                      jun_coordinates,
                      fix_data_types = T,
                      remove_chr = T){
  jun_coordinates_split <- stringr::str_split_fixed(jun_coordinates, ":", n = 3)
  jun_seqnames <- jun_coordinates_split[, 1]
  jun_strands <- jun_coordinates_split[, 3]
  
  jun_positions <- jun_coordinates_split[, 2]
  jun_positions_split <- stringr::str_split_fixed(jun_positions, "-", n = 2)
  jun_start = jun_positions_split[, 1]
  jun_end = jun_positions_split[, 2]
  
  if(remove_chr){
    jun_seqnames <- gsub("chr", "", jun_seqnames)
  }
  GRdata <- tibble::tibble(junID = junID,
                           seqnames = jun_seqnames,
                           start = jun_start,
                           end = jun_end,
                           strand = jun_strands)
  
  if(fix_data_types) {
    GRdata <- GRdata %>%
      dplyr::mutate(start = as.numeric(start),
                    end = as.numeric(end))
  }
  
  return(GRdata)
}

#' Executes the MaxEntScan
#'
#' @param GRdata GRanges class object with the relevant junctions.
#' @param bedtools_path path to the
#'   \href{https://bedtools.readthedocs.io/en/latest/}{bedtools} executable. Can
#'   be left empty if bedtools is in default PATH.
#' @param fasta_path path to the fasta .fa file for the reference genome.
#' @param fordownload_path path to the MaxEntScan pearl scripts. Can be
#'   downloaded from
#'   \href{http://hollywood.mit.edu/burgelab/maxent/download/}{here}.
#'
#' @return junction dataframe with the ss5 and ss3 scores.
#' @export
generateMaxEntScore <- function(GRdata,
                                bedtools_path,
                                fasta_path,
                                fordownload_path) {
  current_wd <- getwd()
  
  ## Add extra positions to calculate the MaxEntScore
  GRdata <- GRdata %>% dplyr::mutate(
    seqnames = seqnames %>% as.character(),
    donorSeqStart = ifelse(strand == "-", end - 6, start - 4),
    donorSeqStop = ifelse(strand == "-", end + 3, start + 5),
    acceptorSeqStart = ifelse(strand == "-", start - 4, end - 20),
    acceptorSeqStop = ifelse(strand == "-", start + 19, end + 3)
  )
  
  ## Generate the dataframes for both donor and acceptors.
  donor_df <- GRdata %>%
    dplyr::select(seqnames,
                  starts = donorSeqStart,
                  end = donorSeqStop,
                  names = junID,
                  strands = strand
    ) %>%
    dplyr::mutate(scores = ".", .before = strands)
  
  acceptor_df <- GRdata %>%
    dplyr::select(seqnames,
                  starts = acceptorSeqStart,
                  end = acceptorSeqStop,
                  names = junID,
                  strands = strand
    ) %>%
    dplyr::mutate(scores = ".", .before = strands)
  
  ## Generate temporary files to store the sequences
  tmp_file <- tempfile()
  tmp_file_seq <- tempfile()
  
  ############ MaxEntScore ############
  ##
  ## 1. For both donor and acceptor dataframes, extract the sequences of the
  ## junctions and add them to the main dataframe.
  ## 2. Calculate the ss5 and ss3 scores (donor and acceptor sequences).
  ## 3. Add the scores to the main dataframe.
  
  ## Donor sequence
  data.table::fwrite(donor_df, file = tmp_file, quote = F, sep = "\t", row.names = F, col.names = F, scipen = 50)
  rm(donor_df)
  system(paste0("bedtools getfasta -name -s -fi ", fasta_path, 
                " -bed ", tmp_file, 
                " -tab -fo ", tmp_file_seq))
  
  donor_sequences_input <- read.delim(tmp_file_seq, header = F)
  stopifnot(identical(
    gsub("\\(\\+\\)", "", gsub("\\(\\*\\)", "", gsub("\\(-\\)", "", as.character(donor_sequences_input$V1)))),
    GRdata$junID %>% as.character()
  ))
  
  ## Acceptor sequence
  data.table::fwrite(acceptor_df, file = tmp_file, quote = F, sep = "\t", row.names = F, col.names = F, scipen = 50)
  rm(acceptor_df)
  system(paste0("bedtools getfasta -name -s -fi ", fasta_path, 
                " -bed ", tmp_file, 
                " -tab -fo ", tmp_file_seq))
  
  acceptor_sequences_input <- read.delim(tmp_file_seq, header = F)
  stopifnot(identical(
    gsub("\\(\\+\\)", "", gsub("\\(\\*\\)", "", gsub("\\(-\\)", "", as.character(acceptor_sequences_input$V1)))),
    GRdata$junID %>% as.character()
  ))
  
  ## Add to dataframe
  GRdata <- GRdata %>%
    dplyr::mutate(donor_sequence = as.character(donor_sequences_input$V2)) %>%
    dplyr::mutate(acceptor_sequence = as.character(acceptor_sequences_input$V2))
  
  ## Generate MaxEntScore
  tmp_file <- tempfile()
  setwd(fordownload_path)
  
  data.table::fwrite(list(GRdata$donor_sequence), file = tmp_file, quote = F, row.names = F, col.names = F, scipen = 50)
  ss5score_vector <- read.delim(pipe(paste0("perl ", fordownload_path, "score5.pl ", tmp_file)), header = F)
  GRdata <- GRdata %>% dplyr::mutate(ss5score = ss5score_vector$V2) 
  rm(ss5score_vector)
  
  data.table::fwrite(list(GRdata$acceptor_sequence), file = tmp_file, quote = F, row.names = F, col.names = F, scipen = 50)
  ss3score_vector <- read.delim(pipe(paste0("perl ", fordownload_path, "score3.pl ", tmp_file)), header = F)
  GRdata <- GRdata %>% dplyr::mutate(ss3score = ss3score_vector$V2)
  rm(ss3score_vector)
  
  setwd(current_wd)
  rm(tmp_file, tmp_file_seq, donor_sequences_input, acceptor_sequences_input)
  invisible(gc())
  return(GRdata)
}

