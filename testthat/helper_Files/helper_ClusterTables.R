## Libraries
library(foreach)
library(doParallel)
suppressWarnings(suppressMessages(library(GenomicRanges)))

## Load the necessary tables from the database
con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
df_novel <- dplyr::tbl(con, "novel") %>% dplyr::collect()
df_intron <- dplyr::tbl(con, "intron") %>% dplyr::collect()
df_transcript <- dplyr::tbl(con, "transcript") %>% dplyr::collect()
df_gene <- dplyr::tbl(con, "gene") %>% dplyr::collect()

if(identical(test_clusters, "all")) test_clusters = getClusters(con)

# test_clusters <- getClusters(con)
# test_clusters <- c("Adipose - Subcutaneous_ADIPOSE_TISSUE", "Brain - Frontal Cortex (BA9)_BRAIN")
# test_clusters <- c("Bladder_BLADDER", "Uterus_UTERUS")


#################### Functions

#' Reads the original split reads data per cluster
#'
#' @param split_read_counts_path Path to split read counts.
#'
#' @return Dataframe with the split read counts for the cluster.
#' @export
getClusterSplitReads <- function(split_read_counts_path){
  if(!exists("split_read_counts")){
    split_read_counts <- readRDS(split_read_counts_path)
    if(is.null(names(split_read_counts))){
      split_read_counts <- split_read_counts %>% tibble::as_tibble(rownames = "junID")
    }
    
    return(split_read_counts)
  }else{
    return(split_read_counts)
  }
}

#' Reads the original annotated intron details data per cluster
#'
#' @param annotated_SR_details_path Path to annotated intron details files.
#'
#' @return Dataframe with the annotated intron details for the cluster.
#' @export
getAnnotatedSR <- function(annotated_SR_details_path){
  if(!exists("annotated_SR_details")){
    annotated_SR_details <- readRDS(annotated_SR_details_path)
    return(annotated_SR_details)
  }else{
    return(annotated_SR_details)
  }
}

#' Reads the original gene TPM data per cluster
#'
#' @param tpm_path Path to Path to the gene TPM files.
#'
#' @return Dataframe with the gene TPM for the cluster.
#' @export
getClusterTPM <- function(tpm_path){
  tpm <- readRDS(tpm_path)
  return(tpm)
}

#' Fix the * strand in annotated_SR_details and split_read_counts
#'
#' The function modifies the two input dataframes in the global enviroment, it
#' doesn't return any of them.
#'
#' @param split_read_counts Dataframe containing the original information for
#'   the split read counts.
#' @param annotated_SR_details Dataframe containing the original annotated
#'   intron details.
#'
#' @return
#' @export
fixUnknownStrand <- function(split_read_counts, annotated_SR_details){
  if(any(stringr::str_detect(string = split_read_counts$junID, pattern = "\\*"))){
    if(!any(stringr::str_detect(string = annotated_SR_details$junID, pattern = "\\*"))){
      rm(annotated_SR_details)
      annotated_SR_details <<- getAnnotatedSR(annotated_SR_details_path)
    }
    
    split_read_counts <<- split_read_counts %>% 
      dplyr::left_join(annotated_SR_details %>% 
                         dplyr::select(junID, seqnames, start, end, strand),
                       by = "junID") %>%
      dplyr::select(-junID) %>%
      dplyr::mutate(junID = paste0("chr", seqnames, ":", start, "-", end, ":", strand)) %>%
      dplyr::select(-seqnames, -start, -end, -strand) %>%
      dplyr::relocate(junID)
  }
  
  if(any(stringr::str_detect(string = annotated_SR_details$junID, pattern = "\\*"))){
    annotated_SR_details <<- annotated_SR_details %>%
      dplyr::mutate(junID = paste0("chr", seqnames, ":", start, "-", end, ":", strand))
  }
}

#' Generates the distance dataframe for a particular novel junction type and
#' direction
#'
#' More information about the parameters in \link[GenomicRanges]{findOverlaps}.
#'
#' @param query GRanges object with junctions to find the overlaps.
#' @param subject GRanges object with junctions to find the overlaps.
#' @param type The type of overlap to look for.
#' @param sample The sample ID being studied.
#' @param junc_type The type of junction being studied.
#'
#' @return Dataframe with the novel junctions and the associated reference
#'   junctions for a given sample and type of novel donor and direction.
#' @export
getDistancesDataFrame <- function(query,
                                  subject,
                                  type,
                                  sample,
                                  junc_type) {
  ## Find the overlaps between the query and the subject GRanges
  overlaps <- GenomicRanges::findOverlaps(
    query = query,
    subject = subject,
    ignore.strand = FALSE,
    type = type
  )
  novel_junctions <- query[S4Vectors::queryHits(overlaps), ]
  ref_junctions <- subject[S4Vectors::subjectHits(overlaps), ]
  
  ## Calculation of the distances according to the direction of the junctions
  if (type == "start") {
    distance <- end(novel_junctions) - end(ref_junctions)
  } else {
    distance <- start(ref_junctions) - start(novel_junctions)
  }
  
  ## Dataframe with all the relevant information
  df <- tibble::tibble(
    sample = sample,
    type = junc_type,
    distance = distance,
    novel_junID = novel_junctions$junID,
    novel_counts = novel_junctions$counts,
    novel_seq = novel_junctions %>% seqnames() %>% as.character(),
    novel_start = novel_junctions %>% start(),
    novel_end = novel_junctions %>% end(),
    novel_strand = novel_junctions %>% strand() %>% as.character(),
    novel_width = novel_junctions %>% width(),
    ref_junID = ref_junctions$junID,
    ref_counts = ref_junctions$counts,
    ref_seq = ref_junctions %>% seqnames() %>% as.character(),
    ref_start = ref_junctions %>% start(),
    ref_end = ref_junctions %>% end(),
    ref_strand = ref_junctions %>% strand() %>% as.character(),
    ref_width = ref_junctions %>% width()
  )
  
  return(df)
}

#' Categorize the reference junctions in "novel", "donor", "acceptor" or "never"
#'
#' @param novel_donor_ratio Mis-splicing ratio of the novel donors.
#' @param novel_acceptor_ratio Mis-splicing ratio of the novel acceptor.
#'
#' @return Reference junction category.
#' @export
missplicingClass <- function(novel_donor_ratio, novel_acceptor_ratio) {
  ref_junction_category = ""
  if (novel_donor_ratio > 0 & novel_acceptor_ratio > 0) {
    ref_junction_category = "both"
  } else if (novel_donor_ratio > 0 & novel_acceptor_ratio == 0) {
    ref_junction_category = "donor"
  } else if (novel_donor_ratio == 0 & novel_acceptor_ratio > 0) {
    ref_junction_category = "acceptor"
  } else {
    ref_junction_category = "never"
  }
  
  return(ref_junction_category)
}
