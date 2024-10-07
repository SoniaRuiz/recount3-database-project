# source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/ENCODE_SR/ENCODE_Metadata_Extraction/Knockdown_Efficiency/TPM/KnockdownEfficiency_TPM_Script.R")

# 1. Load libraries and variables ----
## Required libraries ----
shhh <- suppressPackageStartupMessages
shhh(library(biomaRt))
shhh(library(doSNOW))


DownloadKnockdownEfficiencyTPM <- function(metadata,
                                           results.path) {
  
  
  main_path <- file.path(results.path, "/TPM/RBPs/")
  dir.create(main_path, recursive = T, showWarnings = F)
  
  
  ## Input Files ----

  metadata_path <- paste0(main_path, "metadata_samples.tsv")
  metadata_TPM_output <- paste0(main_path, "metadata_TPM_kEff.tsv")
  
  
  ## Define the algorithm variables ----
  download_only = F
  download_cores = 16
  overwrite_results = T
  
  # 2. Pipeline ----
  ## Generate the variables ----
  #metadata <- readr::read_delim(metadata_path, show_col_types = F) %>% as_tibble()
  
  target_RBPs <- metadata %>%
    #dplyr::filter(if_any(c(Splicing_regulation, Spliceosome, Exon_junction_complex, NMD), ~ . != 0)) %>%
    dplyr::filter(!is.na(gene_quantification_id)) %>%
    dplyr::pull(target_gene) %>%
    unique()
  
  print(target_RBPs)
  
  library("org.Hs.eg.db") # remember to install it if you don't have it already
  ensembl_target_RBPs <- mapIds(org.Hs.eg.db, keys = target_RBPs, keytype = "SYMBOL", column="ENSEMBL") %>% 
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    dplyr::rename("hgnc_symbol" = "rowname",
                  "ensembl_gene_id" = ".")
  
  metadata_filtered <- metadata %>%
    dplyr::filter(target_gene %in% target_RBPs, !is.na(gene_quantification_id)) %>%
    dplyr::select(target_gene, cell_line, experiment_type, sample_id, gene_quantification_id) %>%
    dplyr::mutate(path = paste0(main_path, target_gene, "/", experiment_type, "/")) %>%
    dplyr::left_join(y = ensembl_target_RBPs, 
                     by = c("target_gene" = "hgnc_symbol"), multiple = "all") %>%
    dplyr::relocate(ensembl_gene_id, .before = cell_line)
  
  ## Create the directories ----
  createDirectories(target_RBPs, metadata_filtered)
  
  ## Download the files and add a column with their path ----
  metadata_quantifications <- downloadGeneQuantifications(metadata_filtered,
                                                          download_cores, 
                                                          overwrite_results)
  print(metadata_quantifications)
  ## Extract & save the TPMs ----
  metadata_TPM <- extractTPM(metadata_quantifications,
                             output_file = main_path)
  
  ## Generate and save the knockdown efficiencies ----
  # metadata_kEff <- generateKnockdownEfficiency(metadata_TPM,
  #                                              output_file = metadata_TPM_output)
}




#############################################################
## HELPER FUNCTIONS
#############################################################

#' Creates the required subdirectories
#'
#' @param target_RBPs List of target genes for which to create the
#'   subdirectories.
#' @param metadata_filtered Data.frame containing the information required to
#'   download the gene quantifications per sample.
#'
#' @return
#' @export
createDirectories <- function(target_RBPs, 
                              metadata_filtered){
  for(i in seq(length(target_RBPs))){
    target_RBP <- target_RBPs[i]
    RBP_metadata <- metadata_filtered %>% dplyr::filter(target_gene == target_RBP)
    for (row in seq(nrow(RBP_metadata))) {
      document.path <- RBP_metadata[row, "path", T]
      dir.create(document.path, recursive = T, showWarnings = F)
    }
  }
}

#' Download the gene quantifications files
#'
#' From the samples' metadata previously generated, it download the gene
#' quantification files per sample. This file contains all the gene counts (and
#' TPM) found in the particular sample.
#'
#' @param metadata_filtered Data.frame containing the information required to
#'   download the gene quantifications per sample.
#' @param download_cores (Optional) Number of cores to use to download. Defaults to 1.
#' @param overwrite_results (Optional) Whether to download already found files. Defaults to FALSE.
#' @param silent (Optional) Whether to print progress bar of the download process. Defaults to FALSE.
#'
#' @return Data.frame with the information of the downloaded files.
#' @export
downloadGeneQuantifications <- function(metadata_filtered,
                                        download_cores = 1,
                                        overwrite_results = F,
                                        silent = F){
  if(!silent){
    pb <- txtProgressBar(max=nrow(metadata_filtered), style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
  }else{
    opts <- list()
  }
  
  cl <- makeCluster(download_cores)
  registerDoSNOW(cl)
  metadata_gene_quant <- foreach(row_index = seq(nrow(metadata_filtered)), .export = "getDownloadLinkGeneQuantification", .options.snow=opts) %dopar%{
    
    # row_index <- 1
    row = metadata_filtered[row_index, ]
    
    ## Two possible download paths:
    download_link <- getDownloadLinkGeneQuantification(row)
    
    ## Where to save the file
    file_path <- paste0(row$path, row$gene_quantification_id, ".tsv")
    row$file_path <- file_path
    
    ## If overwrite results is set to TRUE or if the file does not exists or if
    ## the file is too small, try to download it from the different links. If
    ## none it success, remove the resulted file. We also modify the column
    ## "file_path" if the file was successfully created.
    if(overwrite_results || !file.exists(file_path) || file.info(file_path)$size < 10){
      tryCatch({
        download.file(download_link, file_path, method = "wget", extra = "--quiet --no-check-certificate")
      }, error = function(e){
        file.remove(file_path)
      })
    }
    
    if(!file.exists(file_path) || file.info(file_path)$size < 10){
      row$file_path <- NA
    }
    
    return(row)
  } %>% dplyr::bind_rows()
  if(!silent) close(pb)
  stopCluster(cl)
  
  return(metadata_gene_quant)
}

#' Generates the download link
#'
#' @param row Row of the data.frame with the gene quantification ID.
#'
#' @return Download link of the gene quantification file.
#' @export
getDownloadLinkGeneQuantification <- function(row){
  gene_quantification_id <- row$gene_quantification_id
  
  return(paste0("https://www.encodeproject.org/files/", gene_quantification_id, "/@@download/", gene_quantification_id, ".tsv"))
}

#' Converts HGNC gene names to ENSEMBL ID
#'
#' @param gene_list List of genes in HGCN nomenclature
#'
#' @return Data.frame containing two columns: HGCN gene names and ENSEMBL gene
#'   names. In some situations, two or more ENSEMBL names are associated with
#'   the same HGCN name.
#' @export
translateGenes <- function(gene_list){
  ensembl <- biomaRt::useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl")
  
  
  gene_df <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                            filter = "hgnc_symbol",
                            mart = ensembl,
                            values = gene_list)
  
  return(gene_df)
}

#' Extract TPM of a particular
#'
#' Given a Data.frame with the samples' metadata, it reads the downloaded gene
#' quantification files and extract the TPM of the target gene.
#'
#' @param metadata_quantifications Data.frame with the information of the
#'   downloaded files.
#'
#' @return Data.frame with the TPM of the target gene for each sample.
#' @export
extractTPM <- function(metadata_quantifications,
                       output_file){
  
  all_RBPs <- metadata_quantifications$target_gene %>% unique
  
  
  for(RBP in all_RBPs) {
    
    # RBP <- all_RBPs[1]
    
    experiment_types <- metadata_quantifications %>%
      filter(target_gene == RBP) %>%
      pull(experiment_type) %>%
      unique
    
    message(Sys.time(), " - ", RBP, "...")
    
    for (type in experiment_types) {
      
      # type <- experiment_types[1]
      
      metadata_quantifications_local <- metadata_quantifications %>%
        filter(target_gene == RBP,
               experiment_type == type)
      
      metadata_TPM <- foreach(row_index = seq(nrow(metadata_quantifications_local))) %do%{
        suppressPackageStartupMessages(library(tidyverse))
        
        # row_index <- 1
        
        row = metadata_quantifications_local[row_index, ]
        
        message(Sys.time(), " - ", type, "...")
        
        file_path <- row$file_path
        
        read_tpm <- readr::read_delim(file_path, delim = "\t", show_col_types = F) %>%
          dplyr::mutate(sample_id = row$sample_id,
                        gene_quantification_id = row$gene_quantification_id) %>%
          dplyr::select(gene_id, TPM, sample_id) %>%
          filter(TPM > 0)
        
        return(read_tpm)
      } %>% dplyr::bind_rows()
      
      
      metadata_TPM_tidy <- metadata_TPM %>%
        mutate(gene_id = str_extract(gene_id, "[^.]+")) %>%
        group_by(sample_id, gene_id) %>%
        mutate(avg_TPM = TPM %>% mean) %>%
        ungroup() %>%
        dplyr::select(gene_id, sample_id, avg_TPM) %>%
        group_by(sample_id) %>%
        dplyr::distinct(gene_id, .keep_all = T) %>%
        ungroup() %>%
        pivot_wider(id_cols = gene_id, names_from = sample_id, values_from = avg_TPM )
      
      
      metadata_TPM_tidy[is.na(metadata_TPM_tidy)] <- 0
      
      metadata_TPM_tidy %>% nrow() %>% print()
      
      dir_path <- paste0(output_file, RBP, "/tpm/")
      dir.create(path = dir_path, recursive = T)
      
      saveRDS(object = metadata_TPM_tidy,
              file = file.path(dir_path, paste0(RBP, "_", type, "_tpm.rds")))
      
    }
  }
  
}

#' Calculates the knockdown efficiency (kEff)
#'
#' From the metadata data.frame with the TPM information, it groups by target
#' gene and by experiment type (control/case) and calculates the knockdown
#' efficiency first by calculating the averages in both clusters, and then by
#' applying the following equation:
#'
#' $$ kEff=(1-\frac{TPM_case}{TPM_control})*100% $$
#'
#' The same procedure is also executed but grouped by cell line, to also report
#' the efficiency for different cell lines.
#'
#' @param metadata_TPM Data.frame with the TPM of the target gene for each
#'   sample.
#' @param output_file (Optional) Path to where the resulted data.frame will be
#'   stored
#'
#' @return Data.frame with the target gene and the estimated knockdown
#'   efficiency. The efficiency is also calculated by cell line.
#' @export
generateKnockdownEfficiency <- function(metadata_TPM,
                                        output_file = ""){
  kEff_global <- metadata_TPM %>% 
    dplyr::group_by(target_gene, experiment_type) %>%
    dplyr::summarize(TPM_avg = mean(TPM, na.rm = T)) %>%
    tidyr::pivot_wider(id_cols = target_gene, names_from = experiment_type, values_from = TPM_avg) %>%
    dplyr::mutate(kEff = (1 - case/control)*100)
  
  kEff_cell_line <- metadata_TPM %>% 
    dplyr::group_by(target_gene, cell_line, experiment_type) %>%
    dplyr::summarize(TPM_avg = mean(TPM, na.rm = T)) %>%
    tidyr::pivot_wider(id_cols = c(target_gene, cell_line), names_from = experiment_type, values_from = TPM_avg) %>%
    dplyr::mutate(kEff = (1 - case/control)*100) %>%
    tidyr::pivot_wider(id_cols = target_gene, names_from = cell_line, values_from = kEff)
  
  metadata_kEff <- kEff_global %>% 
    dplyr::left_join(kEff_cell_line, by = "target_gene") %>%
    `colnames<-`(c("target_gene", "Avg_TPM_case", "Avg_TPM_control", "kEff", "kEff_HepG2", "kEff_K562"))
  
  if(output_file != ""){
    write.table(metadata_kEff, output_file, sep = "\t", row.names = F, quote = FALSE)
  }
  
  return(metadata_kEff)
}
