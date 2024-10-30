##############################################################
## CODE Adapted from:
## https://github.com/guillermo1996/ENCODE_Metadata_Extraction
##############################################################

library(biomaRt)


#' Title
#' Downloads from the ENCODE platform the WB details regarding the knockdown efficiency of each RBP
#' @param metadata 
#' @param results.path 
#'
#' @return
#' @export
#'
#' @examples
DownloadKnockdownEfficiencyWB <- function(metadata,
                                          results.path,
                                          replace,
                                          num.cores) {
  
  ## Create and load python environment -------------------------
  
  virtualenv <- here::here("virtual_env/ENCODE_python_env")
  reticulate::virtualenv_create(virtualenv)
  reticulate::use_virtualenv(virtualenv)
  
  reticulate::py_available()
  reticulate::py_install(envname = virtualenv, packages = "tesseract")
  reticulate::py_install(envname = virtualenv, packages = "pytesseract") 
  reticulate::py_install(envname = virtualenv, packages = "pymupdf")
  reticulate::py_install(envname = virtualenv, packages = "frontend")
  
  #reticulate::virtualenv_remove(envname = virtualenv, packages = c('tesseract'))
  
  
  pytesseract <- reticulate::import("pytesseract")
  PIL <- reticulate::import("PIL")
  
  ## Set Variables ----------------------------------------------

  results_path <- file.path(results.path, "/WesternBlotting_PCR/RBPs/")
  dir.create(results_path, recursive = T, showWarnings = F)
  metadata_WB_output <- file.path(results_path, "metadata_WB_kEff.tsv")
  metadata_PCR_output <- file.path(results_path, "metadata_PCR_kEff.tsv")
  
  #download_cores = 1
  resize_perc = 0.25
  min_width = 600
  
  target_RBPs <- metadata %>% dplyr::pull(target_gene) %>%  unique()
  
  metadata_filtered <- metadata %>%
    dplyr::select(target_gene, cell_line, experiment_id, biosample, bio_rep, document, biosample_alias) %>%
    dplyr::mutate(path = paste0(results_path, target_gene, "/", experiment_id, "/", bio_rep, "/"))
  
  
  # 2. Pipeline --------------------------------------------------
  
  ## Create the directories
  CreateDirectories(target_RBPs, metadata_filtered)
  
  ## Download the files and add a column with their path
  metadata_documents <- DownloadCharacterizationDocuments(metadata_filtered = metadata_filtered, 
                                                          download_cores = num.cores, 
                                                          overwrite_results = replace, 
                                                          silent = F)
  
  
  ## Check the existence of the files
  for (row_index in seq(nrow(metadata_documents))) {
    row = metadata_documents[row_index, ]
    file_path <- row$file_path
    path <- row$path
    
    if (!file.exists(file_path) || file.info(file_path)$size < 10) {
      logger::log_warn("Error for row ", row_index, "! File path ", path)
    }
  }
  
  ## Extract the images of all files
  metadata_images <- ExtractImages(metadata_documents, overwrite_results = replace, virtualenv_path = virtualenv)
  
  ## Extract text from images ----
  ## This cannot be separated into an external function.
  ## Probably because of some incompatibility with the reticulate library to use
  ## python.
  metadata_kEff <- foreach(row_index = seq(nrow(metadata_images))) %do%{
    
    # row_index<-1
    row = metadata_images[row_index, ]
    path <- row$path
    image_path <- row$image_path
    
    ## If the image does not exists, return the unmodified row
    if (!file.exists(image_path)) return(row)
    
    
    image <- PIL$Image$open(fp = fs::path_expand(image_path))
    image_cropped <- image$crop(list(0, image$height*0.9, 0.8*image$width, image$height))
    image_small <- ResizeImage(image, resize_perc, min_width)
    image_small$save(file.path(fs::path_expand(path), "cropped_image.png"))
    
    # pytesseract$pytesseract$tesseract_cmd <- "/home/sg2173/PROJECTS/SR/recount3-database-project/virtual_env/ENCODE_python_env/lib/python3.11/site-packages/tesseract"
    
    text_df <- pytesseract$image_to_string(image = image_small) %>%
      str_replace_all("\\f", "") %>%
      str_split("\\n", simplify = T) %>%
      str_split(" ", simplify = T) %>%
      .[1:2, ] %>%
      as_tibble()
    
    if (ncol(text_df) == 5) {
      
      kEff_df <- text_df %>% 
        `colnames<-`(c("method", "K562_1", "K562_2", "HepG2_1", "HepG2_2")) %>%
        dplyr::mutate(across(-method, function(x) as.numeric(sub("%", "", x)))) %>%
        dplyr::mutate(K562 = rowMeans(dplyr::select(., K562_1, K562_2)),
                      HepG2 = rowMeans(dplyr::select(., HepG2_1, HepG2_2))) %>%
        dplyr::select(method, K562, HepG2)
      
      row$WB_HepG2 <- kEff_df %>% dplyr::filter(method == "Western") %>% dplyr::pull(HepG2)
      row$WB_K562 <- kEff_df %>% dplyr::filter(method == "Western") %>% dplyr::pull(K562)
      
      row$PCR_HepG2 <- kEff_df %>% dplyr::filter(method != "Western") %>% dplyr::pull(HepG2)
      row$PCR_K562 <- kEff_df %>% dplyr::filter(method != "Western") %>% dplyr::pull(K562)
      
    } else if (ncol(text_df) == 3) {
      
      kEff_df <- text_df %>% 
        `colnames<-`(c("method", "cell_line_1", "cell_line_2")) %>%
        dplyr::mutate(across(-method, function(x) as.numeric(sub("%", "", x)))) %>%
        dplyr::mutate(cell_line = rowMeans(dplyr::select(., cell_line_1, cell_line_2))) %>%
        dplyr::select(method, cell_line)
      
      row$WB_HepG2 <- ifelse(row$cell_line == "HepG2", kEff_df %>% dplyr::filter(method == "Western") %>% dplyr::pull(cell_line), NA)
      row$WB_K562 <- ifelse(row$cell_line == "K562", kEff_df %>% dplyr::filter(method == "Western") %>% dplyr::pull(cell_line), NA)
      
      row$PCR_HepG2 <- ifelse(row$cell_line == "HepG2", kEff_df %>% dplyr::filter(method != "Western") %>% dplyr::pull(cell_line), NA)
      row$PCR_K562 <- ifelse(row$cell_line == "K562", kEff_df %>% dplyr::filter(method != "Western") %>% dplyr::pull(cell_line), NA)
      
    } else {
      
      logger::WARN("Error in row ", row_index, ". Columns are not valid")
      row$WB_HepG2 <- NA
      row$WB_K562 <- NA
      row$PCR_HepG2 <- NA
      row$PCR_K562 <- NA
    }
    
    return(row)
  } %>% dplyr::bind_rows()
  
  ## Test for consistency between the cell lines ----
  for (target_RBP in target_RBPs) {
    metadata_RBP <- metadata_kEff %>%  filter(target_gene == target_RBP)
    WB_HepG2 <- metadata_RBP$WB_HepG2
    WB_K562 <- metadata_RBP$WB_K562
    PCR_HepG2 <- metadata_RBP$PCR_HepG2
    PCR_K562 <- metadata_RBP$PCR_K562
    
    if(length(unique(na.omit(WB_HepG2))) > 1) logger::ERROR("Error in RBP ", target_RBP, " cell line HepG2, method WB.")
    if(length(unique(na.omit(WB_K562))) > 1) logger::ERROR("Error in RBP ", target_RBP, " cell line K562, method WB.")
    
    if(length(unique(na.omit(PCR_HepG2))) > 1) logger::ERROR("Error in RBP ", target_RBP, " cell line HepG2, method PCR.")
    if(length(unique(na.omit(PCR_K562))) > 1) logger::ERROR("Error in RBP ", target_RBP, " cell line K562, method PCR.")
  }
  
  ## Write the knockdown efficiency table to disk ----
  WriteEfficiencyTable(metadata_kEff, "WB", metadata_WB_output)
  WriteEfficiencyTable(metadata_kEff, "PCR", metadata_PCR_output)
}





#' Title
#' Downloads the gene TPM values from each ENCODE shRNA RBP knockdown experiment studied
#' @param metadata 
#' @param results.path 
#'
#' @return
#' @export
#'
#' @examples
DownloadKnockdownEfficiencyTPM <- function(metadata,
                                           results.path,
                                           replace,
                                           num.cores) {
  
  if (replace) {
    ## Set Variables ----------------------------------------------
    
    results_path <- results.path # file.path(results.path, "/TPM/RBPs/")
    dir.create(results_path, recursive = T, showWarnings = F)
    
    main_path <- results.path # file.path(results.path, "/TPM/RBPs/")
    dir.create(main_path, recursive = T, showWarnings = F)
    
    
    ## Input Files ----
    
    metadata_path <- file.path(main_path, "metadata_samples.tsv")
    metadata_TPM_output <- file.path(main_path, "metadata_TPM_kEff.tsv")
    
    
    ## Define the algorithm variables ----
    download_only = F
    
    # 2. Pipeline ----
    ## Generate the variables ----
    #metadata <- readr::read_delim(metadata_path, show_col_types = F) %>% as_tibble()
    
    target_RBPs <- metadata %>%
      dplyr::filter(!is.na(gene_quantification_id)) %>%
      dplyr::pull(target_gene) %>%
      unique()
    
    #print(target_RBPs)
    
    ensembl_target_RBPs <- mapIds(org.Hs.eg.db, keys = target_RBPs, keytype = "SYMBOL", column="ENSEMBL") %>% 
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      dplyr::rename("hgnc_symbol" = "rowname", "ensembl_gene_id" = ".")
    
    
    
    metadata_filtered <- metadata %>%
      dplyr::filter(target_gene %in% target_RBPs, !is.na(gene_quantification_id)) %>%
      dplyr::select(target_gene, cell_line, experiment_type, sample_id, gene_quantification_id) %>%
      dplyr::mutate(path = paste0(main_path, target_gene, "/", experiment_type, "/")) %>%
      dplyr::left_join(y = ensembl_target_RBPs, 
                       by = c("target_gene" = "hgnc_symbol"), multiple = "all") %>%
      dplyr::relocate(ensembl_gene_id, .before = cell_line)
    
    ## Create the directories ----
    CreateDirectories(target_RBPs, metadata_filtered)
    
    ## Download the files and add a column with their path ----
    metadata_quantifications <- downloadGeneQuantifications(metadata_filtered,
                                                            download_cores = num.cores, 
                                                            replace )
    print(metadata_quantifications)
    
    ## Extract & save the TPMs ----
    ExtractTPM(metadata_quantifications, output_file = main_path)
    
    ## Generate and save the knockdown efficiencies ----
    # metadata_kEff <- generateKnockdownEfficiency(metadata_TPM,
    #                                              output_file = metadata_TPM_output)
  }
  
}




## HELPER FUNCTIONS - WESTERN BLOT ---------------------------------------------


#' Creates the required subdirectories
#'
#' @param target_RBPs List of target genes for which to create the
#'   subdirectories.
#' @param metadata_filtered Data.frame containing the information required to
#'   download the biosample characterization document per sample.
#'
#' @return
#' @export
CreateDirectories <- function(target_RBPs, 
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


#' Download the biosample preparation and characterization documents
#'
#' From the samples' metadata previously generated, it download the biosample
#' preparation and characterization documents per sample. This file contains the
#' reported efficiencies in both qRT-PCR and Western Blot.
#'
#' @param metadata_filtered Data.frame containing the information required to
#'   download the biosample preparation and characterization documents per sample.
#' @param download_cores (Optional) Number of cores to use to download. Defaults
#'   to 1.
#' @param overwrite_results (Optional) Whether to download already found
#'   files. Defaults to FALSE.
#' @param silent (Optional) Whether to print progress bar of the download
#'   process. Defaults to FALSE.
#'
#' @return Data.frame with the information of the downloaded files.
#' @export
DownloadCharacterizationDocuments <- function(metadata_filtered,
                                              download_cores,
                                              overwrite_results,
                                              silent){
  
  if (!silent) {
    pb <- txtProgressBar(max=nrow(metadata_filtered), style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
  } else {
    opts <- list()
  }
  
  cl <- makeCluster(download_cores)
  registerDoSNOW(cl)
  metadata_documents <- foreach(row_index = seq(nrow(metadata_filtered)), .options.snow=opts) %dopar%{
    
    # row_index <- 1
    row = metadata_filtered[row_index, ]
    
    ## Two possible download paths:
    download_link <- paste0("https://www.encodeproject.org", row$document, "@@download/attachment/", row$biosample_alias, ".pdf")
    alt_download_link <- paste0("https://www.encodeproject.org", row$document, "@@download/attachment/", row$biosample_alias, "_update.pdf")
    
    ## Where to save the file
    file_path <- paste0(row$path, row$biosample_alias, ".pdf")
    row$file_path <- file_path
    
    ## If overwrite results is set to TRUE or if the file does not exists or if
    ## the file is too small, try to download it from the different links. If
    ## none it success, remove the resulted file. We also modify the column
    ## "file_path" if the file was successfully created.
    if (overwrite_results || !file.exists(file_path) || file.info(file_path)$size < 10) {
      tryCatch({
        download.file(download_link, file_path, method = "wget", extra = "--quiet --no-check-certificate")
      }, error = function(e){
        tryCatch({
          download.file(alt_download_link, file_path, method = "wget", extra = "--quiet --no-check-certificate")
        }, error = function(k){
          file.remove(file_path)
          #message(k)
        })
      })
    } 
    
    if (!file.exists(file_path) || file.info(file_path)$size < 10) {
      row$file_path <- NA
    }
    
    return(row)
  } %>% dplyr::bind_rows()
  
  if (!silent) close(pb)
  stopCluster(cl)
  
  return(metadata_documents)
}

#' Extract the images from documents
#'
#' Given the data.frame with the sample's information and their characterization
#' document path, this function extracts the images from the documents and keeps
#' the last found image (which always contains the reported knockdown
#' efficiency).
#'
#' @param metadata_documents Data.frame with the information of the downloaded
#'   files.
#' @param overwrite_results (Optional) Whether to download already found files.
#'   Defaults to FALSE.
#'
#' @return Data.frame with the information of the extracted images.
#' @export
ExtractImages <- function(metadata_documents,
                          overwrite_results = FALSE,
                          virtualenv_path){
  
  
  metadata_images <- foreach(row_index = seq(nrow(metadata_documents))) %do%{
    
    # row_index <- 1
    row = metadata_documents[row_index, ]
    
    biosample = row$biosample
    file_path <- row$file_path
    image_path <- paste0(row$path, biosample, "_Western_Blot_Analysis.png")
    
    ## If the file does not exists return the unmodified row
    if (!file.exists(file_path)) return(row)
    
    ## If overwrite results is set to TRUE or the image does not exists, execute
    ## a python command to extract the figures. All the figures must have a .png
    ## or .jpeg extension. We remove all images but the last one, which contains
    ## the Western Blotting results.
    if (overwrite_results || !file.exists(image_path)) {
      
      system2(command = path.expand(path = paste0(virtualenv_path,"/bin/python")),
              args = c("-m pymupdf extract -images", file_path, "-output", row$path))
      
      short_images <- list.files(row$path, pattern = ".png|.jpeg")
      long_images <- list.files(row$path, pattern = ".png|.jpeg", full.names = T)
      
      last_image <- sapply(short_images, function(x) str_extract(x, "\\d+") %>% as.numeric()) %>% max(na.rm = T)
      
      for(i in seq(length(short_images))){
        image_name <- short_images[i]
        file <- long_images[i]
        
        if(!grepl(last_image, image_name)){
          file.remove(file)
        }else{
          file.rename(file, image_path)
        }
      }
    }
    
    row$image_path <- image_path 
    row
  } %>% dplyr::bind_rows()
  
  return(metadata_images)
}

#' Resize an image
#'
#' @param image PIL image object.
#' @param resize_perc Final size of the output image in percentage. Defaults to
#'   0.25.
#' @param min_width Minimum width to apply the resize. If the image is too
#'   small, the text recognition algorithm will fail. Defaults to 600px.
#'
#' @return PIL image object.
#' @export
ResizeImage <- function(image, 
                        resize_perc = 0.25,
                        min_width = 600){
  if(image$width > 600){
    image_small <- image_cropped$resize(list((image_cropped$width*resize_perc) %>% as.integer(), 
                                             (image_cropped$height*resize_perc) %>% as.integer()), 
                                        PIL$Image$Resampling$LANCZOS)
  }else{
    image_small <- image
  }
  
  return(image_small)
}

#' Summarize and write to disk the Knockdown Efficiency
#'
#' @param metadata_kEff Data.frame containing the extracted efficiencies for
#'   both methods.
#' @param method Whether to extract the Western blotting or qRT-PCR results.
#'   Valid values are: "WB" and "PCR". Defaults to "WB".
#' @param output_file (Optional) Path to where the resulted data.frame will be
#'   stored
#'
#' @return Data.frame with the summarized knockdown efficiencies reported by
#'   ENCODE with either qRT-PCR or Western blotting.
#' @export
WriteEfficiencyTable <- function(metadata_kEff,
                                 method = "WB",
                                 output_file = ""){
  if(method == "WB"){
    metadata_kEff_output <- metadata_kEff %>%
      dplyr::select(target_gene, WB_HepG2, WB_K562) %>%
      tidyr::pivot_longer(c(WB_HepG2, WB_K562)) %>%
      dplyr::mutate(name = str_replace_all(name, "WB_", ""))
  }else if(method == "PCR"){
    metadata_kEff_output <- metadata_kEff %>%
      dplyr::select(target_gene, PCR_HepG2, PCR_K562) %>%
      tidyr::pivot_longer(c(PCR_HepG2, PCR_K562)) %>%
      dplyr::mutate(name = str_replace_all(name, "PCR_", ""))
  }
  
  metadata_kEff_output <- metadata_kEff_output %>%
    dplyr::group_by(target_gene, name) %>%
    dplyr::summarize(kEff = ifelse(all(is.na(value)), NA, max(value, na.rm = T))) %>%
    dplyr::mutate(kEff_avg = mean(kEff, na.rm = T)) %>%
    dplyr::mutate_all(~ifelse(is.nan(.), NA, .)) %>%
    tidyr::pivot_wider(c(target_gene, kEff_avg), values_from = kEff, names_from = name, names_prefix = "kEff_")
  
  ## Write to disk
  if(output_file != ""){
    write.table(metadata_kEff_output, output_file, sep = "\t", row.names = F, quote = FALSE)
  }
  
  return(metadata_kEff_output)
}






## HELPER FUNCTIONS - TPM ------------------------------------------------------


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
ExtractTPM <- function(metadata_quantifications,
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
      
      dir_path <- file.path(output_file, RBP, "tpm")
      dir.create(path = dir_path, recursive = T)
      
      saveRDS(object = metadata_TPM_tidy,
              file = file.path(dir_path, paste0(RBP, "_", type, "_tpm.rds")))
      
    }
  }
  
}
