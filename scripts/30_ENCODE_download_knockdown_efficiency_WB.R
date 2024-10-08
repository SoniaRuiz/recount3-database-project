##############################################################
## CODE Adapted from:
## https://github.com/guillermo1996/ENCODE_Metadata_Extraction
##############################################################


library(doSNOW)
library(reticulate)


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
                                          results.path) {
  
  
  ## Create and load python environment -------------------------
  
  virtualenv <- "ENCODE_python_env"
  reticulate::virtualenv_create(virtualenv)
  reticulate::use_virtualenv(virtualenv)
  
  py_available()
  py_install("pytesseract") 
  py_install("tesseract")
  py_install("pymupdf")
  py_install("fitz") 
  py_install("frontend") 
  
  pytesseract <- import("pytesseract")
  PIL <- import("PIL")
  
  
  ## Set Variables ----------------------------------------------

  results_path <- file.path(results.path, "/WesternBlotting_PCR/RBPs/")
  dir.create(results_path, recursive = T, showWarnings = F)
  metadata_WB_output <- file.path(results_path, "metadata_WB_kEff.tsv")
  metadata_PCR_output <- file.path(results_path, "metadata_PCR_kEff.tsv")
  
  download_cores = 16
  overwrite_results = T
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
  metadata_documents <- DownloadCharacterizationDocuments(metadata_filtered, download_cores, overwrite_results, silent = F)
  
  ## Check the existence of the files
  for(row_index in seq(nrow(metadata_documents))){
    
    row = metadata_documents[row_index, ]
    file_path <- row$file_path
    path <- row$path
    
    if(!file.exists(file_path) || file.info(file_path)$size < 10){
      logger::log_warn("Error for row ", row_index, "! File path ", path)
    }
    
  }
  
  ## Extract the images of all files
  metadata_images <- ExtractImages(metadata_documents, overwrite_results = overwrite_results, virtualenv_path)
  
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
    if(!file.exists(image_path)) return(row)
    
    image <- PIL$Image$open(image_path)
    image_cropped <- image$crop(list(0, image$height*0.9, 0.8*image$width, image$height))
    image_small <- ResizeImage(image, resize_perc, min_width)
    
    image_small$save(paste0(path, "cropped_image.png"))
    text_df <- pytesseract$image_to_string(image = image_small) %>%
      str_replace_all("\\f", "") %>%
      str_split("\\n", simplify = T) %>%
      str_split(" ", simplify = T) %>%
      .[1:2, ] %>%
      as_tibble()
    
    if(ncol(text_df) == 5){
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
    }else if(ncol(text_df) == 3){
      kEff_df <- text_df %>% 
        `colnames<-`(c("method", "cell_line_1", "cell_line_2")) %>%
        dplyr::mutate(across(-method, function(x) as.numeric(sub("%", "", x)))) %>%
        dplyr::mutate(cell_line = rowMeans(dplyr::select(., cell_line_1, cell_line_2))) %>%
        dplyr::select(method, cell_line)
      
      row$WB_HepG2 <- ifelse(row$cell_line == "HepG2", kEff_df %>% dplyr::filter(method == "Western") %>% dplyr::pull(cell_line), NA)
      row$WB_K562 <- ifelse(row$cell_line == "K562", kEff_df %>% dplyr::filter(method == "Western") %>% dplyr::pull(cell_line), NA)
      
      row$PCR_HepG2 <- ifelse(row$cell_line == "HepG2", kEff_df %>% dplyr::filter(method != "Western") %>% dplyr::pull(cell_line), NA)
      row$PCR_K562 <- ifelse(row$cell_line == "K562", kEff_df %>% dplyr::filter(method != "Western") %>% dplyr::pull(cell_line), NA)
    }else{
      logger::WARN("Error in row ", row_index, ". Columns are not valid")
      row$WB_HepG2 <- NA
      row$WB_K562 <- NA
      row$PCR_HepG2 <- NA
      row$PCR_K562 <- NA
    }
    
    return(row)
  } %>% dplyr::bind_rows()
  
  ## Test for consistency between the cell lines ----
  for(target_RBP in target_RBPs){
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




## HELPER FUNCTIONS -------------------------------------------------------------------------------------


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
    RBP_metadata <- metadata_filtered %>% 
      dplyr::filter(target_gene == target_RBP)
    
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
                                              download_cores = 1,
                                              overwrite_results = FALSE,
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
  metadata_documents <- foreach(row_index = seq(nrow(metadata_filtered)), .options.snow=opts) %dopar%{
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
    if(overwrite_results || !file.exists(file_path) || file.info(file_path)$size < 10){
      tryCatch({
        download.file(download_link, file_path, method = "wget", extra = "--quiet --no-check-certificate")
      }, error = function(e){
        tryCatch({
          download.file(alt_download_link, file_path, method = "wget", extra = "--quiet --no-check-certificate")
        }, error = function(k){
          file.remove(file_path)
        })
      })
    } 
    
    if(!file.exists(file_path) || file.info(file_path)$size < 10){
      row$file_path <- NA
    }
    
    return(row)
  } %>% dplyr::bind_rows()
  if(!silent) close(pb)
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
    row = metadata_documents[row_index, ]
    
    biosample = row$biosample
    file_path <- row$file_path
    image_path <- paste0(row$path, biosample, "_Western_Blot_Analysis.png")
    
    ## If the file does not exists return the unmodified row
    if(!file.exists(file_path)) return(row)
    
    ## If overwrite results is set to TRUE or the image does not exists, execute
    ## a python command to extract the figures. All the figures must have a .png
    ## or .jpeg extension. We remove all images but the last one, which contains
    ## the Western Blotting results.
    if(overwrite_results || !file.exists(image_path)){
      system2(command = path.expand(path = paste0(virtualenv_path,"/bin/python")),
              args = c("-m fitz extract -images", file_path, "-output", row$path))
      
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

## 30