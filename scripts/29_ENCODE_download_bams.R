##############################################################
## CODE Adapted from:
## https://github.com/guillermo1996/ENCODE_Metadata_Extraction
##############################################################


ENCODEDownloadBams <- function(metadata,
                               results.path,
                               regtools.path,
                               samtools.path) {
  

  ## Script parameters ----
  
  ### Whether to write files to disk and to overwrite previously created files. If
  ### rw_disk is set to TRUE, the files exists and overwrite is set to FALSE, the
  ### variables are loaded from disk.
  
  overwrite <- T
  download_cores = 1
  samtools_threads = 1
  samtools_memory = "1G"
  
  
  ## Load metadata ----
  target_RBPs <- metadata %>% dplyr::pull(target_gene) %>% unique()
  
  
  # # 2. Download the BAM files ---
  #logger::log_info("Starting the download BAM files process....")
  
  for (target_RBP in target_RBPs) {
    
    # target_RBP <- target_RBPs[1]
    logger::log_info("Target RBP: ", target_RBP)
    
    ## Target RBP variables
    RBP_metadata <- metadata %>% dplyr::filter(target_gene == target_RBP) %>% dplyr::arrange(experiment_type)
    RBP_path <- file.path(results.path, target_RBP, "/")
    RBP_clusters <- RBP_metadata %>% dplyr::pull(experiment_type) %>% unique()
    
    ## If the junctions from all sample experiments have already been extracted, move on to the next RBP
    if (nrow(CheckDownloadedFiles(RBP.metadata = RBP_metadata, RBP.path = RBP_path)) == nrow(RBP_metadata)) {
      logger::log_info("Ignoring download and extraction. All junctions already extracted!")
      next;
    }
    
    ## Create the subfolders
    CreateSubFolders(RBP.metadata = RBP_metadata,
                     RBP.path = RBP_path,
                     RBP.clusters = RBP_clusters,
                     generate_script = T)
    
    
    DownloadExtractBamFiles(RBP.metadata = RBP_metadata,
                            RBP.path = RBP_path,
                            num.cores = download_cores,
                            samtools.threads = samtools_threads,
                            samtools.memory = samtools_memory,
                            samtools.path = samtools.path,
                            regtools.path = regtools.path,
                            overwrite = overwrite)
  }

  
}


############################################
## FUNCTIONS
############################################

#' Creates the subfolders for a particular RBP
#'
#' @param RBP.metadata Dataframe containing all the metadata for the RBPs/NMDs.
#' @param RBP.path Path to the folder where the generated files will be stored.
#' @param RBP.clusters The different clusters found for the target gene.
#' @param generate_script Whether to generate a download script for the samples
#'   by cluster. By default, TRUE.
#'
#' @return NULL
#' @export
CreateSubFolders <- function(RBP.metadata,
                             RBP.path, 
                             RBP.clusters,
                             generate_script = T) {
  logger::log_info("Generating the folder structure in ", RBP.path)
  ## Loop through the clusters and generate a specific folder for each one.
  for (cluster in RBP.clusters) {
    cluster_path <- paste0(RBP.path, cluster, "/")
    cluster_metadata <- RBP.metadata %>% filter(experiment_type == cluster)
    
    dir.create(cluster_path, recursive = T, showWarnings = F)
    
    ## If generate_script is set to TRUE, a download script will be generated in
    ## the cluster folder. However, it is not recommended to use this script to
    ## download the files.
    if (generate_script == T) {
      GenerateDownloadScript(cluster_metadata, cluster_path)
    }
  }
}

#' Generates the download script
#'
#' Given the metadata for the current target gene and cluster (either case or
#' control), it generates a script to download the BAM files in the
#' corresponding location. The script can be executed to download the BAM files,
#' but it is not recommended.
#'
#' @param cluster_metadata Dataframe containing all the metadata for the current
#'   cluster and target gene.
#' @param cluster_path Path to the cluster folder where the generated files will
#'   be stored.
#'
#' @return NULL
#' @export
GenerateDownloadScript <- function(cluster_metadata,
                                   cluster_path) {
  logger::log_info("Generating the download scripts in ", cluster_path)
  
  ## Generate the script's text
  download_script <- "# Script to download the .bam files\n"
  for (i in 1:nrow(cluster_metadata)) {
    download_link <- GetDownloadLinkMetadata(cluster_metadata[i, ])
    download_script <- paste0(download_script, "wget -c ", download_link, "\n")
  }
  
  ## Writing of the script.
  tryCatch(
    {
      file_path <- paste0(cluster_path, "/download.sh")
      file_conn <- file(file_path)
      writeLines(download_script, file_conn)
      close(file_conn)
      
      system2(command = "chmod", args = c("+x", file_path))
    },
    error = function(e) print(e)
  )
}

#' Generates the download link for a given sample
#'
#' @param cluster_metadata_row Row containing all the metadata for the current
#'   cluster and target gene.
#'
#' @return Download link for the sample.
#' @export
GetDownloadLinkMetadata <- function(cluster_metadata_row){
  sample_id <- cluster_metadata_row %>% pull(sample_id)
  sample_format <- cluster_metadata_row %>% pull(file_format)
  
  return(paste0("https://www.encodeproject.org/files/", sample_id, "/@@download/", sample_id, ".", sample_format))
}

#' Download the BAM files and extracts them into JUNC file
#'
#' Given the metadata for the current cluster (either case or control), this
#' function reads the download link and automatically start the process of
#' downloading and extracting the BAM files. The process requires the path for
#' \href{http://www.htslib.org/}{samtools} and
#' \href{https://regtools.readthedocs.io/en/latest/}{regtools}, but they can be
#' left empty if they are already in the default PATH.
#'
#' It is possible to control certain parameters from samtools, like the number
#' of threads or the maximum memory per core. More information in the official
#' samtools \href{http://www.htslib.org/doc/samtools-sort.html}{documentation}.
#'
#' Tested on samtools 1.16.1 and regtools 0.5.2
#'
#' @param RBP.metadata Dataframe containing all the metadata for the RBPs/NMDs.
#'   It is required to be have the field "file_format", "experiment_type" and
#'   "sample_id" in every row.
#' @param RBP.path Path to the folder where the generated files will be stored.
#' @param num.cores Number of multiprocessing cores to use. Memory requirements
#'   significantly increase with the number of cores.
#' @param samtools.threads Number of threads to use in the samtools sort
#'   command.
#' @param samtools.memory Maximum memory per core to use in the samtools sort
#'   command.
#' @param samtools.path Path to the samtools executable. Can be left empty if
#'   samtools is in default PATH.
#' @param regtools.path Path to the regtools executable. Can be left empty if
#'   regtools is in default PATH.
#' @param overwrite Whether to overwrite previously generated results from the
#'   function. If set to FALSE and 'rw_disk' is set to TRUE, the function looks
#'   for the files in memory and loads them if possible. By default, FALSE.
#'
#' @return NULL
#' @export
DownloadExtractBamFiles <- function(RBP.metadata,
                                    RBP.path,
                                    num.cores,
                                    samtools.threads,
                                    samtools.memory,
                                    samtools.path,
                                    regtools.path,
                                    overwrite = F) {
  
  
  
  ## Only download BAM files for the samples not yet downloaded and junction extracted
  RBPs_BAM_to_download <- CheckBAMFilesToDownload(RBP.metadata, RBP.path)
  
  if (nrow(RBPs_BAM_to_download) > 0) {
    
    logger::log_info("Starting the download process....")
    
    ## Multiprocessing generation. Add argument "output.file" to
    ## parallel::makeCluster() if you want to output information about the parallel
    ## execution
    cl <- parallel::makeCluster(num.cores, outfile = "")
    doParallel::registerDoParallel(cl)
    
    ## Multiprocessing loop
    metrics <- foreach(i = 1:nrow(RBPs_BAM_to_download), .export = "GetDownloadLinkMetadata", .packages = "tidyverse") %dopar% {
      #i <- 1
      ## Definition of the variables
      sample_id <- RBPs_BAM_to_download[i, ] %>% pull(sample_id)
      sample_cluster <- RBPs_BAM_to_download[i, ] %>% pull(experiment_type)
      sample_target_gene = RBPs_BAM_to_download[i, ] %>% pull(target_gene)
      download_link <- GetDownloadLinkMetadata(RBPs_BAM_to_download[i, ])
      file_path <- paste0(RBP.path, sample_cluster, "/", sample_id, ".bam")
      
      ## Check if file BAM file is already extracted
      if(!overwrite & file.exists(file_path)) {
        return(paste0("Ignoring extraction and download. Junction file for sample ", sample_id, 
                      " (", sample_target_gene, " - ", sample_cluster, ") already found!"))
      }
      
      ## Download the BAM file
      #logger::log_info("Starting download of sample ", sample_id, " (", sample_target_gene, " - ", sample_cluster, ")")
      download.file(download_link, file_path, method = "wget", extra = "-c --quiet --no-check-certificate")
      
      
    }
    ## Stop the parallel cluster
    parallel::stopCluster(cl)
    
    ## Print the metrics
    #invisible(lapply(metrics, function(x) logger::log_info(x[1])))
    #invisible(lapply(metrics, function(x) if (!is.na(x[2])) logger::log_info(x[2])))
  }
  
  RBPs_jxn_to_extract <- CheckJxnFilesToExtract(RBP.metadata, RBP.path) %>% arrange(sample_id)
  
  if (nrow(RBPs_jxn_to_extract) > 0) {
    
    logger::log_info("Starting the extraction process...")
    
    ## Multiprocessing generation. Add argument "output.file" to
    ## parallel::makeCluster() if you want to output information about the parallel
    ## execution
    cl <- parallel::makeCluster(num.cores, outfile = "")
    doParallel::registerDoParallel(cl)
    
    ## Multiprocessing loop
    metrics <- foreach(i = 1:nrow(RBPs_jxn_to_extract), .packages = "tidyverse") %dopar% {
      #i <- 1
      ## Definition of the variables
      sample_id <- RBPs_jxn_to_extract[i, ] %>% pull(sample_id)
      sample_cluster <- RBPs_jxn_to_extract[i, ] %>% pull(experiment_type)
      sample_target_gene = RBPs_jxn_to_extract[i, ] %>% pull(target_gene)
      
      file_path <- paste0(RBP.path, sample_cluster, "/", sample_id, ".bam")
      sort_path <- paste0(RBP.path, sample_cluster, "/", sample_id, ".bam.sort")
      junc_path <- paste0(RBP.path, sample_cluster, "/", sample_id, ".bam.sort.s0.junc")

      if (!file.exists(file_path)) { stop(".bam file does not exist!")}
      
      ## Writing of the script.
      tryCatch(
        {
          ## Sort the BAM file using samtools
          logger::log_info("Starting sorting process of sample ", sample_id, " (", sample_target_gene, " - ", sample_cluster, ")")
          
          system2(command = file.path(samtools.path, "samtools"), args = c(
            "sort", file_path,
            "-o", sort_path,
            "--threads", samtools.threads))#,
            #"-m", samtools.memory
          #))
          
          if(!file.exists(sort_path)) stop("Sorting failed")
          
          logger::log_info("Starting the indexing process...")
          system2(command = file.path(samtools.path, "samtools"), args = c(
            "index", sort_path,
            "--threads", samtools.threads
          ))
          
          if(!file.exists(sort_path)) stop("Indexing failed")
          
          logger::log_info("Starting the junction extraction process...")
          system2(command = file.path(regtools.path, "regtools"), args = c(
            "junctions extract ", sort_path,
            "-m 25",
            "-M 1000000",
            "-s XS",
            "-o", junc_path
          ))

        },
        error = function(e) message("Error: ", e)
      )
      
      if (!file.exists(junc_path)) {
        return(c(paste0("Error extracting sample ", sample_id, " (", sample_target_gene, " - ", sample_cluster, ")", ".")))
      } else {
        ## Remove the files that are not necessary
        system2(command = "rm", args = c(file_path))
        system2(command = "rm", args = c(sort_path))
        system2(command = "rm", args = c(paste0(sort_path, ".bai")))
        return(paste0("Successfully extracted the junctions. All intermediary files (BAM included) removed!"))
      }
      
    }
    ## Stop the parallel cluster
    parallel::stopCluster(cl)
    
    ## Print the metrics
    invisible(lapply(metrics, function(x) logger::log_info(x[1])))
    invisible(lapply(metrics, function(x) if (!is.na(x[2])) logger::log_info(x[2])))
  }

}

#' Checks if all the JUNC files in the input dataframe already exists
#'
#' @param RBP.metadata Dataframe containing all the metadata for the RBPs/NMDs.
#' @param RBP.path Path to the folder where the generated files should be
#'   stored.
#'
#' @return Whether the JUNC files for the input metadata exists or not.
#' @export
CheckDownloadedFiles <- function(RBP.metadata,
                                 RBP.path){
  
  file_names <- apply(RBP.metadata, 1, function(x) {
    paste0(RBP.path, x["experiment_type"], "/", x["sample_id"], ".bam.sort.s0.junc")
  })
  
  return(RBP.metadata[which(file.exists(file_names)),])
}

#' Checks BAM files to download from the ENCODE platform
#'
#' @param RBP.metadata Dataframe containing all the metadata for the RBPs/NMDs.
#' @param RBP.path Path to the folder where the generated files should be
#'   stored.
#'
#' @return Whether the JUNC files for the input metadata exists or not.
#' @export
CheckBAMFilesToDownload <- function(RBP.metadata,
                                    RBP.path){
  
  file_names <- apply(RBP.metadata, 1, function(x) {
    paste0(RBP.path, x["experiment_type"], "/", x["sample_id"], ".bam")
  })
  
  return(RBP.metadata[which(!file.exists(file_names)),])
}

CheckJxnFilesToExtract <- function(RBP.metadata,
                                   RBP.path){
  
  file_names <- apply(RBP.metadata, 1, function(x) {
    paste0(RBP.path, x["experiment_type"], "/", x["sample_id"], ".bam.sort.s0.junc")
  })
  
  return(RBP.metadata[which(!file.exists(file_names)),])
}
