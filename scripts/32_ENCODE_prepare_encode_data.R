##############################################################
## CODE Adapted from:
## https://github.com/guillermo1996/ENCODE_Metadata_Extraction
##############################################################

#' Title
#' Prepares the junction extraction, QC and annotation from each shRNA RBP knockdown experiment from the ENCODE platform
#' @param metadata 
#' @param RBP.source.path 
#' @param results.path 
#' @param database.path 
#' @param gtf.path 
#' @param gtf.version 
#' @param blacklist.path 
#' @param ENCODE.silencing.series 
#' @param num.cores 
#'
#' @return
#' @export
#'
#' @examples
PrepareEncodeData <- function(metadata,
                              RBP.source.path,
                              results.path,
                              database.path,
                              gtf.path,
                              gtf.version,
                              blacklist.path,
                              ENCODE.silencing.series,
                              num.cores = 8) {
  
  
  logger::log_info(paste0(Sys.time(), "\t\t starting junction reading for the RPBs..."))
  
  JunctionReading(metadata = metadata,
                  RBP.source.path = RBP.source.path,
                  results.path = results.path,
                  database.path = database.path,
                  ENCODE.silencing.series,
                  num.cores = 8)
  gc()
  
  logger::log_info(paste0(Sys.time(), "\t\t starting the creation of Level QC1 file..."))
  CreateLevelQ1File(database.path = database.path,
                    blacklist.path,
                    gtf.path = gtf.path,
                    gtf.version = gtf.version)
  gc()
  
  
  logger::log_info(paste0(Sys.time(), "\t\t starting the creation of the base data per RPB..."))
  CreateBaseDataPerRBP(metadata = metadata,
                       results.path = results.path,
                       gtf.version = gtf.version,
                       database.path = database.path)
  gc()
  
}




## HELPER FUNCTIONS ----------------------------------------------------------------------



#' Read JUNC files and extract the junctions
#'
#' From the .junc files generated using
#' \href{https://regtools.readthedocs.io/en/latest/commands/junctions-extract/}{regtools
#' junctions extract}, this functions process the columns into the desired
#' format and joins the information from junctions across all samples.
#'
#' One file is generated at the end, which contains information about all the
#' junctions found inside the samples' .junc files. It also contains the reads
#' per sample of each junction.
#'
#' @param metadata Dataframe containing all the metadata.
#' @param main_samples_path Path to where the JUNC files are stored.
#' @param num.cores Number of multiprocessing cores to use. Memory requirements
#'   significantly increase with the number of cores.
#' @param rw_disk Whether to store the results in disk. By default, TRUE.
#' @param overwrite Whether to overwrite previously generated results from the
#'   function. If set to FALSE and 'rw_disk' is set to TRUE, the function looks
#'   for the files in memory and loads them if possible. By default, FALSE.
#'
#' @return Data.frame containing every junction found in all samples. The
#'   dataframe will contain information about the junction and the reads of that
#'   junction in each sample.
#' @export
JunctionReading <- function(metadata,
                            RBP.source.path,
                            results.path,
                            database.path,
                            ENCODE.silencing.series,
                            num.cores = 4){
  
  logger::log_info("\t Starting the junction reading process.")
  
  
  
  target_genes <- metadata %>% pull(target_gene) %>% unique
  
  all_junc <- map_df(target_genes, function(RBP) {
    
    
    # RBP <- "ARHGEF2"
    # RBP <- target_genes[1]
    message(RBP)
    
    ## Save GLOBAL BASE DATA for the current RBP
    local_results_path <- file.path(results.path, RBP, "/base_data")
    dir.create(path = local_results_path, recursive = T)
    
    c("case", "control") %>% 
      saveRDS(file = file.path(local_results_path, paste0(RBP, "_clusters_used.rds")))
    
    metadata %>% 
      filter(target_gene == RBP) %>%
      mutate(cluster = experiment_type,
             SRA_project = paste0(ENCODE.silencing.series, " ENCODE")) %>%
      distinct(target_gene, SRA_project, cluster, sample_id, .keep_all = T) %>%
      saveRDS(file = file.path(local_results_path, paste0(RBP, "_samples_metadata.rds")))
    
    
    ## GET and save LOCAL BASE DATA per cluster within the current RBP
    map_df(c("case", "control"), function(cluster) {
      
      # cluster <- c("case", "control")[1]
      
      message(cluster)
      
      
      cluster_metadata <- metadata %>% 
        filter(target_gene == RBP,
               experiment_type == cluster)
      
      ## Multiprocessing loop
      cl <- parallel::makeCluster(num.cores)
      doParallel::registerDoParallel(cl)
      logger::log_info("\t\t Reading all extracted BAM files (num.cores = ", num.cores, ").")
      
      all_junc <- foreach(i = 1:nrow(cluster_metadata %>% distinct(sample_id)), 
                          .export=c("RBP","RBP.source.path"),
                          .combine = 'rbind', .packages = c("tidyverse")) %dopar% {
                            
                            ## i = 1
                            ## Definition of the variables
                            sample_id <- cluster_metadata[i, ] %>% pull(sample_id)
                            junc_path <- file.path(RBP.source.path, RBP, cluster, paste0(sample_id, ".bam.sort.s0.junc"))
                            
                            logger::log_info("\t Reading junctions from ", RBP)
                            
                            if(!file.exists(junc_path)) return(tibble())
                            
                            ## Read the junction file into a tibble using readr::read_table()
                            ## The "locale" argument is to read comma separated values (i.e. 25,12)
                            junc <- readr::read_table(
                              junc_path,
                              col_names = F,
                              show_col_types = F,
                              progress = F,
                              locale = readr::locale(grouping_mark = "")
                            )
                            
                            ## Transformations of the junctions
                            junc <- junc %>%
                              dplyr::select(
                                chr = X1,
                                start = X2,
                                stop = X3,
                                junID = X4,
                                reads = X5,
                                strand = X6,
                                blockSizes = X11
                              ) %>%
                              dplyr::mutate(strand = ifelse(strand == "?", "*", strand)) %>%
                              dplyr::mutate(sampleID = sample_id) %>%
                              tidyr::separate(col = blockSizes, sep = ",", c("blockSizesStart", "blockSizesEnd"), conver = T) %>%
                              dplyr::mutate(start = start + blockSizesStart + 1, stop = stop - blockSizesEnd) %>%
                              GenomicRanges::GRanges() %>%
                              diffloop::rmchr() 
                            
                            tryCatch(
                              eval(
                                junc %>%
                                  GenomeInfoDb::keepSeqlevels(value = c(seq(1, 22) %>% as.character(), "X", "Y"), pruning.mode = "tidy") %>%
                                  tibble::as_tibble() %>%
                                  dplyr::mutate(
                                    seqnames = as.character(seqnames),
                                    strand = as.character(strand)
                                  ) %>%
                                  dplyr::select(-junID, -blockSizesStart, -blockSizesEnd) %>%
                                  return()
                              ),
                              error = function(e) {
                                return(tibble())
                              }
                            )
                          }
      ## Stop the parallel cluster
      parallel::stopCluster(cl)
      
      # print(all_junc %>% head)
      
      if ( all_junc %>% nrow() > 0 ) {
        
        all_junc %>%
          dplyr::group_by(seqnames, start, end) %>%
          dplyr::mutate(junID = dplyr::cur_group_id(), .before = seqnames) %>%
          dplyr::ungroup() %>%
          saveRDS(file = file.path(local_results_path, paste0(RBP, "_", cluster, "_all_split_reads_raw.rds")))
        
        ## Save the samples used
        samples_used <- all_junc %>%
          distinct(sampleID) %>%
          pull()
        samples_used %>% saveRDS(file = file.path(local_results_path, paste0(RBP, "_", cluster, "_samples_used.rds")))
        
        
        
        return(all_junc)
      } else {
        logger::log_info("\t ", RBP, " does not have junctions.")
        return(NULL)
      }
      
    })
  })
  
  
  ## Save the unique junctions with ID
  all_junc %>%
    distinct(seqnames, start, end, .keep_all = T) %>%
    dplyr::group_by(seqnames, start, end) %>%
    dplyr::mutate(junID = dplyr::cur_group_id(), .before = seqnames) %>%
    dplyr::ungroup() %>%
    dplyr::select(-sampleID) %>%
    saveRDS(paste0(database.path, "/all_split_reads_raw.rds"))
  
}



CreateLevelQ1File <- function(database.path,
                              blacklist.path,
                              gtf.path,
                              gtf.version){
  
  logger::log_info("\t Creating the 'all_split_reads_qc_level1.rds' file...")
  
  
  all_junc_combined <- readRDS(paste0(database.path, "all_split_reads_raw.rds")) %>% as_tibble()
  all_junc_combined %>% head()
  
  
  ## 1. Remove junctions shorter than 25bp
  all_junc_combined <- RemoveShortJunctions(all_junc_combined)
  logger::log_info("\t Split reads shorter than 25bp removed!")
  gc()
  
  
  
  ## 2. Convert to a GRanges before blacklist and annotating
  all_junc_combined <- all_junc_combined %>% GenomicRanges::GRanges()
  
  all_junc_combined <- RemoveEncodeBlacklistRegions(GRdata = all_junc_combined, blacklist.path = blacklist.path)
  logger::log_info("\t Split reads overlapping blacklist regions removed!")
  
  
  
  ## 3. Annotate Dasper
  edb <- LoadEdb(gtf.path)
  all_junc_combined_annotated <- AnnotateDasper(GRdata = all_junc_combined, edb) %>% tibble::as_tibble()
  logger::log_info("\t Split reads annotation finished!")
  
  
  
  
  ## 4. Remove uncategorized junctions
  all_junc_combined_annotated_tidy <- RemoveUncategorizedJunctions(input.SR.details = all_junc_combined_annotated)
  logger::log_info("\t Uncategorised split reads removed!")
  
  
  
  
  ## 5. Remove junctions with ambiguous genes
  dir.create(path = file.path(database.path, gtf.version), recursive = T)
  all_junc_combined_annotated_tidy <- RemoveAmbiguousJunctions(all_junc_combined_annotated_tidy, 
                                                               database.folder = file.path(database.path, gtf.version))
  logger::log_info("\t Ambiguous split reads removed!")
  
  
  
  
  ## 6. Save the junctions
  all_junc_combined_annotated_tidy %>% 
    dplyr::select("junID", "seqnames", "start", "end", "width", "strand",
                  gene_id = "gene_id_junction", "in_ref", "type", "tx_id_junction") %>%
    dplyr::bind_rows() %>% 
    dplyr::distinct() %>% 
    mutate(junID = paste0("chr", seqnames,":", start, "-", end, ":", strand)) %>%
    saveRDS(file = file.path(database.path, gtf.version, "all_split_reads_qc_level1.rds"))
  
  logger::log_info("\t 'all_split_reads_qc_level1.rds' file saved!")
  
}



CreateBaseDataPerRBP <- function(metadata,
                                 results.path,
                                 gtf.version,
                                 database.path) {
  
  
  all_split_reads_qc_level1 <- readRDS(file = file.path(database.path, gtf.version, "all_split_reads_qc_level1.rds")) 
  target_genes <- metadata %>% pull(target_gene) %>% unique()
  
  for (RBP in target_genes) {
    
    # RBP <- target_genes[1]
    logger::log_info(paste0("\t\t Working with '", RBP, "' gene..."))
    
    local_results_path <- file.path(results.path, RBP,  "base_data")
    dir.create(path = local_results_path, recursive = T)
    
    for (cluster in c("case", "control")) {
      
      # cluster <- c("case", "control")[1]
      
      if (file.exists(file.path(local_results_path, paste0(RBP, "_", cluster, "_all_split_reads_raw.rds")))) {
        
        logger::log_info("\t Working with ", RBP, " junctions....")
        
        cluster_metadata <- metadata %>% 
          filter(target_gene == RBP, experiment_type == cluster)
        
        ## Read the samples used
        samples_used <- readRDS(file = file.path(local_results_path, paste0(RBP, "_", cluster, "_samples_used.rds")))
        
        ## Load junctions current RBP
        all_junc <- readRDS(file = file.path(local_results_path, paste0(RBP, "_", cluster, "_all_split_reads_raw.rds")))
        
        ## Filter local junctions to get those from the global object passing level1 QC
        all_junc_tidy <- all_junc %>% 
          distinct(seqnames, start, end, sampleID, reads, .keep_all = T) %>%
          dplyr::select(-strand, -junID) %>%
          inner_join(y = all_split_reads_qc_level1 %>% 
                       dplyr::select("junID", "seqnames", "start", "end", "strand", "gene_id", "in_ref", "type", "tx_id_junction"),
                     by = c("seqnames", "start", "end")) %>%
          dplyr::relocate(strand, .after="end") %>%
          dplyr::relocate(junID)
        
        
        ## Save the split read counts for the current sample cluster
        
        ## 1. Novel donor, novel acceptor and annotated introns
        split_read_counts <- all_junc_tidy  %>% 
          filter(sampleID %in% samples_used,
                 type %in% c("annotated", "novel_donor", "novel_acceptor")) %>%
          dplyr::select(junID, reads, sampleID) %>%
          tidyr::pivot_wider(values_from = reads, names_from = sampleID) %>% 
          mutate(across(everything(), .fns = ~replace_na(.,0)))
        split_read_counts %>%
          saveRDS(file = file.path(local_results_path, paste0(RBP, "_", cluster, "_split_read_counts.rds")))
        ## Save the details about the split reads for the current sample cluster
        all_junc_tidy %>%
          filter(junID %in% split_read_counts$junID) %>%
          dplyr::bind_rows() %>% 
          dplyr::distinct() %>%
          dplyr::distinct(junID, .keep_all = T) %>%
          dplyr::select(-sampleID) %>%
          saveRDS(file = file.path(local_results_path, paste0(RBP, "_", cluster, "_all_split_reads.rds")))
        
        ## 2. Novel combos
        split_read_counts_combos <- all_junc_tidy  %>% 
          filter(sampleID %in% samples_used, type == "novel_combo") %>%
          dplyr::select(junID, reads, sampleID) %>%
          tidyr::pivot_wider(values_from = reads, names_from = sampleID) %>% 
          mutate(across(everything(), .fns = ~replace_na(.,0)))
        split_read_counts_combos %>%
          saveRDS(file = file.path(local_results_path, paste0(RBP, "_", cluster, "_split_read_counts_combos.rds")))
        ## Save the details about the split reads for the current sample cluster
        all_junc_tidy %>%
          filter(junID %in% split_read_counts_combos$junID) %>%
          dplyr::bind_rows() %>% 
          dplyr::distinct() %>%
          dplyr::distinct(junID, .keep_all = T) %>%
          dplyr::select(-sampleID) %>%
          saveRDS(file = file.path(local_results_path, paste0(RBP, "_", cluster, "_all_split_reads_combos.rds")))
        
        
      } else {
        logger::log_info("\t ", RBP, " does not have junctions.")
      }
    }
  }
}

#32