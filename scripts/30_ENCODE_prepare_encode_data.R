prepare_encode_data <- function(metadata,
                                
                                RBP_source_path,
                                results_path,
                                database_path,
                                
                                gtf_path,
                                gtf_version,
                                
                                ENCODE_silencing_series,
                                num_cores = 8) {
  
  
  logger::log_info(paste0(Sys.time(), "\t\t starting junction reading for the RPBs..."))

  junctionReading(metadata = metadata,

                  RBP_source_path = RBP_source_path,
                  results_path = results_path,
                  database_path = database_path,
                  ENCODE_silencing_series,
                  num_cores = 8)
  gc()

  logger::log_info(paste0(Sys.time(), "\t\t starting the creation of Level QC1 file..."))
  createLevelQ1File(database_path = database_path,
                    gtf_path = gtf_path,
                    gtf_version = gtf_version)
  gc()
  
  
  logger::log_info(paste0(Sys.time(), "\t\t starting the creation of the base data per RPB..."))
  createBaseDataPerRBP(metadata = metadata,
                       results_path = results_path,
                       gtf_version = gtf_version,
                       database_path = database_path)
  gc()
  
}




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
#' @param num_cores Number of multiprocessing cores to use. Memory requirements
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
junctionReading <- function(metadata,
                            RBP_source_path,
                            results_path,
                            database_path,
                            ENCODE_silencing_series,
                            num_cores = 4){
  
  logger::log_info("\t Starting the junction reading process.")
  
  
  
  target_genes <- metadata %>% pull(target_gene) %>% unique
  
  all_junc <- map_df(target_genes, function(RBP) {
    
    
    # RBP <- "ARHGEF2"
    # RBP <- target_genes[1]
    message(RBP)
    
    ## Save GLOBAL BASE DATA for the current RBP
    local_results_path <- file.path(results_path, RBP, "/base_data")
    dir.create(path = local_results_path, recursive = T)
    
    c("case", "control") %>% 
      saveRDS(file = file.path(local_results_path, paste0(RBP, "_clusters_used.rds")))
    
    metadata %>% 
      filter(target_gene == RBP) %>%
      mutate(cluster = experiment_type,
             SRA_project = paste0(ENCODE_silencing_series, " ENCODE")) %>%
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
      cl <- parallel::makeCluster(num_cores)
      doParallel::registerDoParallel(cl)
      logger::log_info("\t\t Reading all extracted BAM files (num_cores = ", num_cores, ").")
      
      all_junc <- foreach(i = 1:nrow(cluster_metadata %>% distinct(sample_id)), 
                          .export=c("RBP","RBP_source_path"),
                          .combine = 'rbind', .packages = c("tidyverse")) %dopar% {
        
        ## i = 1
        ## Definition of the variables
        sample_id <- cluster_metadata[i, ] %>% pull(sample_id)
        junc_path <- file.path(RBP_source_path, RBP, cluster, paste0(sample_id, ".bam.sort.s0.junc"))
        
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
    saveRDS(paste0(database_path, "/all_split_reads_raw.rds"))
  
}



createLevelQ1File <- function(database_path,
                              gtf_path,
                              gtf_version){
  
  logger::log_info("\t Creating the 'all_split_reads_qc_level1.rds' file...")
  

  all_junc_combined <- readRDS(paste0(database_path, "all_split_reads_raw.rds")) %>% as_tibble()
  all_junc_combined %>% head()
  
  
  ## Remove junctions shorter than 25bp
  all_junc_combined <- removeShortJunctions(all_junc_combined)
  logger::log_info("\t Split reads shorter than 25bp removed!")
  gc()
  
  
  ## Convert to a GRanges before blacklist and annotating
  all_junc_combined <- all_junc_combined %>%
    GenomicRanges::GRanges()
  
  
  blacklist_path <- "/home/grocamora/RytenLab-Research/Additional_files/hg38-blacklist.v2.bed"
  encode_blacklist_hg38 <- loadEncodeBlacklist(blacklist_path)
  all_junc_combined <- removeEncodeBlacklistRegions(all_junc_combined, encode_blacklist_hg38)
  logger::log_info("\t Split reads overlapping blacklist regions removed!")
  
  ## Annotate Dasper
  
  edb <- loadEdb(gtf_path)
  all_junc_combined_annotated <- annotateDasper(GRdata = all_junc_combined, edb) %>% 
    tibble::as_tibble()
  logger::log_info("\t Split reads annotation finished!")
  
  ## Remove junctions with ambiguous genes
  all_junc_combined_annotated <- removeAmbiguousGenes(all_junc_combined_annotated)
  logger::log_info("\t Ambiguous split reads removed!")
  
  ## Remove uncategorized junctions
  all_junc_combined_annotated <- removeUncategorizedJunctions(input_SR_details = all_junc_combined_annotated)
  logger::log_info("\t Uncategorised split reads removed!")
  
  ## Save the junctions
  
  dir.create(path = file.path(database_path, gtf_version), recursive = T)
  all_junc_combined_annotated %>% 
    dplyr::select("junID", "seqnames", "start", "end", "width", "strand",
                  gene_id = "gene_id_junction", "in_ref", "type", "tx_id_junction") %>%
    dplyr::bind_rows() %>% 
    dplyr::distinct() %>% 
    mutate(junID = paste0("chr", seqnames,":", start, "-", end, ":", strand)) %>%
    saveRDS(file = file.path(database_path, gtf_version, "all_split_reads_qc_level1.rds"))
  
  logger::log_info("\t 'all_split_reads_qc_level1.rds' file saved!")
}



createBaseDataPerRBP <- function(metadata,
                                 results_path,
                                 gtf_version,
                                 database_path) {
  
 
  all_split_reads_qc_level1 <- readRDS(file = file.path(database_path, gtf_version, "all_split_reads_qc_level1.rds")) 
  
  target_genes <- metadata %>% pull(target_gene) %>% unique()
  
  for(RBP in target_genes) {
    
    # RBP <- target_genes[3]
    logger::log_info(paste0("\t\t Working with '", RBP, "' gene..."))
    
    
    local_results_path <- file.path(results_path, RBP,  "/base_data")
    dir.create(path = local_results_path, recursive = T)
    
    
    for(cluster in c("case", "control")) {
      
      # cluster <- c("case", "control")[1]
      
      
      
      if ( file.exists(file.path(local_results_path, paste0(RBP, "_", cluster, "_all_split_reads_raw.rds"))) ) {
        
        logger::log_info("\t Working with ", RBP, " junctions....")
        
        
        cluster_metadata <- metadata %>% 
          filter(target_gene == RBP,
                 experiment_type == cluster)
        
        
        
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
                     by = c("seqnames",
                            "start",
                            "end")) %>%
          dplyr::relocate(strand, .after="end") %>%
          dplyr::relocate(junID)
        
        
        ## Save the split read counts
        split_read_counts <- all_junc_tidy  %>% 
          filter(sampleID %in% samples_used) %>%
          dplyr::select(junID, reads, sampleID) %>%
          tidyr::pivot_wider(values_from = reads, names_from = sampleID) %>% 
          mutate(across(everything(), .fns = ~replace_na(.,0))) %>%
          saveRDS(file = file.path(local_results_path, paste0(RBP, "_", cluster, "_split_read_counts.rds")))
        
        
        all_junc_tidy %>%
          dplyr::bind_rows() %>% 
          dplyr::distinct() %>%
          dplyr::distinct(junID, .keep_all = T) %>%
          dplyr::select(-sampleID) %>%
          saveRDS(file = file.path(local_results_path, paste0(RBP, "_", cluster, "_all_split_reads.rds")))
        
      } else {
        logger::log_info("\t ", RBP, " does not have junctions.")
      }
      
      
    }
  }
  
}



#' Loads the reference genome into memory
#'
#' The different versions can be downloaded
#' \href{http://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/}{here}.
#'
#' @param gtf_path Path to the reference genome GTF file.
#'
#' @return Connection to the reference genome DB.
#' @export
loadEdb <- function(gtf_path) {
  if (!exists("edb")) {
    logger::log_info("\t\t Loading the reference genome.")
    edb <<- ensembldb::ensDbFromGtf(gtf_path, outfile = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
    edb <<- ensembldb::EnsDb(x = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
  } else {
    logger::log_info("\t\t Variable 'edb' already loaded!")
  }
  
  return(edb)
}

#' Loads the ENCODE blacklisted regions into memory
#'
#' The different versions can be downloaded
#' \href{https://github.com/Boyle-Lab/Blacklist/tree/master/lists}{here}.
#'
#' @param blacklist_path Path to the ENCODE blacklisted regions BED file.
#'
#' @return The ENCODE blacklisted region in GRanges object.
#' @export
loadEncodeBlacklist <- function(blacklist_path) {
  if (!exists("encode_blacklist_hg38")) {
    logger::log_info("\t\t Loading the v2 ENCODE blacklisted regions.")
    encode_blacklist_hg38 <<- rtracklayer::import(blacklist_path) %>% diffloop::rmchr()
  } else {
    logger::log_info("\t\t Variable 'encode_blacklist_hg38' is already loaded!")
  }
  
  return(encode_blacklist_hg38)
}



#' Remove the junctions from the ENCODE blacklisted regions
#'
#' @param GRdata GRanges class object with the relevant junctions.
#' @param encode_blacklist_hg38 GRanges class object with the blacklisted regions.
#'
#' @return Junctions that do not overlap with the blacklisted regions.
#' @export
removeEncodeBlacklistRegions <- function(GRdata,
                                         encode_blacklist_hg38) {
  logger::log_info("\t\t Removing junctions from ENCODE blacklisted regions.")
  ## Look fot the overlaps between GRdata and the ENCODE blacklisted region
  overlaps <- GenomicRanges::findOverlaps(
    query = encode_blacklist_hg38,
    subject = GRdata,
    ignore.strand = F,
    type = "any"
  )
  
  idxs <- S4Vectors::subjectHits(overlaps)
  
  ## If an overlap is found, remove the junctions
  if (length(idxs) > 0) {
    #logger::log_info(paste0("A total of ", length(unique(idxs)), " junctions overlap with an ENCODE blacklisted region."))
    GRdata <- GRdata[-idxs, ]
    #logger::log_info(paste0("A total of ", length(GRdata), " are left after the removal."))
  } else {
    logger::log_info("No junctions overlapping with an ENCODE blacklist region.")
  }
  
  return(GRdata)
}

#' Annotate and categorize every junction
#'
#' @param GRdata GRanges class object with the relevant junctions.
#' @param edb The connection to the reference genome DB.
#'
#' @return Annotated input by dasper.
#' @export 
annotateDasper <- function(GRdata, edb) {
  logger::log_info("\t\t Annotating using dasper::junction_annot().")
  GRdata <- dasper::junction_annot(GRdata,
                                   ref = edb,
                                   ref_cols = c("gene_id", "gene_name", "symbol", "tx_id"),
                                   ref_cols_to_merge = c("gene_id", "gene_name", "tx_id")
  )
  
  return(GRdata)
}

#' Remove the junctions shorter than 25bp.
#'
#' @param input_SR_details Dataframe object with the relevant junctions.
#'
#' @return Junctions bigger than 25bp.
#' @export
removeShortJunctions <- function(input_SR_details) {
  logger::log_info("\t\t Removing junctions shorter than 25bp.")
  output_SR_details <- input_SR_details %>%
    dplyr::filter(width >= 25)
  
  return(output_SR_details)
}

#' Remove the junctions with uncategorized classification.
#'
#' @param input_SR_details Dataframe object with the relevant junctions.
#'
#' @return Junctions categorized as "annotated", "novel_donor" or
#'   "novel_acceptor".
#' @export
removeUncategorizedJunctions <- function(input_SR_details) {
  logger::log_info("\t\t Removing junctions not classified as novel_acceptor, novel_donor or annotated.")
  output_SR_details <- input_SR_details %>%
    dplyr::filter(type %in% c("novel_acceptor", "novel_donor", "annotated"))
  
  return(output_SR_details)
}

#' Remove the junctions from ambiguous genes.
#'
#' @param input_SR_details Dataframe object with the relevant junctions.
#'
#' @return Junctions assigned to only one gene.
#' @export
removeAmbiguousGenes <- function(input_SR_details) {
  logger::log_info("\t\t Removing junctions associated to more than one gene.")
  
  
  output_SR_details <- input_SR_details[which(sapply(input_SR_details$gene_id_junction, length) == 1), ]
  
  return(output_SR_details)
}


