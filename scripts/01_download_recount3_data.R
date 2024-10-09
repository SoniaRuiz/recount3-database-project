#' Title
#' downloads, quality-control and annotates the split read data from a given recount3 project
#' @param recount3.project.IDs array with the recount3 projects to download
#' e.g. "SRP009615"
#' @param project.name name given locally to the recount3 project.
#' A given recount3 project can be separated in multiple independent projects in recount3.
#' For instance, GTEx is stored in recount3 in multiple independent projects ID (e.g. BRAIN, SKIN, BLOOD, etc), 
#' but all of them belong to GTEx. Hence, it is useulf to have a folder named "project.name" that will contain
#' multiple subfolders named as the elements contained within the object 'recount3.project.IDs'
#' @param gtf.version Ensembl version, e.g. "105"
#' @param data.source source of the data within recount3. It can be:
#' - "data_sources/sra"  
#' - "data_sources/gtex" 
#' - "data_sources/tcga"
#'
#' @return
#' @export
#'
#' @examples
DownloadRecount3Data <- function (recount3.project.IDs,
                                  project.name,
                                  gtf.version,
                                  blacklist.path,
                                  gtf.path,
                                  data.source,
                                  database.folder,
                                  results.folder) {
  
  
  ##########################################################
  ## Read all the split reads and return them by tissue
  ##########################################################
  
  if (file.exists(file.path(dirname(database.folder), "all_split_reads_raw.rds"))) {
    
    logger::log_info("\t Project '", project.name, "' --> Loading 'all_split_reads_raw.rds' file...")
    all_split_reads_raw <- readRDS(file = file.path(dirname(database.folder), "all_split_reads_raw.rds"))
    
  } else {
   
    doParallel::registerDoParallel(2)
    all_split_reads_raw <- foreach(j = seq(length(recount3.project.IDs)), .combine = "rbind") %do% {
      
      # j <- 1
      project_id <- recount3.project.IDs[j]
      
      # project_id <- recount3.project.IDs[1]
      
      logger::log_info("\t --> Getting data from '", project_id, "' recount3 project...")
      
      folder_root <- file.path(results.folder, project_id, "base_data/")
      dir.create(file.path(folder_root), recursive = TRUE, showWarnings = T)
      
      #############################################################################
      
      ## THERE ARE OTHER RECOMMENDED WAYS OF DOWNLOADING DATA FROM RECOUNT3
      ## HOWEVER, THIS ALTERNATIVE METHOD IS FOLLOWED IN ORDER TO REDUCE MEMORY USAGE
      ## PARTICULARLY USEFUL WITH LARGE DATASETS SUCH AS GTEX
      ## FOR OTHER RECOMMENDED METHODS, SEE (https://bioconductor.org/packages/release/bioc/manuals/recount3/man/recount3.pdf) - ACCESSED 08/07/2023
      
      #############################################################################
      
      
      jxn_files <- recount3::locate_url(
        project = project_id,
        project_home = data.source,
        type = "jxn",
        organism = "human",
        annotation = "gencode_v29",
        #jxn_format = c("UNIQUE"),
        recount3_url = getOption("recount3_url", "http://duffel.rail.bio/recount3")
      )
      
      feature_info <- utils::read.delim(recount3::file_retrieve(
        url = jxn_files[grep("\\.RR\\.gz$", jxn_files)],
        bfc = recount3::recount3_cache(),
        verbose = getOption("recount3_verbose", TRUE)
      ))
      
      ## Convert unstranded junctions from '*' to '?' to facilitate later conversion to GRanges
      feature_info$strand[feature_info$strand == "?"] <- "*"
      
      all_split_reads <- data.frame(junID = feature_info %>% GRanges() %>% as.character(),
                                    chr = feature_info$chromosome,
                                    start = feature_info$start,
                                    end = feature_info$end,
                                    strand = feature_info$strand,
                                    width = feature_info$length,
                                    annotated = feature_info$annotated,
                                    left_motif = feature_info$left_motif,
                                    right_motif = feature_info$right_motif) %>% as_tibble()
      
      logger::log_info("\t ", project_id, " --> ", all_split_reads %>% nrow(), " split reads.")

      return(all_split_reads %>% distinct(junID, .keep_all = T) %>% mutate(recount_project = project_id))
    }
   
    logger::log_info(all_split_reads_raw %>% nrow(), " getting unique split reads across all projects....")
    
    ## To keep track of project frequency per jxn
    all_split_reads_raw <- if (length(all_split_reads_raw$recount_project %>% unique) > 1) {
      all_split_reads_raw %>% 
        dplyr::group_by(junID) %>%
        mutate(n_projects = n()) %>%
        ungroup() %>%
        distinct(junID, .keep_all = T) %>%
        dplyr::select(-recount_project)
    } else {
      all_split_reads_raw %>% 
        mutate(n_projects = 1) %>%
        distinct(junID, .keep_all = T) %>%
        dplyr::select(-recount_project)
    }
    
    logger::log_info(all_split_reads_raw %>% nrow(), " initial number of split reads across all projects.")
    gc()
    
    ## Save data
    dir.create(file.path(database.folder), recursive = TRUE, showWarnings = T)
    saveRDS(object = all_split_reads_raw, file = file.path(dirname(database.folder), "all_split_reads_raw.rds"))
    
  }
  
  #######################################################################
  ## 1. Discard all split reads shorter than 25 bp
  #######################################################################
  
  ## Remove junctions shorter than 25bp
  all_split_reads_raw_tidy <- RemoveShortJunctions(all_split_reads_raw)
  logger::log_info("\t Split reads shorter than 25bp removed!")
  
  rm(all_split_reads_raw)
  
  #######################################################################
  ## 2. Remove split reads located in unplaced sequences in the chromosomes
  #######################################################################
  
  logger::log_info("\t Removing unplaced sequences genome...")
  
  all_split_reads_raw_tidy_gr <- all_split_reads_raw_tidy %>% filter(!(strand=="?")) %>% GenomicRanges::GRanges() %>% diffloop::rmchr()
  
  rm(all_split_reads_raw_tidy)
  
  all_split_reads_raw_tidy_gr %>% length()
  all_split_reads_raw_tidy_gr <- GenomeInfoDb::keepSeqlevels(x = all_split_reads_raw_tidy_gr,
                                                             value = intersect(all_split_reads_raw_tidy_gr %>% GenomeInfoDb::seqnames() %>% levels(), 
                                                                               c( "1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                                                                  "2", "20", "21", "22", "3", "4", "5", "6", "7", "8", "9", "X", "Y")), 
                                                             pruning.mode = "tidy")
  ## Data check
  all_split_reads_raw_tidy_gr %>% head()
  logger::log_info(all_split_reads_raw_tidy_gr %>% length(), " split reads after removing those not aligning the hg38")

  
  ##########################################################
  ## 3. Remove split reads overlapping the ENCODE backlist
  ##########################################################
  
  logger::log_info("\t Removing blacklist sequences...")
  all_split_reads_raw_tidy_gr <- RemoveEncodeBlacklistRegions(GRdata = all_split_reads_raw_tidy_gr, blacklist.path = blacklist.path)
  
  
  ## Log number of split reads
  all_split_reads_raw_tidy_gr %>% as_tibble() %>% distinct(junID, .keep_all = T) %>% nrow() %>% logger::log_info()
  logger::log_info("\t Split reads overlapping blacklist regions removed!")
  
  #######################################################
  ## 4. Anotate using the R package 'dasper'
  #######################################################
  
  logger::log_info("\t Annotating dasper...")
  
  ## Annotate Dasper
  edb <- LoadEdb(gtf.path)
  all_split_reads_details_w_symbol <- AnnotateDasper(GRdata = all_split_reads_raw_tidy_gr %>% GRanges(), edb) %>% tibble::as_tibble()
  logger::log_info("\t Split reads annotation finished!")
  
  rm(all_split_reads_raw_tidy_gr)
  rm(edb)
  
  ################################################################################
  ## 5. Discard all junctions that are not annotated, novel donor or novel acceptor
  ################################################################################
  
  ## Remove uncategorized junctions
  all_split_reads_details_w_symbol <- RemoveUncategorizedJunctions(input.SR.details = all_split_reads_details_w_symbol)
  logger::log_info("\t Uncategorised split reads removed!")
  
  logger::log_info("\t ", all_split_reads_details_w_symbol %>% length(), " annotated, novel donor, novel acceptor and novel combo junctions.")
  
  
  ############################################################################
  ## 6. Discard all ambiguous split reads (i.e. assigned to multiple genes)
  ############################################################################
  
  logger::log_info("\t Discarding ambiguous split reads ...")
  
  ## Remove junctions with ambiguous genes
  all_split_reads_details_w_symbol <- RemoveAmbiguousGenes(input.SR.details = all_split_reads_details_w_symbol, database.folder)
  logger::log_info("\t Ambiguous split reads removed!")
  
  
  #####################
  ## 7. SAVE RESULTS
  #####################
  
  all_split_reads_details_w_symbol <- all_split_reads_details_w_symbol %>% dplyr::filter(ambiguous == F) %>% dplyr::select(-ambiguous)
  
  logger::log_info("\t ", ambiguous_introns %>% nrow(), " ambiguous split reads discarded.")
  saveRDS(object = all_split_reads_details_w_symbol %>% dplyr::rename(gene_id = gene_id_junction),
          file = file.path(database.folder, "all_split_reads_qc_level1.rds"))
  
  
  ## FREE UP SOME MEMORY
  rm(ambiguous_introns)
  rm(all_split_reads_details_w_symbol)
  gc()
  
}