#####################################
## FUNCTIONS - PREPARE RECOUNT3 DATA
#####################################

#' Title
#' downloads, quality-control and annotates the split read data from a given recount3 project
#' @param recount3.project.IDs array with the recount3 projects to download
#' e.g. "SRP009615"
#' @param project.name name given locally to the recount3 project.
#' A given recount3 project can be separated in multiple independent projects in recount3.
#' For instance, GTEx is stored in recount3 in multiple independent projects ID (e.g. BRAIN, SKIN, BLOOD, etc), 
#' but all of them belong to GTEx. Hence, it is useulf to have a folder named "project.name" that will contain
#' multiple subfolders named as the elements contained within the object 'recount3.project.IDs'
#' @param gtf.version Ensembl version (tested using Ensembl v105)
#' e.g. "105"
#' @param data.source source of the data within recount3. It can be:
#' - "data_sources/sra"  
#' - "data_sources/gtex" 
#' - "data_sources/tcga"
#'
#' @return
#' @export
#'
#' @examples
init_recount3_data <- function (recount3.project.IDs,
                                project.name,
                                gtf.version,
                                data.source,
                                database.folder,
                                results.folder) {
  
  
  ##########################################################
  ## Read all the split reads and return them by tissue
  ##########################################################
  
  if ( file.exists(paste0(database.folder, "/all_split_reads_raw.rds")) ) {
    
    print(paste0(Sys.time(), " - loading 'all_split_reads_raw.rds' file..."))
    all_split_reads_raw_tidy <- readRDS(file = paste0(database.folder, "/all_split_reads_raw.rds"))
    
  } else {
    
    ## We access recount3 files directly to save memory
    all_split_reads_raw <- map_df(recount3.project.IDs, function(project_id) {
      
      # project_id <- recount3.project.IDs[1]
      # project_id <- "KIDNEY"
      print(paste0(Sys.time(), " - getting data from '", project_id, "' recount3 project..."))
      
      folder_root <- paste0(results.folder, "/", project_id, "/base_data/")
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
      
      ## Convert unstranded junctions from '*' to '?' to facilitate later conversion to genome:ranges
      feature_info$strand[feature_info$strand == "?"] <- "*"
      #feature_info <- GenomicRanges::GRanges(feature_info)
      
      all_split_reads <- data.frame(chr = feature_info$chromosome,
                                    start = feature_info$start,
                                    end = feature_info$end,
                                    strand = feature_info$strand,
                                    width = feature_info$length,
                                    annotated = feature_info$annotated,
                                    junID = paste0(feature_info$chromosome, ":",
                                                   feature_info$start, "-", feature_info$end, ":",
                                                   feature_info$strand))
  
      
      saveRDS(object = all_split_reads %>% data.table::as.data.table(),
              file = paste0(folder_root, "/all_split_reads_raw.rds"))
 
      return(all_split_reads %>%
               distinct(junID, .keep_all = T))
    })
    
    
    all_split_reads_raw_tidy <- all_split_reads_raw %>%
      distinct(junID, .keep_all = T)
    
    print(paste0(Sys.time(), " - ", all_split_reads_raw_tidy %>% nrow(), " initial number of split reads"))
    gc()
    
    ## Save data
    dir.create(file.path(database.folder), recursive = TRUE, showWarnings = T)
    saveRDS(object = all_split_reads_raw_tidy %>% data.table::as.data.table(),
            file = paste0(database.folder, "/all_split_reads_raw.rds"))
    
  }
  
  all_split_reads_raw_tidy %>% 
    distinct(junID) %>%
    nrow()
  
  #######################################################################
  ## 1. Discard all split reads shorter than 25 bp
  #######################################################################
  
  all_split_reads_raw_tidy %>% head()
  all_split_reads_raw_tidy %>% nrow()
  
  all_split_reads_raw_tidy <- all_split_reads_raw_tidy %>%
    filter(width >= 25)
  
  all_split_reads_raw_tidy %>% head()
  all_split_reads_raw_tidy %>% nrow()
  
  #######################################################################
  ## 2. Remove split reads located in unplaced sequences in the chromosomes
  #######################################################################
  
  print(paste0(Sys.time(), " - removing unplaced sequences genome..."))
  
  all_split_reads_raw_tidy_gr <- all_split_reads_raw_tidy %>%
    filter(!(strand=="?")) %>%
    GenomicRanges::GRanges() %>%
    diffloop::rmchr()
  rm(all_split_reads_raw_tidy)
  gc()
  
  all_split_reads_raw_tidy_gr %>% length()
  all_split_reads_raw_tidy_gr <- GenomeInfoDb::keepSeqlevels(x = all_split_reads_raw_tidy_gr,
                                                             value = intersect(all_split_reads_raw_tidy_gr %>% 
                                                                                 GenomeInfoDb::seqnames() %>% 
                                                                                 levels(), 
                                                                               c( "1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                                                                  "2", "20", "21", "22", "3", "4", "5", "6", "7", "8", "9", "X", "Y")), 
                                                             pruning.mode = "tidy")
  ## Data check
  all_split_reads_raw_tidy_gr %>% head()
  all_split_reads_raw_tidy_gr %>%
    as_tibble() %>%
    distinct(junID, .keep_all = T) %>%
    nrow() %>%
    print()
  
  #######################################################
  ## 3. Remove split reads overlapping the ENCODE backlist
  #######################################################
  
  print(paste0(Sys.time(), " - removing blacklist sequences..."))
  
  blacklist_path <- paste0(dependencies_folder, "/hg38-blacklist.v2.bed")
  all_split_reads_raw_tidy_gr <- remove_encode_blacklist_regions(GRdata = all_split_reads_raw_tidy_gr,
                                                                 blacklist_path = blacklist_path)
  
  ## Data check
  all_split_reads_raw_tidy_gr %>%
    as_tibble() %>%
    distinct(junID, .keep_all = T) %>%
    nrow() %>%
    print()
  
  
  
  #######################################################
  ## 4. Anotate using 'dasper'
  #######################################################
  
  print(paste0(Sys.time(), " - annotating dasper..."))
  
  if ( file.exists(paste0(dependencies_folder, "/Homo_sapiens.GRCh38.", gtf.version, ".chr.gtf")) ) {
    gtf_path <- paste0(dependencies_folder, "/Homo_sapiens.GRCh38.", gtf.version, ".chr.gtf")
  } else {
    gtf_path <- paste0(dependencies_folder, "/Homo_sapiens.GRCh38.", gtf.version, ".gtf")
  }
  edb <- ensembldb::ensDbFromGtf(gtf_path, outfile = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
  edb <- ensembldb::EnsDb(x = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
  
  all_split_reads_details_w_symbol <- dasper::junction_annot(junctions = all_split_reads_raw_tidy_gr %>% GRanges(), 
                                                             ref = edb,
                                                             ref_cols = c("gene_id", "gene_name", "symbol", "tx_id"), 
                                                             ref_cols_to_merge = c("gene_id", "gene_name", "tx_id"))
  all_split_reads_details_w_symbol %>% head
  rm(all_split_reads_raw_tidy_gr)
  rm(edb)
  gc()
  
  
  ################################################################################
  ## 5. Discard all junctions that are not annotated, novel donor or novel acceptor
  ################################################################################
  
  print(paste0(Sys.time(), " - removing split reads not classified as 'annotated', 'novel_donor' or 'novel_acceptor'..."))
  
  ## Subset columns
  all_split_reads_details_w_symbol <- all_split_reads_details_w_symbol[,c("junID", "gene_id_junction", "in_ref", "type", "tx_id_junction")]
  
  ## Only use annotated introns, novel donor and novel acceptor junctions
  all_split_reads_details_w_symbol <- all_split_reads_details_w_symbol[(elementMetadata(all_split_reads_details_w_symbol)[,"type"] %in% c("annotated", "novel_donor", "novel_acceptor"))]
  
  
  all_split_reads_details_w_symbol %>% length()
  
  
  ############################################
  ## 6. Make sure none of the split reads are shorter than 25bp
  ############################################
  
  print(paste0(Sys.time(), " - removing split reads shorter than 25bp..."))
  
  all_split_reads_details_w_symbol <- all_split_reads_details_w_symbol[all_split_reads_details_w_symbol %>% width() >= 25,]
  
  
  

  #saveRDS(all_split_reads_details_w_symbol,
  #        file = paste0(database.folder, "/all_split_reads_recount3_qc.rds"))
  
  
  
  
  ############################################################################
  ## 7. Discard all introns assigned to multiple genes (i.e. ambiguous introns)
  ############################################################################
  
  #if ( !exists("all_split_reads_details_w_symbol") ) {
  #  folder_path <- paste0(getwd(), "/database/v", gtf.version, "/")
  #  all_split_reads_details_w_symbol <- readRDS(file = paste0(folder_path, "/all_split_reads_recount3_", gtf.version, ".rds"))
  #}
  
  print(paste0(Sys.time(), " - discarding ambiguous introns..."))
  
  all_split_reads_details_w_symbol <- all_split_reads_details_w_symbol %>%
    as_tibble() %>%
    distinct(junID, .keep_all = T) %>% 
    rowwise() %>%
    mutate(ambiguous = ifelse(gene_id_junction %>% unlist() %>% length() > 1, T, F))
  
  ambiguous_introns <- all_split_reads_details_w_symbol %>%
    dplyr::filter(ambiguous == T)
  
  ambiguous_introns %>%
    distinct(junID, .keep_all = T) %>%
    dplyr::count(type) %>% 
    print()
  
  ambiguous_introns %>% 
    distinct(junID, .keep_all = T) %>%
    print()
  
  saveRDS(object = ambiguous_introns,
          file = paste0(database.folder, "/all_ambiguous_jxn.rds"))
  
  all_split_reads_details_w_symbol <- all_split_reads_details_w_symbol %>%
    dplyr::filter(ambiguous == F) %>%
    dplyr::select(-ambiguous)
  
  all_split_reads_details_w_symbol %>% nrow()
  
  saveRDS(object = all_split_reads_details_w_symbol %>% dplyr::rename(gene_id = gene_id_junction),
          file = paste0(database.folder, "/all_split_reads_qc_level1.rds"))
  
  
  ## FREE UP SOME MEMORY
  rm(edb)
  rm(ambiguous_introns)
  rm(all_split_reads_raw_tidy_gr)
  rm(all_split_reads_details_w_symbol)
  rm(all_split_reads_details_w_symbol_reduced)
  rm(all_split_reads_details_w_symbol_reduced_keep_gr)
  rm(all_split_reads_details_w_symbol_reduced_keep)
  rm(all_split_reads_details_w_symbol_reduced_discard)
  gc()
}


#' Title
#' Separates the quality-controlled split-read data by sample cluster for a given recount3 project
#' (e.g. by case/control cluster, by tissue, etc. )
#' @param recount3.project.IDs array with the recount3 projects to download
#' e.g. "SRP009615"
#' @param project.name name given locally to the recount3 project.
#' A given recount3 project can be separated in multiple independent projects in recount3.
#' For instance, GTEx is stored in recount3 in multiple independent projects ID (e.g. BRAIN, SKIN, BLOOD, etc), 
#' but all of them belong to GTEx. Hence, it is useulf to have a folder named "project.name" that will contain
#' multiple subfolders named as the elements contained within the object 'recount3.project.IDs'
#' @param gtf.version Ensembl version (tested using Ensembl v105)
#' e.g. "105"
#' @param data.source source of the project in recount3. It can be:
#' - "data_sources/sra"  
#' - "data_sources/gtex" 
#' - "data_sources/tcga"
#'
#' @return
#' @export
#'
#' @examples
tidy_sample_cluster_data <- function(recount3.project.IDs, 
                                     project.name,
                                     gtf.version,
                                     data.source,
                                     database.folder,
                                     results.folder) {
  
  
  message(Sys.time()," loading tidy split reads ID from ", project.name)
  
  all_split_reads_details_w_symbol_reduced_keep_gr <- 
    readRDS(file = paste0(database.folder, "/all_split_reads_qc_level1.rds"))
  all_split_reads_details_w_symbol_reduced_keep_gr %>%
    nrow()
  gc()
  
  # ## Generate the raw file
  for (project_id in recount3.project.IDs) {
    
    ## If this is GTEx, the recount3.project.IDs object can contain up to 54 different values
    ## e.g. KIDNEY, BRAIN, BLOOD, etc
    # project_id <- recount3.project.IDs[1]
    
    local_folder_results <- paste0(results.folder, "/", project_id, "/base_data/")
    dir.create(file.path(local_folder_results), recursive = TRUE, showWarnings = T)
    print(paste0(Sys.time(), " - getting junction data from recount3 - '", project_id, "' ..."))
    
    #############################################################################
    
    ## THERE ARE OTHER RECOMMENDED WAYS OF DOWNLOADING DATA FROM RECOUNT3
    ## HOWEVER, THIS ALTERNATIVE METHOD IS FOLLOWED IN ORDER TO REDUCE MEMORY USAGE
    ## PARTICULARLY USEFUL WITH LARGE DATASETS SUCH AS GTEX
    ## FOR OTHER RECOMMENDED METHODS, SEE (https://bioconductor.org/packages/release/bioc/manuals/recount3/man/recount3.pdf) - ACCESSED 08/07/2023
    
    #############################################################################
    
    bfc <- recount3::recount3_cache()
    recount3_url <- getOption("recount3_url", "http://duffel.rail.bio/recount3")
    verbose <- getOption("recount3_verbose", TRUE)
    
    message(Sys.time()," loading metadata info.")
    metadata.info <- recount3::read_metadata(recount3::file_retrieve(
      url = recount3::locate_url(
        project = project_id,
        project_home = data.source,
        organism = "human",
        annotation = "gencode_v29",
        type = "metadata",
        recount3_url = recount3_url
      ),
      bfc = bfc,
      verbose = verbose
    ))
    
    jxn_files <- recount3::locate_url(
      project = project_id,
      project_home = data.source,
      organism = "human",
      annotation = "gencode_v29",
      #jxn_format = c("UNIQUE"),
      type = "jxn",
      recount3_url = recount3_url
    )
    
    feature_info <- 
      readRDS(file = paste0(local_folder_results, "/all_split_reads_raw.rds")) %>%
      pull(junID)
    gc()
    
    feature_info %>% length()
    
    
    if (verbose) {
      message(
        Sys.time(),
        " matching exon-exon junction counts with the metadata."
      )
    }
    
    ## The samples in the MM jxn table are not in the same order as the
    ## metadata!
    jxn_rail <- read.delim(recount3::file_retrieve(
      url = jxn_files[grep("\\.ID\\.gz$", jxn_files)],
      bfc = bfc,
      verbose = verbose
    ))
    m <- match(metadata.info$rail_id, jxn_rail$rail_id)
    stopifnot(
      "Metadata rail_id and exon-exon junctions rail_id are not matching." =
        !all(is.na(m))
    )
    
    
    saveRDS(object = metadata.info, 
            file = paste0(local_folder_results, "/", project_id, "_samples_raw_metadata.rds"))
    
    #################################
    ## GET SAMPLE CLUSTERS
    #################################
    
    # metadata.info <- readRDS(file = paste0(folder_results, "/", project_id, "_samples_raw_metadata.rds"))
    metadata_tidy <- separate_clusters(project.metadata = metadata.info, data.source)
    metadata_tidy$all_mapped_reads %>% min
    metadata_tidy %>% head()

    #################################
    
    set.seed(12)
    ## From the samples obtained, only keep the ones matching similar RIN numbers
    m.out <- MatchIt::matchit(cluster ~ rin_score,
                              data = metadata_tidy,
                              distance = metadata_tidy$rin_score,
                              method = "nearest",
                              caliper = c(rin_score = 0),
                              std.caliper = F)
    metadata_tidy_filter <- MatchIt::match.data(m.out) %>%
      dplyr::select(-c("distance", "weights", "subclass"))
    metadata_tidy_filter %>% 
      distinct(external_id, .keep_all = T) %>%
      count(cluster)
    metadata_tidy_filter %>% 
      distinct(external_id, .keep_all = T) %>%
      as_tibble()
    
    
    saveRDS(object = metadata_tidy_filter, 
            file = paste0(local_folder_results, "/", project_id, "_samples_metadata.rds"))
    
    metadata_tidy_filter$all_mapped_reads %>% min
    metadata_tidy_filter %>% as_tibble()
    metadata_tidy_filter$external_id %>% sort()
    
    #################################
    ## GET SPLIT READS AND COUNTS  
    #################################
    
    
    
    clusters_ID <- metadata_tidy$cluster %>% unique()
    clusters_used <- NULL
    gc()
    
    for (cluster_id in clusters_ID) {
      
      # cluster_id <- clusters_ID[2]
      print(paste0(Sys.time(), " - filtering junction data by cluster - '", cluster_id, "' tissue"))
      
      message(Sys.time()," loading count matrix.")
      counts <- Matrix::readMM(recount3::file_retrieve(
        url = jxn_files[grep("\\.MM\\.gz$", jxn_files)],
        bfc = bfc,
        verbose = verbose
      ))
      counts %>% nrow()
      #counts <- counts[, m, drop = FALSE]
      message(Sys.time()," ordering count matrix.")
      colnames(counts) <- metadata_tidy$external_id[m]
      rownames(counts) <- feature_info
      gc()
      
      colnames(counts) %>% length()
      
      ################
      ## Get clusters
      ################
      
      cluster_samples <- metadata_tidy_filter %>% 
        filter(cluster == cluster_id) %>%
        distinct(external_id) %>% 
        pull()
      
      saveRDS(object = cluster_samples, 
              file = paste0(local_folder_results, "/", project_id, "_", cluster_id, "_samples_used.rds"))
      clusters_used <- c(clusters_used, cluster_id)
      
      
      ################
      ## Get counts
      ################
      
      print(paste0(Sys.time(), " - filtering split read counts matrix for the current matrix."))
      
      local_counts <- counts[,(colnames(counts) %in% cluster_samples)]
      local_counts %>% nrow()
      
      local_counts <- local_counts[(rownames(local_counts) %in% 
                                      all_split_reads_details_w_symbol_reduced_keep_gr$junID),]
      local_counts %>% nrow()
      
      print(paste0(Sys.time(), " - converting the split read counts matrix into tibble."))
      
      local_counts %>% nrow()
      local_counts %>% head()
      
      
      
      local_counts <- local_counts %>% as.matrix()
      local_counts <- local_counts[rowSums(local_counts) > 1, ]
      local_counts %>% nrow()
      
      #local_counts_linux <- read.csv(file = "~/PROJECTS/recount3-database-project/database/SRP100948/105/file.csv")
      #local_counts_linux %>% head
      #setdiff(local_counts_linux$junID,local_counts$junID)
      #local_counts %>% nrow()
      #gc()
      
      #local_counts["chr1:43422633-43422820:+", ]
      
      local_counts <- local_counts %>% as_tibble(rownames = "junID")
      # local_counts %>%
      #   filter(junID == "chr1:43422633-43422820:+") %>%
      #   as.data.frame()
      
      print(paste0(Sys.time(), " - saving data."))
      print(object.size(local_counts), units = "GB")
      
      stopifnot(
        "Still there are split reads with less than 2 supportive reads" =
          local_counts %>% 
          mutate(sumCounts = rowSums(select(., !contains("junID")))) %>%
          filter(sumCounts <= 1) %>% 
          nrow() == 0
      )
      
      saveRDS(object = local_counts,
              file = paste0(local_folder_results, "/", project_id, "_", cluster_id, "_split_read_counts.rds"))
      gc()
      
      #######################
      ## Get split reads ID
      #######################
      
      print(paste0(Sys.time(), " - getting split read IDs"))
      all_split_reads <- local_counts %>%
        dplyr::select(junID) %>%
        data.table::as.data.table() %>%
        inner_join(y = all_split_reads_details_w_symbol_reduced_keep_gr,
                   by = "junID")
      
      #######################
      ## qc and save results
      #######################
      
      if (any(all_split_reads$width < 25) |
          any(str_detect(string = all_split_reads$chr, pattern = "random")) |
          any(str_detect(string = str_to_lower(all_split_reads$chr), pattern = "u")) |
          any(!(all_split_reads$type %in% c("annotated", "novel_acceptor", "novel_donor")))) {
        print("ERROR! The split reads do not meet the minimum filtering criteria.")
        break;
      }
      
      ## Check how much memory this object uses
      print(object.size(all_split_reads %>% data.table::as.data.table()), units = "GB")
      
      print(paste0(Sys.time(), " - saving data."))
      ## Save the object
      saveRDS(object = all_split_reads %>% data.table::as.data.table(),
              file = paste0(local_folder_results, "/", project_id, "_", cluster_id, "_all_split_reads.rds"))
      
      
      ## Free up some memory
      rm(all_split_reads)
      rm(local_counts)
      rm(cluster_samples)
      rm(counts)
      gc()
      
    } 
    
    if ( !is.null(clusters_used) ) {
      saveRDS(object = clusters_used, 
              file = paste0(local_folder_results, "/", project_id, "_clusters_used.rds"))
    }
    
    #rm(rse)
    rm(metadata.info)
    rm(clusters_used)
    rm(clusters_ID)
    
    rm(feature_info)
    rm(jxn_files)
    rm(jxn_rail)
    rm(m)
    gc()
    
  }
  
}


#' Title
#' Each recount3 project has a different data source, hence the metadata is structured differently
#' This function extracts the metadata (clusters, RIN, etc) from each recount3 project depending on its data source
#' @param project.metadata raw metadata object as provided by recount3 (https://rdrr.io/bioc/recount3/man/read_metadata.html)
#' @param data.source source of the project in recount3. It can be:
#' - "data_sources/sra"  
#' - "data_sources/gtex" 
#' - "data_sources/tcga"
#'
#' @return
#' an standardised metadata object 
#' @export
#'
#' @examples
separate_clusters <- function(project.metadata,
                              data.source) {
  
  
  if (data.source == "data_sources/gtex") {
    
    project_metadata_tidy <- project.metadata %>%
      as_tibble() %>%
      filter(gtex.smafrze != "EXCLUDE",
             gtex.smrin >= 6.0) %>%
      dplyr::rename(cluster = gtex.gtex.smtsd,
                    rin = gtex.smrin %>% as.double())
    
    
    ## TODO same standard column names as with data_sources/sra projects
    project_metadata_tidy <- data.frame(age = project_metadata_tidy$gtex.age %>% as.character(),
                                        rin = project_metadata_tidy$gtex.smrin %>% as.character(),
                                        gender = project_metadata_tidy$gtex.sex %>% as.character(),
                                        cluster = project_metadata_tidy$gtex.smtsd,
                                        smnabtcht = project_metadata_tidy$gtex.smnabtcht,
                                        sample_id = project_metadata_tidy %>% distinct(gtex.sampid) %>% nrow(),
                                        smafrze = project_metadata_tidy$gtex.smafrze,
                                        avg_read_length = project_metadata_tidy$recount_seq_qc.avg_len,
                                        mapped_read_count = project_metadata_tidy$recount_qc.star.all_mapped_reads,
                                        SRA_project = project_metadata_tidy$recount_project.project) 
    
  } else if (data.source == "data_sources/sra") {
    
    project_metadata_tidy <- project.metadata %>%
      dplyr::select(external_id, 
                    sra.experiment_title, 
                    sra.sample_attributes, 
                    all_mapped_reads = recount_qc.star.all_mapped_reads) %>%
      mutate(rn = row_number()) %>%
      separate_rows(sra.sample_attributes, sep = "\\|\\s*") %>%
      separate(sra.sample_attributes, into = c('col1', 'col2'), sep = ";;") %>% 
      pivot_wider(names_from = col1, values_from = col2) %>% 
      dplyr::select(-rn) %>%
      as.data.frame() %>%
      mutate(cluster = ifelse(test = str_detect(sra.experiment_title, pattern="AD"),
                              yes = "AD",
                              no = "control"),
             cluster = cluster %>%as.factor(),
             rin_score = rin_score %>% as.double()) 
  } else {
    ## TODO "data_sources/tcga"
  }
  
  
  return(project_metadata_tidy)
  
}

#' Title
#' Function to pair novel junctions with annotated introns across the samples of each sample cluster 
#' (i.e. case/control, tissue, etc)
#' @param recount3.project.IDs array with the recount3 projects to download
#' e.g. "SRP009615" 
#' @param project.name name given locally to the recount3 project.
#' A given recount3 project can be separated in multiple independent projects in recount3.
#' For instance, GTEx is stored in recount3 in multiple independent projects ID (e.g. BRAIN, SKIN, BLOOD, etc), 
#' but all of them belong to GTEx. Hence, it is useulf to have a folder named "project.name" that will contain
#' multiple subfolders named as the elements contained within the object 'recount3.project.IDs'
#' @param gtf.version Ensembl version (tested using Ensembl v105)
#' e.g. "105"
#'
#' @return
#' @export
#'
#' @examples
junction_pairing <- function(recount3.project.IDs, 
                             project.name,
                             gtf.version,
                             database.folder,
                             results.folder) {
  
  
  
  for (project_id in recount3.project.IDs) {
    
    # project_id <- recount3.project.IDs[1]
    
    folder_root <- paste0(results.folder, "/", project_id)
    folder_base_data <- paste0(folder_root, "/base_data/")
    
    print(paste0(Sys.time(), " - getting data from '", project_id, "' project..."))
    
    ## Load clusters
    metadata.info <- readRDS(file = paste0(folder_base_data, "/", project_id, "_samples_metadata.rds"))
    clusters_ID <- readRDS(file = paste0(folder_base_data, "/", project_id, "_clusters_used.rds"))
    
    
    for (cluster_id in clusters_ID) {
      
      # cluster_id <- clusters_ID[1]
      # cluster_id <- clusters_ID[2]
      
      print(paste0(Sys.time(), " - loading '", cluster_id, "' source data ..."))
      
      ############################################
      ## LOAD DATA FOR THE CURRENT PROJECT
      ############################################
      
      ## Load samples
      samples_used <- readRDS(file = paste0(folder_base_data, "/", 
                                            project_id, "_", cluster_id, "_samples_used.rds"))
      
      ## Load split read data
      all_split_reads_details <- readRDS(file = paste0(folder_base_data, "/", project_id, "_", cluster_id, 
                                                       "_all_split_reads.rds")) %>% as_tibble()
      
      ## Load split read counts
      split_read_counts <- readRDS(file = paste0(folder_base_data, "/", project_id, "_", cluster_id, "_",
                                                 "split_read_counts.rds")) %>% as_tibble()
      
      if ( !identical(all_split_reads_details$junID, split_read_counts$junID) ) {
        print("ERROR! The number of junctions considered is not correct.")
        break;
      }
      
      ## TODO ONLY CONSIDER NOVEL JUNCTIONS AND INTRONS WITH > 1 SUPPORTING READS
      
      ############################################
      ## DISTANCES SUITE OF FUNCTIONS
      ############################################
      
      folder_pairing_results <- paste0(folder_root, "/junction_pairing/", cluster_id, "/")
      dir.create(file.path(folder_pairing_results), recursive = TRUE, showWarnings = T)
      
      get_distances(cluster = cluster_id,
                    samples = samples_used,
                    split_read_counts = split_read_counts,
                    all_split_reads_details = all_split_reads_details,
                    folder_name = folder_pairing_results)
      gc()
      
      
      extract_distances(cluster = cluster_id,
                        samples = samples_used,
                        folder_name = folder_pairing_results)
      gc()
      
      
      get_never_misspliced(cluster = cluster_id,
                           samples = samples_used,
                           split_read_counts = split_read_counts,
                           all_split_reads_details = all_split_reads_details,
                           folder_name = folder_pairing_results)
      
      rm(all_split_reads_details)
      rm(split_read_counts)
      gc()
      
      #}
    }
  }
}



#' Title
#' Removes the ambiguous junctions and prepares the data prior generation of the SQL database
#' @param recount3.project.IDs 
#' @param project.name 
#' @param gtf.version 
#' @param all.clusters 
#'
#' @return
#' @export
#'
#' @examples
tidy_data_pior_sql <- function (recount3.project.IDs,
                                project.name,
                                gtf.version,
                                all.clusters = NULL,
                                database.folder,
                                results.folder) {
  
  
  print(paste0(Sys.time(), " - loading base data ..."))
  

  
  ## Load base recount3 object containing the split reads passing the QC criteria
  all_split_reads_details_qc_level1 <- readRDS(file = paste0(database.folder, "/all_split_reads_qc_level1.rds")) %>%
    as_tibble()
  all_split_reads_details_qc_level1 %>% nrow()
  
  ############################################
  ## Discard all junctions from 
  ## EXCLUDE ME samples
  ############################################
  
  all_split_reads_details_qc_level2 <- readRDS(file = paste0(database.folder, "/all_split_reads_qc_level2.rds"))  %>%
    as_tibble()
  
  all_split_reads_details_qc_level2 %>% nrow()
  all_split_reads_details_qc_level2 %>% head()
  
  ## This should be zero
  if ( setdiff(all_split_reads_details_qc_level2$junID, 
               all_split_reads_details_qc_level1$junID) %>% 
       length() > 0) {
    print("ERROR! Some of the annotated split reads are not found within the base tidied recount3 object.")
    break;
  }
  
  # ## These are the number of split reads from the samples excluded
  setdiff(all_split_reads_details_qc_level1$junID, 
          all_split_reads_details_qc_level2$junID) %>% unique %>% length()
  
   
  
  ############################################
  ## QC
  ############################################ 
  
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = all_split_reads_details_qc_level1$junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    all_split_reads_details_qc_level1[ind, "junID"] <- 
      str_replace(string = all_split_reads_details_qc_level1[ind, "junID"]$junID, 
                  pattern = "\\*", 
                  replacement = all_split_reads_details_qc_level1[ind, "strand"]$strand %>% as.character() )
    if (any(str_detect(all_split_reads_details_qc_level1$junID, pattern = "\\*"))) {
      print("ERROR!")
      break;
    }
  }
  
  
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = all_split_reads_details_qc_level2$junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    all_split_reads_details_qc_level2[ind, "junID"] <- str_replace(string = all_split_reads_details_qc_level2[ind, "junID"]$junID, 
                                                                     pattern = "\\*", 
                                                                     replacement = all_split_reads_details_qc_level2[ind, "strand"]$strand %>% as.character() )
    if (any(str_detect(all_split_reads_details_qc_level2$junID, pattern = "\\*")) ) {
      print("ERROR!")
      break;
    }
  }
  
  
  ##########################################
  ## Load all distances pairings
  ##########################################
  
  df_all_jxn_pairings <- readRDS(file = paste0(database.folder, "/all_jxn_pairings.rds"))
  
  ## QC
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = df_all_jxn_pairings$ref_junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    df_all_jxn_pairings[ind, "ref_junID"] <- str_replace(string = df_all_jxn_pairings[ind, "ref_junID"]$ref_junID, 
                                                                   pattern = "\\*", 
                                                                   replacement = df_all_jxn_pairings[ind, "ref_strand"]$ref_strand %>% as.character())
    if( any(str_detect(df_all_jxn_pairings$ref_junID, pattern = "\\*")) ) {
      print("ERROR!")
      break;
    }
  }
  ## Remove potential * in the junID of the novel junctions
  ind <- which(str_detect(string = df_all_jxn_pairings$novel_junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    df_all_jxn_pairings[ind, "novel_junID"] <- str_replace(string = df_all_jxn_pairings[ind, "novel_junID"]$novel_junID, 
                                                                     pattern = "\\*", 
                                                                     replacement = df_all_jxn_pairings[ind, "novel_strand"]$novel_strand  %>% as.character())
    if (any(str_detect(df_all_jxn_pairings$novel_junID, pattern = "\\*")) ) {
      print("ERROR!")
      break;
    }
  }
  
  
  ##########################################
  ## Get never mis-spliced
  #########################################
  
  df_never_misspliced <- get_intron_never_misspliced(recount3.project.IDs = recount3.project.IDs,
                                                     all.clusters = all.clusters,
                                                     project.name = project.name,
                                                     database.folder,
                                                     results.folder)
  
  if ( any(str_detect(df_never_misspliced$ref_junID, pattern = "\\*")) ) {
    
    print("ERROR! Some junctions still contain a '*' in their IDs")
    break;
  }
  
  ## Remove the introns paired with novel junctions (i.e. mis-spliced)
  df_never_misspliced_tidy <- df_never_misspliced %>%
    dplyr::filter(!(ref_junID %in% df_all_jxn_pairings$ref_junID)) %>%
    as_tibble()
  
  df_never_misspliced_tidy %>% distinct(ref_junID) %>% as_tibble()
  
  
  ############################################
  ## GET all not paired
  ############################################ 
  
  # df_all_novel_raw_tidy
  df_not_paired <- all_split_reads_details_qc_level2 %>%
    data.table::as.data.table() %>%
    dplyr::filter(!(junID %in% c(df_all_jxn_pairings$ref_junID,
                                 df_all_jxn_pairings$novel_junID)))
  
  ## These are all the non-paired, including the never mis-spliced. Thus, GTEx 'introverse' data: [768,646 - 38,521 = 730125]
  df_not_paired %>% 
    distinct(junID)
  df_not_paired %>%
    distinct(junID, .keep_all = T) %>%
    dplyr::count(type)
  
  
  if (any(str_detect(df_never_misspliced_tidy$ref_junID, pattern = "\\*")) |
      any(str_detect(df_not_paired$junID, pattern = "\\*"))) {
    print("ERROR!")
  }
  
  ## All never mis-spliced should be categorised as not paired.
  ## Thus, this should be zero
  setdiff(df_never_misspliced_tidy$ref_junID, 
          df_not_paired %>% dplyr::filter(type == "annotated") %>% pull(junID))
  
  df_not_paired_tidy <- df_not_paired %>%
    dplyr::filter(!(junID %in% df_never_misspliced_tidy$ref_junID)) %>%
    distinct(junID, .keep_all = T) %>%
    as_tibble()
  
  df_not_paired_tidy %>%
    distinct(junID, .keep_all = T) %>%
    dplyr::count(type)
  df_not_paired_tidy %>% distinct(junID)
  
  
  ## This should be zero
  if ( intersect(df_not_paired_tidy$junID, 
                 df_all_jxn_pairings$novel_junID) %>% length() > 0 ) {
    print("ERROR!")
  }
  
  
  df_all_jxn_pairings %>% distinct(novel_junID) %>% nrow() +
    df_all_jxn_pairings %>% distinct(ref_junID) %>% nrow() +
    df_never_misspliced_tidy %>% distinct(ref_junID) %>% nrow()
  
  ##########################################
  ## Remove ambiguous junctions
  ##########################################
  
  
  ## All these should be zero
  if( intersect(df_not_paired_tidy$junID, df_all_jxn_pairings$novel_junID) %>% length() > 0 |
      intersect(df_not_paired_tidy$junID, df_all_jxn_pairings$ref_junID) %>% length() > 0 |
      intersect(df_not_paired_tidy$junID, df_never_misspliced_tidy$ref_junID) %>% length() > 0) {
    print("ERROR!")
  }
  
  
  ## 1. Obtain the ambiguous junctions
  
  df_ambiguous_novel <- df_all_jxn_pairings %>%
    dplyr::filter(!(novel_junID %in% df_not_paired_tidy$junID),
                  !(ref_junID %in% df_not_paired_tidy$junID),
                  !(ref_junID %in% df_never_misspliced_tidy$ref_junID)) %>%
    dplyr::group_by(novel_junID) %>%
    mutate(distances_sd = distance %>% sd()) %>%
    dplyr::filter(distances_sd > 0)
  
  
  df_ambiguous_novel %>% ungroup() %>% distinct(novel_junID)
  df_ambiguous_novel %>% ungroup() %>% distinct(ref_junID)
  
  
  ## 2. Remove ambiguous junctions
  
  df_all_jxn_pairings_tidy <- df_all_jxn_pairings %>%
    dplyr::filter(!(novel_junID %in% df_ambiguous_novel$novel_junID)) %>%
    data.table::as.data.table() %>%
    distinct(novel_junID, ref_junID, .keep_all = T) %>%
    mutate(ref_strand = ref_strand %>% as.character(),
           novel_strand = novel_strand %>% as.character()) 
  
  
  ## 3. Get ambiguous figures and stats
  
  df_all_jxn_pairings_tidy %>%
    dplyr::distinct(novel_junID)
  df_all_jxn_pairings_tidy %>%
    dplyr::distinct(ref_junID)
  
  
  (df_ambiguous_novel %>% 
      ungroup() %>%
      distinct(ref_junID) %>% nrow()) - (intersect(c(df_all_jxn_pairings_tidy$novel_junID,
                                                     df_all_jxn_pairings_tidy$ref_junID),
                                                   df_ambiguous_novel %>% 
                                                     ungroup() %>%
                                                     distinct(ref_junID) %>% pull) %>% length())
  
  (df_all_jxn_pairings_tidy %>%
      dplyr::distinct(novel_junID) %>% 
      nrow()) + 
    (df_all_jxn_pairings_tidy %>%
       dplyr::distinct(ref_junID) %>% 
       nrow()) + 
    (df_never_misspliced_tidy %>% 
       distinct(ref_junID) %>% 
       nrow())
  
  
  ##############################################################################
  ## SAVE FINAL OBJECT
  ##############################################################################
  
  ## 1. DISTANCES PAIRINGS
  
  if (any(str_detect(string = df_all_jxn_pairings_tidy$ref_junID, pattern = "\\*")) |
      any(str_detect(string = df_all_jxn_pairings_tidy$novel_junID, pattern = "\\*")) ) {
    print("ERROR! Some junctions still have a '*' within their IDs!")
  }
  df_all_jxn_pairings_tidy <- df_all_jxn_pairings_tidy %>%
    inner_join(y = all_split_reads_details_qc_level2 %>% 
                 dplyr::select(junID, gene_id, tx_id_junction),
               by = c("ref_junID" = "junID"))
  saveRDS(object = df_all_jxn_pairings_tidy,
          file = paste0(database.folder, "/all_jxn_correct_pairings.rds"))
  
  
  ## 2. NEVER MIS-SPLICED
  
  if (any(str_detect(string = df_never_misspliced_tidy$ref_junID, pattern = "\\*")) ) {
    print("ERROR! Some NEVER MIS-SPLICED junctions still have a '*' in their IDs!")
  }
  df_never_misspliced_tidy <- df_never_misspliced_tidy %>%
    inner_join(y = all_split_reads_details_qc_level2 %>% 
                 dplyr::select(junID, seqnames, start, end, width, strand, gene_id, tx_id_junction),
               by = c("ref_junID" = "junID"))
  
  saveRDS(object = df_never_misspliced_tidy,
          file = paste0(database.folder, "/all_jxn_never_misspliced.rds"))
  
  
  ## 3. AMBIGUOUS JUNCTIONS
  
  if (any(str_detect(string = df_ambiguous_novel$ref_junID, pattern = "\\*")) |
      any(str_detect(string = df_ambiguous_novel$novel_junID, pattern = "\\*")) ) {
    print("ERROR! Some junctions still have a * in their IDs!")
  }
  saveRDS(df_ambiguous_novel,
          file = paste0(database.folder, "/all_jxn_ambiguous_pairings.rds"))
}


#' Title
#' SQL database helper function to create the different tables sequentially
#' In case of removing different tables, useful to control which tables are removed 
#' (i.e. only child tables or also master tables) 
#' @param database.path 
#' @param recount3.project.IDs 
#' @param project.name 
#' @param gtf.version 
#' @param remove.all 
#'
#' @return
#' @export
#'
#' @examples
sql_database_generation <- function(database.path,
                                    recount3.project.IDs, 
                                    project.name,
                                    gtf.version,
                                    remove.all = NULL,
                                    database.folder,
                                    results.folder) {
  
  
  print(paste0(Sys.time(), " --> ", database.path, "..."))
  
  con <- dbConnect(RSQLite::SQLite(), database.path)
  tables <- DBI::dbListTables(conn = con)
  
  print("Database tables:")
  tables %>% print()
  
  if ( !is.null(remove.all) ) {
    
    remove_tables(database.path = database.path, all = remove.all)
    
  }
  
  
  
  
  if ( !any(tables == 'master') ) {
    create_metadata_table(database.path = database.path,
                          project.name = project.name,
                          gtf.version = gtf.version,
                          recount3.project.IDs = recount3.project.IDs,
                          results.folder = results_folder)
  }
  
  
  if ( ! any(tables %in% c('intron', 'novel', 'gene', 'transcript')) ) {
    create_master_tables(database.path = database.path,
                         project.name = project.name,
                         gtf.version = gtf.version,
                         database.folder,
                         results.folder)
  }
  
  print(paste0(Sys.time(), " - creating cluster tables ..."))

  create_cluster_tables(database.path = database.path,
                        gtf.version = gtf.version,
                        recount3.project.IDs = recount3.project.IDs,
                        project.name = project.name,
                        database.folder,
                        results.folder)
  
}


get_all_annotated_split_reads <- function(recount3.project.IDs,
                                          gtf.version,
                                          all.clusters = NULL,
                                          project.name,
                                          database.folder,
                                          results.folder) {
  
  
  
  #############################################
  ## These are all split reads from all tissues
  ## obtained from 'USE ME' SAMPLES and samples with more than 6 RIN
  #############################################
  
  ## The sample filter in IntroVerse corresponded to 
  
  all_split_reads_details_all_sample_clusters <- map_df(recount3.project.IDs, function(project_id) {
    
    # project_id <- recount3.project.IDs[1]
    # project_id <- "BONE_MARROW"
    
    folder_results_root <- paste0(results.folder, "/", project_id, "/")
    
    ## Get the metadata and clusters info
    metadata.info <- readRDS(file = paste0(folder_results_root, "/base_data/", 
                                           project_id, "_samples_metadata.rds"))
    
    if (is.null(all.clusters)) {
      all_clusters <-  readRDS(file = paste0(folder_results_root, "/base_data/", 
                                             project_id, "_clusters_used.rds"))
    }
    
    
    all_jxn_qc <- map_df(all_clusters, function(cluster) {
      
      # cluster <- all_clusters[1]
      print(paste0(Sys.time(), " - ", project_id, " loading '", cluster, "'  data ..."))
      
      if ( file.exists(paste0(folder_results_root, "/base_data/", 
                              project_id, "_", cluster, "_all_split_reads.rds")) ) {
        
        all_split_reads_details_105 <- readRDS(file = paste0(folder_results_root, "/base_data/", project_id, "_",
                                                             cluster, "_all_split_reads.rds"))
        
        all_split_reads_details_tidy <- all_split_reads_details_105 %>% 
          distinct(junID, .keep_all = T) %>% 
          as_tibble() %>%
          return()
        
      } else {
        return(NULL)
      }
      
    })
    
    if (all_jxn_qc %>% nrow() > 0 ) {
      
      all_jxn_qc %>%
        distinct(junID, .keep_all = T) %>% 
        return()
      
    } else {
      return(NULL)
    }
    
  })
  
  
  print(paste0(Sys.time(), " - saving 'all_annotated_split_reads' for the database!"))
  
  all_split_reads_details_all_sample_clusters <- all_split_reads_details_all_sample_clusters %>%
    distinct(junID, .keep_all = T)
  
  ## Save data
  saveRDS(object = all_split_reads_details_all_sample_clusters,
          file = paste0(database.folder, "/all_split_reads_qc_level2.rds") )
  
}


get_all_raw_distances_pairings <- function(recount3.project.IDs,
                                           gtf.version,
                                           all.clusters = NULL,
                                           project.name,
                                           database.folder,
                                           results.folder) {
  
  
  ## LOOP THROUGH PROJECTS
  df_all_distances_pairings_raw <- map_df(recount3.project.IDs, function(project_id) {
    
    # project_id <- recount3.project.IDs[1]
    # project_id <- "KIDNEY"
    
    print(paste0(Sys.time(), " --> Working with '", project_id, "' DataBase..."))
    folder_results_root <- paste0(results.folder, "/", project_id, "/")
    
    
    ## Read all clusters considered from the current tissue
    metadata.info <- readRDS(file = paste0(folder_results_root, "/base_data/", 
                                           project_id, "_samples_metadata.rds"))
    
    
    if ( is.null(all.clusters) ) {
      all.clusters <-  readRDS(file = paste0(folder_results_root, "/base_data/", 
                                             project_id, "_clusters_used.rds"))
    }
    
    map_df(all.clusters, function(cluster) {
      
      # cluster <- all.clusters[1]
      
      print(paste0(Sys.time(), " - ", project_id, " loading '", cluster, "'  data ..."))
      
      ## Load samples
      if ( file.exists(paste0(folder_results_root, "/base_data/", project_id, "_", cluster, "_samples_used.rds")) ) {
        
        samples <- readRDS(file = paste0(folder_results_root, "/base_data/", project_id, "_", cluster,  "_samples_used.rds"))
        
        if ( samples %>% length() > 0 ) {
          
          folder_cluster_pairings <- paste0(folder_results_root, "/junction_pairing/", cluster, "/")
          
          if ( !file.exists(paste0(folder_cluster_pairings, "/", cluster, "_raw_distances_tidy.rds")) ) {
            
            ## Obtain the distances across all samples
            df_all <- map_df(samples, function(sample) { 
              
              # sample <- samples[1]
              
              file_name <- paste0(folder_cluster_pairings, "/", cluster, "_", sample, "_distances.rds")
              
              
              if ( file.exists(file_name) ) {
                print(paste0(cluster, " - ", sample))
                df <- readRDS(file = file_name)
                
                return(df)
              } else {
                return(NULL)
              }
              
            })
            
            if ( nrow(df_all) > 0 ) {
              saveRDS(object = df_all %>%
                        distinct(novel_junID, ref_junID, .keep_all = T) %>%
                        mutate(tissue = cluster),
                      file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds"))
            }
            
          } else {
            print(paste0("File '", cluster, "_raw_distances_tidy.rds' already exists!"))
            df_all <- readRDS( file = paste0(folder_cluster_pairings, "/", cluster, "_raw_distances_tidy.rds") )
          }
          
          
          if ( nrow(df_all) > 0 ) {
            
            df_all %>%
              distinct(novel_junID, ref_junID, .keep_all = T) %>%
              mutate(project = project_id) %>%
              return()
            
          } else {
            return(NULL)
          }          
          # df_all2 <- readRDS(file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds"))
          
        } else {
          print(paste0("ERROR: no samples available for the tissue: ", project_id))
          return(NULL)
        }
        
      } else {
        return(NULL)
      }
    })  
  })  
  
  
  print(paste0(Sys.time(), " - saving 'df_all_distances_pairings' for the database!"))
  
  
  saveRDS(object = df_all_distances_pairings_raw %>% distinct(project) %>% pull(),
          file = paste0(results.folder, "/all_final_projects_used.rds"))
  
  saveRDS(object = df_all_distances_pairings_raw %>%
            distinct(novel_junID, ref_junID, .keep_all = T),
          file = paste0(database.folder,"/all_jxn_pairings.rds"))
  
}

get_intron_never_misspliced <- function (recount3.project.IDs,
                                         all.clusters = NULL,
                                         project.name,
                                         database.folder,
                                         results.folder) {
  
  
  df_never <- map_df(recount3.project.IDs, function(project_id) {
    
    # project_id <- recount3.project.IDs[1]
    
    print(paste0(Sys.time(), " --> Working with '", project_id, "' DataBase..."))
    base_folder <- paste0(results.folder, "/", project_id)
    
    if ( is.null(all.clusters) && file.exists(paste0(base_folder, "/base_data/", 
                                                     project_id, "_samples_metadata.rds")) ) {
      metadata.info <- readRDS(file = paste0(base_folder, "/base_data/", 
                                             project_id, "_samples_metadata.rds"))
      all_clusters <-  readRDS(file = paste0(base_folder, "/base_data/", 
                                             project_id, "_clusters_used.rds"))
      
    }
    
    
    map_df(all_clusters, function(cluster) { 
      
      # cluster <- all_clusters[1]
      
      print(paste0(Sys.time(), " --> ", cluster))
      if ( file.exists(paste0(base_folder, "/junction_pairing/", cluster, 
                              "/not-misspliced/", cluster, "_all_notmisspliced.rds")) ) {
        df_introns_never <- readRDS(file = paste0(base_folder, "/junction_pairing/", cluster, 
                                                  "/not-misspliced/", 
                                                  cluster, "_all_notmisspliced.rds")) %>% as_tibble()
        return(data.frame(ref_junID = df_introns_never$value))
      } else {
        return(NULL)
      }
      
    })
  })
  
  df_never %>%
    distinct(ref_junID) %>%
    return()
}


get_mean_coverage <- function(split_read_counts,
                              samples,
                              junIDs) {
  
  split_read_counts_intron <- split_read_counts %>%
    dplyr::filter(junID %in% junIDs) %>%
    dplyr::select(junID, all_of(samples %>% as.character())) 
  
  split_read_counts_intron[,"n_individuals"] <- (matrixStats::rowCounts(split_read_counts_intron[, -c(1)] > 0, na.rm = T)) 
  split_read_counts_intron <- split_read_counts_intron %>% as.data.frame()
  
  split_read_counts_intron[,"sum_counts"] <- Matrix::rowSums(split_read_counts_intron[,-c(split_read_counts_intron %>% ncol(),1)], na.rm = T)
  split_read_counts_intron <- split_read_counts_intron %>% as.data.frame()
  
  split_read_counts_intron <- split_read_counts_intron[, c(1,(split_read_counts_intron %>% ncol() - 1),(split_read_counts_intron %>% ncol()))]
  
  
  if (any(split_read_counts_intron[, "n_individuals"] < 1)) {
    print("Error: some ref junctions do not present any read across any of the samples.")
    break;
  }
  
  split_read_counts_intron %>% return()
}


generate_transcript_biotype_percentage <- function(recount3.project.IDs,
                                                   project.name,
                                                   gtf.version,
                                                   database.folder,
                                                   results.folder) {
  
  
  #######################################
  ## GET THE TRANSCRIPT BIOTYPE
  #######################################
  
  print(paste0(Sys.time(), " - loading the human reference transcriptome ... "))
  
  ## Import HUMAN REFERENCE transcriptome
  homo_sapiens_v105 <- rtracklayer::import(con = paste0(dependencies_folder, "/Homo_sapiens.GRCh38.",gtf.version,".chr.gtf")) %>% 
    as_tibble()
  
  ## Get v105 transcripts
  transcripts_v105 <- homo_sapiens_v105 %>%
    filter(type == "transcript") %>% 
    dplyr::select(transcript_id, transcript_biotype, gene_id)
  
  transcripts_v105 %>% head()
  
  
  #######################################
  ## LOAD ALL SPLIT READS 
  #######################################
  
  print(paste0(Sys.time(), " - loading the recount3 split reads ... "))
  
  ## LOAD the all split reads from all recount3 GTEx projects
  all_split_reads_details_all_tissues <- readRDS(file = paste0(database.folder, "/all_split_reads_qc_level2.rds") )
  
  all_split_reads_details_all_tissues %>% head()
  all_split_reads_details_all_tissues %>% nrow()
  
  
  
  
  #######################################
  ## EXPLORE AND TIDY THE RESULT
  #######################################
  
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = all_split_reads_details_all_tissues$junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    all_split_reads_details_all_tissues[ind, "junID"] <- str_replace(string = all_split_reads_details_all_tissues[ind, "junID"]$junID, 
                                                                     pattern = "\\*", 
                                                                     replacement = all_split_reads_details_all_tissues[ind, "strand"]$strand %>% as.character() )
    any(str_detect(all_split_reads_details_all_tissues$junID, pattern = "\\*")) %>% print()
  }
  
  print(object.size(all_split_reads_details_all_tissues), units = "Gb")
  
  ## Merge datasets to add transcript biotype
  print(paste0(Sys.time(), " --> adding transcript biotype..."))
  
  ## Merging using data table structures saves time
  df_all_junctions <- all_split_reads_details_all_tissues %>% unnest(tx_id_junction) %>% data.table::as.data.table()
  transcripts_v105 <- transcripts_v105 %>% data.table::as.data.table()
  
  df_all_junctions %>% head()
  transcripts_v105 %>% head()
  
  df_all_junctions_tx <- df_all_junctions %>% 
    left_join(y = transcripts_v105,
              by = c("tx_id_junction" = "transcript_id"))
  
  print(object.size(df_all_junctions_tx), units = "Gb")
  
  
  #######################################
  ## CALCULATE THE TRANSCRIPT PERCENTAGE
  #######################################
  
  
  print(paste0(Sys.time(), " --> starting protein-coding percentage calculation ..."))
  print(paste0(Sys.time(), " --> ", df_all_junctions$junID %>% unique() %>% length(), " total number of junctions."))
  
  
  ## Only keep not ambiguous jxn (i.e. junctions belonging to multiple genes)
  junID_OK <- df_all_junctions %>% 
    dplyr::group_by(junID) %>% 
    distinct(gene_id, .keep_all = T) %>% 
    dplyr::count() %>% 
    filter(n == 1) %>%
    pull(junID)
  
  # junID_OK %>% unique() %>% length()
  
  ## Calculate the biotype percentage
  df_all_percentage <- df_all_junctions_tx %>% 
    filter(junID %in% junID_OK) %>%
    dplyr::group_by(junID, transcript_biotype) %>%
    distinct(tx_id_junction, .keep_all = T) %>% 
    summarise(n = n()) %>% 
    mutate(percent = (n / sum(n)) * 100) %>%
    ungroup()
  
  df_all_percentage %>% head()
  saveRDS(object = df_all_percentage,
          file = paste0(results.folder, "/all_split_reads_qc_level2_biotype.rds"))
  
  ## Only filter by the protein-coding biotype
  df_all_percentage_tidy_PC <- df_all_percentage %>% 
    dplyr::group_by(junID) %>%
    rowwise() %>%
    mutate(percent = ifelse (transcript_biotype == "protein_coding", percent, 0)) %>%
    ungroup() %>%
    dplyr::group_by(junID) %>%
    filter(percent == max(percent)) %>%
    dplyr::select(-transcript_biotype, -n) %>%
    distinct(junID, .keep_all = T) %>%
    ungroup()
  
  ## Only filter by the lncRNA biotype
  df_all_percentage_tidy_lncRNA <- df_all_percentage %>% 
    dplyr::group_by(junID) %>%
    rowwise() %>%
    mutate(percent = ifelse (transcript_biotype == "lncRNA", percent, 0)) %>%
    ungroup() %>%
    dplyr::group_by(junID) %>%
    filter(percent == max(percent)) %>%
    dplyr::select(-transcript_biotype, -n) %>%
    distinct(junID, .keep_all = T) %>%
    ungroup()
  
  
  df_all_percentage_tidy_merged <- merge(x = df_all_percentage_tidy_PC %>% dplyr::rename(percent_PC = percent),
                                         y = df_all_percentage_tidy_lncRNA %>% dplyr::rename(percent_lncRNA = percent),
                                         by = "junID")
  df_all_percentage_tidy_merged %>% nrow()
  
  if (df_all_percentage_tidy_merged %>% filter(percent_PC == 100) %>% distinct(percent_lncRNA) %>% pull() != 0) {
    print("ERROR! some only protein-coding introns have been also classified as lncRNAs!")
  }
  if (df_all_percentage_tidy_merged %>% filter(percent_lncRNA == 100) %>% distinct(percent_PC) %>% pull() != 0) {
    print("ERROR! some only lncRNA introns have been also classified as protein-coding!")
  }
  
  print(object.size(df_all_percentage_tidy_merged), units = "Gb")
  
  saveRDS(object = df_all_percentage_tidy_merged,
          file = paste0(results.folder, "/all_split_reads_qc_level2_PC_biotype.rds"))
  print(paste0(Sys.time(), " - files saved!"))
  
  ##########################################
  ## FREE UP SOME MEMORY 
  ##########################################
  
  rm(all_split_reads_details_all_tissues)
  rm(df_all_junctions)
  rm(df_all_junctions_tx)
  rm(homo_sapiens_v105)
  rm(transcripts_v105)
  gc()
  
}



calculate_tpm <- function(rse, ref_tidy) {
  
  
  # Remove anything after . in ensembl id
  rownames(rse) <- rownames(rse) %>% 
    str_remove("\\..*")
  
  
  # Convert to tpm, which is calculated by:
  # 1. Divide the read counts by the length of each gene in kilobases (i.e. RPK)
  # 2. Count up all the RPK values in a sample and divide this number by 1,000,000.
  # 3. Divide the RPK values by the “per million” scaling factor.
  
  
  srp_rpk <- 
    rse %>% 
    SummarizedExperiment::assay() %>%
    as_tibble(rownames = "gene") %>% 
    tidyr::pivot_longer( ## equivalent to gather
      cols = -c("gene"),
      names_to = "recount_id",
      values_to = "counts"
    ) %>% 
    dplyr::inner_join(
      ref_tidy %>% 
        as_tibble() %>% 
        dplyr::select(gene_id, width),
      by = c("gene" = "gene_id")
    ) %>% # 1. Divide the read counts by the length of each gene in kilobases (i.e. RPK)
    dplyr::mutate(
      rpk = counts/width
    ) 
  srp_rpk %>% head()
  
  
  tpm <- 
    srp_rpk %>% 
    # 2. Count up all the RPK values in a sample and divide this number by 1,000,000.
    dplyr::inner_join(
      srp_rpk %>% 
        dplyr::group_by(recount_id) %>% 
        dplyr::summarise(
          scaling_factor = sum(rpk)/1e6
        ),
      by = "recount_id"
    ) %>% # 3. Divide the RPK values by the “per million” scaling factor.
    dplyr::mutate(
      tpm = rpk/scaling_factor
    )   
  
  tpm <- tpm %>%
    dplyr::select(gene, recount_id, tpm) %>% 
    spread(key = recount_id, value = tpm)
  
  
  return(tpm)
}


generate_recount3_tpm <- function(recount3.project.IDs,
                                  gtf.version,
                                  project.name,
                                  data.source,
                                  database.folder,
                                  results.folder) {
  
  # Reference gtf
  ref <- rtracklayer::import(con = paste0(dependencies_folder, "/Homo_sapiens.GRCh38.",gtf.version,".chr.gtf"))
  ref <- ref %>% GenomeInfoDb::keepSeqlevels(c(1:22), pruning.mode = "coarse") 
  
  ## Per each gene, calculate the length of its coding sequence
  ref_tidy <- ref %>%
    as_tibble() %>% 
    filter(type == "gene") %>%
    distinct(gene_id, .keep_all = T)
  
  for (project_id in recount3.project.IDs) {
    
    # project_id <- recount3.project.IDs[1]
    # project_id <- recount3.project.IDs[8]
    
    
    ## 1. Get expression data from recount3 and transform raw counts
    rse <- recount3::create_rse_manual(
      project = project_id,
      project_home = data.source,
      organism = "human",
      annotation = "gencode_v29",
      type = "gene")
    
    SummarizedExperiment::assays(rse)$counts <- recount3::transform_counts(rse)
    
    recount_tpm <- calculate_tpm(rse = rse, ref_tidy)
    
    rm(rse)
    gc()
    
    
    
    ## 2. For each tissue within the current project, filter the RSE by its samples
    
    results_folder_local <- paste0(results.folder, "/", project_id, "/")
    
    if ( file.exists(paste0(results_folder_local, "/base_data/", project_id, "_samples_metadata.rds")) ) {
      
      metadata.info <- readRDS(file = paste0(results_folder_local, "/base_data/", project_id, "_samples_metadata.rds"))
      
      clusters_ID <- readRDS(file = paste0(results_folder_local, "/base_data/", project_id, "_clusters_used.rds")) %>% unique()
      
      for ( cluster_id in clusters_ID ) {
        
        # cluster_id <- clusters_ID[1]
        samples <- readRDS(file = paste0(results_folder_local, "/base_data/", project_id, "_",cluster_id,"_samples_used.rds")) %>% unique()
        
        
        cluster_samples <- readRDS(file = paste0(results_folder_local, "/base_data/",
                                                 project_id, "_", cluster_id, "_samples_used.rds"))
        
        ## Filter the object for the set of samples corresponding to the current cluster
        recount_tpm_local <- recount_tpm %>%
          dplyr::select(c("gene", all_of(cluster_samples)))
        
        ## Save results
        results_folder_local_tpm <- paste0(results_folder_local, "/tpm/")
        dir.create(file.path(results_folder_local_tpm), recursive = TRUE, showWarnings = T)
        saveRDS(object = recount_tpm_local,
                file = paste0(results_folder_local_tpm, project_id, "_", cluster_id, "_tpm.rds"))
        
        
        rm(recount_tpm_local)
        rm(cluster_samples)
        gc()
        
        #}
        
      }
    }
    
    print(paste0(Sys.time(), " - ", project_id, " finished!"))
    
    rm(folder_root)
    rm(clusters_ID)
    rm(recount_tpm)
    gc()
    
  }
  
}