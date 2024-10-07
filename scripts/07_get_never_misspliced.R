#' #' Get Never Mis-spliced Junctions
#' #'
#' #' This function filters junctions to keep those that are never mis-spliced,
#' #' discarding any junctions that overlap with mis-spliced junctions.
#' #'
#' #' @param cluster A character string indicating the cluster name.
#' #' @param samples A character vector of sample identifiers.
#' #' @param split.read.counts A data frame containing split read counts.
#' #' @param all.split.reads.details A data frame with details of all split reads.
#' #' @param folder.name A character string indicating the folder to save results.
#' #' @param replace A logical indicating whether to replace existing results.
#' #'
#' #' @return A list containing the never mis-spliced and mis-spliced junctions.
#' #' @export
#' #'
#' #' @examples
#' #' GetNeverMisspliced(cluster = "cluster1", samples = c("sample1", "sample2"),
#' #'                    split.read.counts = split_counts_df,
#' #'                    all.split.reads.details = all_details_df,
#' #'                    folder.name = "results", replace = TRUE)
#' GetNeverMisspliced <- function(cluster, 
#'                                samples,
#'                                split.read.counts,
#'                                all.split.reads.details,
#'                                folder.name,
#'                                replace) {
#'   
#'   logger::log_info(paste0(Sys.time(), " - Filtering junctions that are potentially not mis-spliced..."))
#'   
#'   # Check if required file exists
#'   distances_file <- file.path(folder.name, paste0(cluster, "_raw_distances_tidy.rds"))
#'   if (!file.exists(distances_file)) {
#'     stop(paste0("File: '", distances_file, "' doesn't exist."))
#'   }
#'   
#'   # Check for missing junction IDs
#'   if (any(!all.split.reads.details$junID %in% split.read.counts$junID)) {
#'     stop("ERROR: Some junctions in all.split.reads.details are missing in split.read.counts.")
#'   }
#'   
#'   # Load mis-spliced junctions
#'   df_all_misspliced <- readRDS(distances_file) %>%
#'     dplyr::distinct(ref_junID, .keep_all = TRUE)
#'   
#'   # Get all junctions not paired on the first round - these have the potential of being non-mispliced
#'   all_not_misspliced <- all.split.reads.details %>%
#'     dplyr::filter(!(junID %in% df_all_misspliced$ref_junID))
#'   
#'   # Check for duplicates
#'   if (any(duplicated(df_all_misspliced$ref_junID)) || 
#'       any(duplicated(all_not_misspliced$junID))) {
#'     stop("ERROR: Duplicates found in junctions.")
#'   }
#'   
#'   # Main processing loop for samples
#'   
#'   introns_not_misspliced <- character()
#'   
#'   if (replace || 
#'       (!file.exists(file.path(folder.name, paste0(cluster, "_all_notmisspliced.rds"))) && 
#'        !file.exists(file.path(folder.name, paste0(cluster, "_all_misspliced_not_paired.rds"))))) {
#'     
#'     for (sample in samples) {
#'       
#'       introns_misspliced <- character()
#'       
#'       # sample <- samples[3]
#'       if (!sample %in% colnames(split.read.counts)) {
#'         next
#'       }
#'       
#'       # Process sample
#'       split_read_counts_sample <- split.read.counts %>%
#'         dplyr::select(junID, all_of(sample))
#'       
#'       # Merge details with counts
#'       all_not_misspliced_details_sample <- all_not_misspliced %>%
#'         dplyr::inner_join(split_read_counts_sample, by = "junID") %>%
#'         dplyr::rename(counts = !!sym(sample))
#'       
#'       
#'       # Quality control check
#'       sample_check <- all_not_misspliced_details_sample[sample(1:nrow(all_not_misspliced_details_sample), 1), ]
#'       if (!identical(sample_check$counts, split_read_counts_sample[split_read_counts_sample$junID == sample_check$junID,2][[1]])) {
#'         stop("ERROR: QC failed!")
#'       }
#'       
#'       logger::log_info(paste0(cluster, " - processing sample '", sample, "'..."))
#'       
#'       # Filtering annotated, donor, and acceptor junctions
#'       all_annotated <- all_not_misspliced_details_sample %>%
#'         dplyr::filter(type == "annotated") %>%
#'         GenomicRanges::GRanges()
#'       
#'       all_donor <- all_not_misspliced_details_sample %>%
#'         dplyr::filter(type == "novel_donor") %>%
#'         GenomicRanges::GRanges()
#'       
#'       all_acceptor <- all_not_misspliced_details_sample %>%
#'         dplyr::filter(type == "novel_acceptor") %>%
#'         GenomicRanges::GRanges()
#'       
#'       
#'       PairNovelJunctionsIntrons <- function(novel, annotated, pair_type) {
#'         overlaps <- GenomicRanges::findOverlaps(query = novel, subject = annotated, ignore.strand = FALSE, type = pair_type)
#'         introns_misspliced <- annotated[subjectHits(overlaps),]$junID
#'         return(introns_misspliced)
#'       }
#'       
#'       # Pairing logic for donors
#'       if (length(all_donor) > 0) {
#'         introns_misspliced <- c(PairNovelJunctionsIntrons(novel = all_donor[all_donor@strand == "+"], 
#'                                                           annotated = all_annotated[all_annotated@strand == "+"], pair_type = "end"),
#'                                 
#'                                 PairNovelJunctionsIntrons(novel = all_donor[all_donor@strand == "-"], 
#'                                                           annotated = all_annotated[all_annotated@strand == "-"], pair_type = "start"))
#'       }
#'       logger::log_info(paste0("Processed sample '", sample, "': ", length(introns_misspliced), " mis-spliced introns at the donor splice site..."))
#'       
#'       # Pairing logic for acceptors
#'       if (length(all_acceptor) > 0) {
#'         introns_misspliced <- c(introns_misspliced, 
#'                                 PairNovelJunctionsIntrons(all_acceptor[all_acceptor@strand == "+"], all_annotated[all_annotated@strand == "+"], pair_type = "start" ),
#'                                 PairNovelJunctionsIntrons(all_acceptor[all_acceptor@strand == "-"], all_annotated[all_annotated@strand == "-"], pair_type = "end"))
#'       }
#'       logger::log_info(paste0("Processed sample '", sample, "': ", length(introns_misspliced), " mis-spliced introns at the donor & acceptor splice sites..."))
#'       
#'       # Collect non-mis-spliced junctions
#'       introns_not_misspliced <- unique(c(introns_not_misspliced, 
#'                                          all_annotated$junID[!all_annotated$junID %in% introns_misspliced]))
#'       logger::log_info(paste0("Processed sample '", sample, "': ", length(introns_not_misspliced), " not mis-spliced introns..."))
#'       
#'       
#'     }
#'     
#'     # Final checks and save results
#'     introns_not_misspliced <- setdiff(introns_not_misspliced, introns_misspliced)
#'     logger::log_info(paste0("Processed sample '", sample, "': ", length(introns_not_misspliced), " not mis-spliced introns..."))
#'     
#'     if (length(intersect(introns_not_misspliced, introns_misspliced)) > 0) {
#'       stop("ERROR: Overlapping not-mispliced and mis-spliced junctions found!")
#'     }
#'     
#'     logger::log_info(paste0(Sys.time(), " - ", length(introns_not_misspliced), " not mis-spliced junctions!"))
#'     folder_path <- file.path(folder.name, "not-misspliced")
#'     dir.create(folder_path, showWarnings = FALSE)
#'     
#'     saveRDS(introns_not_misspliced, file.path(folder_path, paste0(cluster, "_all_notmisspliced.rds")))
#'     saveRDS(junc_ignore, file.path(folder_path, paste0(cluster, "_all_misspliced_not_paired.rds")))
#'     
#'     logger::log_info(paste0(Sys.time(), " - Results saved!"))
#'   }
#' }
#' 






#' Title
#'
#' @param cluster 
#' @param samples 
#' @param split.read.counts 
#' @param all.split.reads.details 
#' @param folder.name 
#'
#' @return
#' @export
#'
#' @examples
GetNeverMisspliced <- function(cluster, 
                               samples,
                               split.read.counts,
                               all.split.reads.details,
                               folder.name,
                               replace) {
  
  #######################################
  ## GET ALL JUNCTIONS THAT WERE NOT 
  ## MIS-SPLICED IN THE 1ST JXN PAIRING
  #######################################
  
  logger::log_info(paste0(Sys.time(), " - filtering junctions that are potentially not mis-spliced..."))
  
  if ( file.exists(paste0(folder.name, "/", cluster, "_raw_distances_tidy.rds")) ) {
    
    # ## This should be zero
    if ( setdiff(all.split.reads.details$junID, split.read.counts$junID) %>% length() > 0) {
      logger::log_info("ERROR")
      break;
    }
    
    all.split.reads.details$junID %>% length()
  
    ## Get mis-spliced and junctions that were not mis-spliced in the first round
    df_all_misspliced <- readRDS(file = paste0(folder.name, "/", cluster, "_raw_distances_tidy.rds")) %>%
      data.table::as.data.table() %>%
      dplyr::distinct(ref_junID, .keep_all = T)
    
    all_not_misspliced <- all.split.reads.details %>%
      data.table::as.data.table() %>%
      dplyr::filter( !(junID %in% df_all_misspliced$ref_junID) )
    
    if (intersect(all_not_misspliced$junID, df_all_misspliced$ref_junID) %>% length() > 0){
      logger::log_info("Error: some of the never misspliced introns appear within the 'distances' file!")
      break;
    }
    if (any(df_all_misspliced$ref_junID %>% duplicated()) || any(all_not_misspliced$junID %>% duplicated())) {
      logger::log_info("Error: none of the mis-spliced introns should within the never mis-spliced introns!")
      break;
    }
    
  } else {
    logger::log_info(paste0("File: '", folder.name, "/", cluster, "_distances_raw.rds' doesn't exist."))
    break;
  }
  
  
  
  #######################################
  ## FROM THE JUNCTIONS NOT PAIRED IN THE FIRST ROUND,
  ## WE TRY TO PAIR THEM AGAIN.
  ## THE JUNCTIONS THAT CANNOT BE PAIRED IN THIS SECOND ROUND, ARE 
  ## CONSIDERED AS NEVER MIS-SPLICED
  #######################################
  
  num_sample <- 1
  junc_ignore <- NULL
  junc_not_misspliced <- NULL
  
  
  if ( replace || 
       (!file.exists(paste0(folder.name, "/", cluster, "_all_notmisspliced.rds")) && 
        !file.exists(paste0(folder.name, "/", cluster, "_all_misspliced_not_paired.rds"))) ) {
    
    ## Per each sample from the current cluster, we obtain all junctions, counts and ratios
    for (sample in samples) { 
      
      # sample <- samples[1]
      # sample <- samples[16]
      
      # Every sampleID is unique, so the result of this comparison should be equal to 1
      
      if ( length(which((colnames(split.read.counts) == sample) == TRUE)) == 1 ) { 
      
        
        ## Obtain the split reads from the current sample
        split_read_counts_sample <- split.read.counts %>%
          as_tibble()%>%
          dplyr::select(junID, all_of(sample %>% as.character())) 
        
        
        ## Discard junctions with zero reads within the current sample
        split_read_counts_sample[split_read_counts_sample == 0] <- NA
        split_read_counts_sample <- split_read_counts_sample %>% drop_na()
        
        
        ## Add count data to yet-never-misspliced junctions
        split_reads_details_sample <- all_not_misspliced %>% 
          inner_join(y = split_read_counts_sample, by = "junID") %>%
          dplyr::rename( counts = all_of(sample %>% as.character()) ) %>% 
          as_tibble()
        
        
        ## Do some QC
        index <- runif(n = 1, 1, split_reads_details_sample %>% nrow()) %>% as.integer()
        data <- split_reads_details_sample[index, c("junID", "counts")]
        data_col <- split_read_counts_sample %>%
          dplyr::filter(junID == data$junID)
        if (data_col[[2]] != data$counts) {
          logger::log_info("Error: QC failed!")
          break;
        }
        
        ## Print an informative message
        logger::log_info(cluster, "...")
        
        ###################################################
        ##########    ANNOTATED JUNC   ####################
        ###################################################
        
        all_annotated <- split_reads_details_sample %>%
          dplyr::filter(type == "annotated") %>% 
          GenomicRanges::GRanges()
        
        logger::log_info("sample '", num_sample, "': ", length(all_annotated$junID), " initial annotated junctions.")
        
        all_annotated_forward <- all_annotated[all_annotated@strand == "+"]
        all_annotated_reverse <- all_annotated[all_annotated@strand == "-"]
        
        
        ###################################################
        ##########    NOVEL DONOR      ####################
        ###################################################
        
        all_donor <- split_reads_details_sample %>%
          dplyr::filter(type == "novel_donor") %>% 
          GenomicRanges::GRanges()
        
        logger::log_info("sample '", num_sample, "': ", length(all_donor$junID), " initial novel donor junctions.")
        
        all_donor_forward <- all_donor[all_donor@strand == "+"]
        all_donor_reverse <- all_donor[all_donor@strand == "-"]
        
        
        ###################################################
        ##########    NOVEL ACCEPTOR      #################
        ###################################################
        
        all_acceptor <- split_reads_details_sample %>%
          dplyr::filter(type == "novel_acceptor") %>% 
          GenomicRanges::GRanges()
        
        logger::log_info("sample '", num_sample, "': ", length(all_acceptor$junID), " initial novel acceptor junctions.")
        
        all_acceptor_forward <- all_acceptor[all_acceptor@strand == "+"]
        all_acceptor_reverse <- all_acceptor[all_acceptor@strand == "-"]
        
        
        
        ###################################################
        ##########    NOVEL DONOR PAIRINGS      ###########
        ###################################################
        
        if ( all_donor_forward %>% length() > 0 && all_donor_reverse %>% length() > 0 &&
             all_annotated_forward %>% length() > 0 && all_annotated_reverse %>% length() > 0 ) {
          

          ##  NOVEL-DONOR & FORWARD STRAND ---------------------------------
          
          overl_df <- GenomicRanges::findOverlaps(query = all_donor_forward,
                                                  subject = all_annotated_forward,
                                                  ignore.strand = FALSE,
                                                  type = "end")
          ignore_df <- all_annotated_forward[subjectHits(overl_df),]$junID
          ref_annotated_d_forward <- all_annotated_forward[-subjectHits(overl_df),]$junID
          
          

          ## NOVEL-DONOR & REVERSE STRAND -------------------------------------
          
          overl_dr <- GenomicRanges::findOverlaps(query = all_donor_reverse,
                                                  subject = all_annotated_reverse,
                                                  type = "start",
                                                  ignore.strand = FALSE)
          ignore_dr <- all_annotated_reverse[subjectHits(overl_dr),]$junID
          ref_annotated_d_reverse <- all_annotated_reverse[-subjectHits(overl_dr),]$junID
          
          
        } else {
          ignore_df <- NULL
          ref_annotated_d_forward <- NULL
          ignore_dr <- NULL
          ref_annotated_d_reverse <- NULL
        }
        

        
        ###################################################
        ##########    NOVEL ACCEPTOR PAIRINGS      ########
        ###################################################

        
        if (all_acceptor_forward %>% length() > 0 &&
            all_acceptor_reverse %>% length() > 0 &&
            all_annotated_forward %>% length() > 0 &&
            all_annotated_reverse %>% length() > 0) {
          
          
          ##  NOVEL-ACCEPTOR & FORWARD STRAND  -------------------------------------
          
          overl_af <- GenomicRanges::findOverlaps(query = all_acceptor_forward,
                                                  subject = all_annotated_forward,
                                                  type = "start",
                                                  ignore.strand = FALSE)
          ignore_af <- all_annotated_forward[subjectHits(overl_af),]$junID
          ref_annotated_a_forward <- all_annotated_forward[-subjectHits(overl_af),]$junID
          
          

          ##  NOVEL-ACCEPTOR & REVERSE STRAND  -------------------------------------

          
          overl_ar <- GenomicRanges::findOverlaps(query = all_acceptor_reverse,
                                                  subject = all_annotated_reverse,
                                                  type = "end",
                                                  ignore.strand = FALSE)
          ignore_ar <- all_annotated_reverse[subjectHits(overl_ar),]$junID
          ref_annotated_a_reverse <- all_annotated_reverse[-subjectHits(overl_ar),]$junID
          
        } else {
          ignore_af <- NULL
          ref_annotated_a_forward <- NULL
          ignore_ar <- NULL
          ref_annotated_a_reverse <- NULL
        }
        
        
        ## These are the annotated introns that have been paired, and therefore, are potentially mis-spliced
        
        junc_ignore <- c(junc_ignore,
                         ignore_af, 
                         ignore_ar,
                         ignore_df, 
                         ignore_dr) %>% unique()
        
        
        ## These are the annotated introns that have not paired in these second round, and therefore, are NOT mis-spliced
        
        junc_not_misspliced <- c(junc_not_misspliced,
                                 ref_annotated_d_forward, 
                                 ref_annotated_d_reverse,
                                 ref_annotated_a_forward, 
                                 ref_annotated_a_reverse) %>% unique()
        
        
        
        
        logger::log_info(paste0(Sys.time(), " --> Sample '", num_sample, "' processed // ", 
                     c(ref_annotated_d_forward, ref_annotated_d_reverse,
                       ref_annotated_a_forward, ref_annotated_a_reverse) %>% unique() %>% length(), 
                     " junctions not mis-spliced // ", c(ignore_af,
                                                         ignore_ar,
                                                         ignore_df,
                                                         ignore_dr) %>% unique %>% length(), " mis-spliced junctions."))
        
        
        num_sample <- num_sample + 1
        
        rm(ignore_af)
        rm(ignore_ar)
        rm(ignore_df)
        rm(ignore_dr)
        rm(ref_annotated_d_forward)
        rm(ref_annotated_d_reverse)
        rm(ref_annotated_a_forward) 
        rm(ref_annotated_a_reverse)
        rm(overl_df)
        rm(overl_dr)
        rm(overl_af)
        rm(overl_ar)
        rm(split_read_counts_sample)
        
        rm(split_reads_details_sample)
        rm(all_annotated)
        rm(all_annotated_forward)
        rm(all_annotated_reverse)
        rm(all_donor)
        rm(all_donor_forward)
        rm(all_donor_reverse)
        rm(all_acceptor)
        rm(all_acceptor_forward)
        rm(all_acceptor_reverse)
        rm(data)
        rm(data_col)
        gc()
        
      }
    }
    
    ###########################################
    ## ENSURE NEVER MIS-SPLICED JUNCTIONS
    ## ARE NEVER MIS-SPLICED ACROSS ALL SAMPLES
    ###########################################
    
    
    ## It is possible that some introns not mis-spliced in a particular sample, have been mis-spliced in other samples
    ## We discard them as they have the potential to be mis-spliced and, therefore, they cannot be considered
    ## as never mis-spliced introns
    
    if ( (which(junc_not_misspliced %in% junc_ignore) %>% length()) > 0) {
      junc_not_misspliced <- junc_not_misspliced[-which(junc_not_misspliced %in% junc_ignore)]  
    }
    
    
    ###########################################
    ## SAVE THE RESULTS AND RELEASE MEMORY
    ###########################################
    
    
    if ( intersect(junc_not_misspliced, junc_ignore) %>% length() > 0) {
      
      logger::log_info(paste0(Sys.time(), " - Error! There are overlapping not-misspliced and misspliced junctions!"))
      break;
      
    } else { 
      
      logger::log_info(paste0(Sys.time(), " - ", junc_not_misspliced %>% length(), " not mis-spliced junctions!"))
      folder_path <- paste0(folder.name, "/not-misspliced/")
      dir.create(file.path(folder_path), showWarnings = F)
      
      saveRDS(object = junc_not_misspliced, file = paste0(folder_path, "/", cluster, "_all_notmisspliced.rds"))
      saveRDS(object = junc_ignore, file = paste0(folder_path, "/", cluster, "_all_misspliced_not_paired.rds"))
      
      logger::log_info(paste0(Sys.time(), " - results saved!"))
      
    }
    
    
    
    
    ## RELEASE SOME MEMORY
    rm(num_sample)
    rm(junc_ignore)
    rm(junc_not_misspliced)
    #rm(df_junc_counts)
    rm(split_reads_details_sample)
    rm(all_annotated)
    rm(all_annotated_forward)
    rm(all_annotated_reverse)
    rm(all_donor)
    rm(all_donor_forward)
    rm(all_donor_reverse)
    rm(all_acceptor)
    rm(all_acceptor_forward)
    rm(all_acceptor_reverse)
    rm(overl_df)
    rm(ignore_df)
    rm(ref_annotated_d_forward)
    rm(overl_dr)
    rm(ignore_dr)
    rm(ref_annotated_d_reverse)
    rm(overl_af)
    rm(ignore_af)
    rm(ref_annotated_a_forward)
    rm(overl_ar)
    rm(ignore_ar)
    rm(ref_annotated_a_reverse)
    gc()
    
  }
  
}

