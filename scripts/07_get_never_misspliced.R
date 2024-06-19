#' Title
#'
#' @param cluster 
#' @param samples 
#' @param split_read_counts 
#' @param all_split_reads_details 
#' @param folder_name 
#'
#' @return
#' @export
#'
#' @examples
get_never_misspliced <- function(cluster, 
                                 samples,
                                 split_read_counts,
                                 all_split_reads_details,
                                 folder_name,
                                 replace) {
  
  #######################################
  ## GET ALL JUNCTIONS THAT WERE NOT 
  ## MIS-SPLICED IN THE 1ST JXN PAIRING
  #######################################
  
  logger::log_info(paste0(Sys.time(), " - filtering junctions that are potentially not mis-spliced..."))
  
  if ( file.exists(paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds")) ) {
    
    # ## This should be zero
    if ( setdiff(all_split_reads_details$junID, split_read_counts$junID) %>% length() > 0) {
      logger::log_info("ERROR")
      break;
    }
    
    all_split_reads_details$junID %>% length()
  
    ## Get mis-spliced and junctions that were not mis-spliced in the first round
    df_all_misspliced <- readRDS(file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds")) %>%
      data.table::as.data.table() %>%
      dplyr::distinct(ref_junID, .keep_all = T)
    
    all_not_misspliced <- all_split_reads_details %>%
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
    logger::log_info(paste0("File: '", folder_name, "/", cluster, "_distances_raw.rds' doesn't exist."))
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
       (!file.exists(paste0(folder_name, "/", cluster, "_all_notmisspliced.rds")) && 
        !file.exists(paste0(folder_name, "/", cluster, "_all_misspliced_not_paired.rds"))) ) {
    
    ## Per each sample from the current cluster, we obtain all junctions, counts and ratios
    for (sample in samples) { 
      
      # sample <- samples[1]
      # sample <- samples[16]
      
      # Every sampleID is unique, so the result of this comparison should be equal to 1
      
      if ( length(which((colnames(split_read_counts) == sample) == TRUE)) == 1 ) { 
      
        
        ## Obtain the split reads from the current sample
        split_read_counts_sample <- split_read_counts %>%
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
      folder_name <- paste0(folder_name, "/not-misspliced/")
      dir.create(file.path(folder_name), showWarnings = F)
      
      saveRDS(object = junc_not_misspliced, file = paste0(folder_name, "/", cluster, "_all_notmisspliced.rds"))
      saveRDS(object = junc_ignore, file = paste0(folder_name, "/", cluster, "_all_misspliced_not_paired.rds"))
      
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

