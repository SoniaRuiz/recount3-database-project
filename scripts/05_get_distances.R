

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
get_distances <- function(project_id,
                          cluster, 
                          samples,
                          split_read_counts,
                          all_split_reads_details,
                          folder_name,
                          replace) {
  
  ## DECLARE SOME VARIABLES
  num_sample <- 1
  
  ## Per each sample from the current cluster, we obtain all junctions, counts and ratios
  for (sample in samples) { 
    
    # sample <- samples[1]
    # sample <- samples[2]
    # sample <- "54285"
    
    if ( !file.exists(paste0(folder_name, "/", cluster, "_", sample, "_distances.rds")) || replace ) {
      
      # Every sampleID is unique, so the result should be equal to 1
      if ( length(which((colnames(split_read_counts) == sample) == TRUE)) == 1 ) { 
        
       
        print(paste0(Sys.time(), " - ", project_id, " - ", cluster, " sample '", num_sample, "'"))
        
        ############################################################
        ## OBTAINING ALL JUNCTIONS AND COUNTS FOR THE CURRENT SAMPLE
        ############################################################
        
        split_read_counts_sample <- split_read_counts %>%
          as_tibble() %>%
          dplyr::select(junID, all_of(sample %>% as.character())) %>%
          drop_na() 
        
        split_read_counts_sample <- split_read_counts_sample[(split_read_counts_sample[,2] > 0),]
        
        split_read_counts %>% nrow()
        split_read_counts_sample %>% nrow()
        split_read_counts_sample %>% head()
        
        
        split_reads_details_sample <- all_split_reads_details %>%
          inner_join(y = split_read_counts_sample,
                     by = "junID")  %>%
          dplyr::rename(counts = all_of(sample %>% as.character())) %>%
          as_tibble()
        
        split_reads_details_sample %>% head()
        
        ## Do some QC
        index <- runif(n = 1, 1, split_reads_details_sample %>% nrow()) %>% as.integer()
        data <- split_reads_details_sample[index, c("junID", "counts")]
        data_col <- split_read_counts_sample %>%
          dplyr::filter(junID == data$junID)
        if (data_col[, colnames(data_col) == sample] != data$counts) {
          print("Error: QC failed!")
          break;
        }
        
        ###################################################
        ##########    ANNOTATED JUNC   ####################
        ###################################################
        
        all_annotated <- split_reads_details_sample %>%
          dplyr::filter(type == "annotated") %>% 
          GenomicRanges::GRanges()
        
        all_annotated_forward <- all_annotated[all_annotated@strand == "+"]
        all_annotated_reverse <- all_annotated[all_annotated@strand == "-"]
        
        
        ###################################################
        ##########    NOVEL DONOR      ####################
        ###################################################
        
        all_donor <- split_reads_details_sample %>%
          dplyr::filter(type == "novel_donor") %>% 
          GenomicRanges::GRanges()
        
        all_donor_forward <- all_donor[all_donor@strand == "+"]
        all_donor_reverse <- all_donor[all_donor@strand == "-"]
        
        
        ###################################################
        ##########    NOVEL ACCEPTOR      #################
        ###################################################
        
        all_acceptor <- split_reads_details_sample %>%
          dplyr::filter(type == "novel_acceptor") %>% 
          GenomicRanges::GRanges()
        
        all_acceptor_forward <- all_acceptor[all_acceptor@strand == "+"]
        all_acceptor_reverse <- all_acceptor[all_acceptor@strand == "-"]
        
        
        ####################################################
        ########  NOVEL-DONOR & FORWARD STRAND  ############
        ####################################################
        
        overlaps <- GenomicRanges::findOverlaps(query = all_donor_forward,
                                                subject = all_annotated_forward,
                                                ignore.strand = FALSE,
                                                type = "end")
        
        novel_junctions <- all_donor_forward[queryHits(overlaps),]
        ref_junctions <- all_annotated_forward[subjectHits(overlaps),]
        
        distance <- start(ref_junctions) - start(novel_junctions)
        
        if ( distance %>% length() > 0 ) {
          
          df_nd_f <- data.frame(sample = sample,
                                type = "novel_donor", 
                                distance = distance,
                                novel_junID = novel_junctions$junID,
                                novel_counts = novel_junctions$counts,
                                novel_seq = novel_junctions %>% seqnames() %>% as.character(),
                                novel_start = novel_junctions %>% start() %>% as.character(),
                                novel_end = novel_junctions %>% end() %>% as.character(),
                                novel_strand = novel_junctions %>% strand() %>% as.character(),
                                novel_width = novel_junctions %>% width() %>% as.character(),
                                ref_junID = ref_junctions$junID,
                                ref_counts = ref_junctions$counts,
                                ref_seq = ref_junctions %>% seqnames(),
                                ref_start = ref_junctions %>% start(),
                                ref_strand = ref_junctions %>% strand(),
                                ref_end = ref_junctions %>% end(),
                                ref_width = ref_junctions %>% width(),
                                stringsAsFactors = FALSE) %>% as_tibble()
          
          ## In case the novel junction has been attached to many ref junctions, the ref junction chosen will be 
          ## the closest one. In case many ref junctions present the same number of counts, 
          ## the one with the highest number of supporting split reads will be chosen.
          ## We do not consider novel junctions as ambiguous at this stage, given that the pairing has not been done yet.
          df_nd_f <- df_nd_f %>% 
            dplyr::group_by(novel_junID) %>% 
            dplyr::filter(abs(distance) == min(abs(distance))) %>%
            dplyr::filter(ref_counts == max(ref_counts)) %>%
            ungroup() 
          
        } else {
          df_nd_f <- NULL
        }
        
        ####################################################
        ########  NOVEL-DONOR & REVERSE STRAND  ############
        ####################################################
        
        overlaps <- GenomicRanges::findOverlaps(query = all_donor_reverse,
                                                subject = all_annotated_reverse,
                                                ignore.strand = FALSE,
                                                type = "start")
        
        novel_junctions <- all_donor_reverse[queryHits(overlaps),]
        ref_junctions <- all_annotated_reverse[subjectHits(overlaps),]
        
        distance <- end(novel_junctions) - end(ref_junctions) 
        
        if (distance %>% length() > 0) {
          df_nd_r <- data.frame(sample = sample,
                                type = "novel_donor", 
                                distance = distance,
                                novel_junID = novel_junctions$junID,
                                novel_counts = novel_junctions$counts,
                                novel_seq = novel_junctions %>% seqnames(),
                                novel_start = novel_junctions %>% start(),
                                novel_end = novel_junctions %>% end(),
                                novel_strand = novel_junctions %>% strand(),
                                novel_width = novel_junctions %>% width(),
                                ref_junID = ref_junctions$junID,
                                ref_counts = ref_junctions$counts,
                                ref_seq = ref_junctions %>% seqnames(),
                                ref_start = ref_junctions %>% start(),
                                ref_end = ref_junctions %>% end(),
                                ref_strand = ref_junctions %>% strand(),
                                ref_width = ref_junctions %>% width(),
                                stringsAsFactors = FALSE)
          
          
          ## In case the novel junction has been attached to many ref junctions, the ref junction chosen will be 
          ## the closest one. In case many ref junctions present the same number of counts, 
          ## the one with the highest number of supporting split reads will be chosen.
          ## We do not consider novel junctions as ambiguous at this stage, given that the pairing has not been done yet.
          df_nd_r <- df_nd_r %>% 
            dplyr::group_by(novel_junID) %>%
            dplyr::filter(abs(distance) == min(abs(distance))) %>%
            dplyr::filter(ref_counts == max(ref_counts)) %>%
            ungroup() 
          
        } else {
          df_nd_r <- NULL
        }
        
        
        ####################################################
        ########  NOVEL-ACCEPTOR & FORWARD STRAND  #########
        ####################################################
        
        
        overlaps <- GenomicRanges::findOverlaps(query = all_acceptor_forward,
                                                subject = all_annotated_forward,
                                                ignore.strand = FALSE,
                                                type = "start")
        
        novel_junctions <- all_acceptor_forward[queryHits(overlaps),]
        ref_junctions <- all_annotated_forward[subjectHits(overlaps),]
        
        distance <- end(novel_junctions) - end(ref_junctions)
        
        if (distance %>% length() > 0) {
          
          df_na_f <- data.frame(sample = sample,
                                type = "novel_acceptor", 
                                distance = distance,
                                novel_junID = novel_junctions$junID,
                                novel_counts = novel_junctions$counts,
                                novel_seq = novel_junctions %>% seqnames(),
                                novel_start = novel_junctions %>% start(),
                                novel_end = novel_junctions %>% end(),
                                novel_strand = novel_junctions %>% strand(),
                                novel_width = novel_junctions %>% width(),
                                ref_junID = ref_junctions$junID,
                                ref_counts = ref_junctions$counts,
                                ref_seq = ref_junctions %>% seqnames(),
                                ref_start = ref_junctions %>% start(),
                                ref_end = ref_junctions %>% end(),
                                ref_strand = ref_junctions %>% strand(),
                                ref_width = ref_junctions %>% width(),
                                stringsAsFactors = FALSE)
          
          
          ## In case the novel junction has been attached to many ref junctions, the ref junction chosen will be 
          ## the closest one. In case many ref junctions present the same number of counts, 
          ## the one with the highest number of supporting split reads will be chosen.
          ## We do not consider novel junctions as ambiguous at this stage, given that the pairing has not been done yet.
          df_na_f <- df_na_f %>% 
            dplyr::group_by(novel_junID) %>%
            dplyr::filter(abs(distance) == min(abs(distance))) %>%
            dplyr::filter(ref_counts == max(ref_counts)) %>%
            ungroup() 
          
        } else {
          df_na_f <- NULL
        }
        
        
        ####################################################
        ########  NOVEL-ACCEPTOR & REVERSE STRAND  #########
        ####################################################
        
        overlaps <- GenomicRanges::findOverlaps(query = all_acceptor_reverse,
                                                subject = all_annotated_reverse,
                                                ignore.strand = FALSE,
                                                type = "end")
        
        
        novel_junctions <- all_acceptor_reverse[queryHits(overlaps),]
        ref_junctions <- all_annotated_reverse[subjectHits(overlaps),]
        
        distance <- start(ref_junctions) - start(novel_junctions)
        
        if (distance %>% length() > 0) {
          df_na_r <- data.frame(sample = sample,
                                type = "novel_acceptor", 
                                distance = distance,
                                
                                novel_junID = novel_junctions$junID,
                                novel_counts = novel_junctions$counts,
                                novel_seq = novel_junctions %>% seqnames(),
                                novel_start = novel_junctions %>% start(),
                                novel_end = novel_junctions %>% end(),
                                novel_strand = novel_junctions %>% strand(),
                                novel_width = novel_junctions %>% width(),
                                
                                ref_junID = ref_junctions$junID,
                                ref_counts = ref_junctions$counts,
                                ref_seq = ref_junctions %>% seqnames(),
                                ref_start = ref_junctions %>% start(),
                                ref_end = ref_junctions %>% end(),
                                ref_strand = ref_junctions %>% strand(),
                                ref_width = ref_junctions %>% width(),
                                
                                stringsAsFactors = FALSE)
          
          ## In case the novel junction has been attached to many ref junctions, the ref junction chosen will be 
          ## the closest one. In case many ref junctions present the same number of counts, 
          ## the one with the highest number of supporting split reads will be chosen.
          ## We do not consider novel junctions as ambiguous at this stage, given that the pairing has not been done yet.
          df_na_r <- df_na_r %>% 
            dplyr::group_by(novel_junID) %>%
            dplyr::filter(abs(distance) == min(abs(distance))) %>%
            dplyr::filter(ref_counts == max(ref_counts)) %>%
            ungroup() 
          
        } else {
          df_na_r <- NULL
        }
        
        #################################
        ########  SAVE RESULTS  #########
        #################################
        
        distances <- rbind(df_nd_f, df_nd_r,
                           df_na_f, df_na_r) %>% dplyr::select(-novel_width, -ref_width, -sample)
        
        saveRDS(object = distances,
                file = paste0(folder_name, "/", cluster, "_", sample, "_distances.rds"))
        
        num_sample <- num_sample + 1
        
        ## FREE-UP SOME MEMORY - once every 50 samples processed
        if (num_sample %% 50 == 0) {
          rm(split_reads_details_sample)
          rm(split_read_counts_sample)
          rm(all_annotated)
          rm(all_annotated_forward)
          rm(all_annotated_reverse)
          rm(all_donor)
          rm(all_donor_forward)
          rm(all_donor_reverse)
          rm(all_acceptor)
          rm(all_acceptor_forward)
          rm(all_acceptor_reverse)
          rm(overlaps)
          rm(distances)
          rm(df_nd_f)
          rm(df_nd_r)
          rm(df_na_f)
          rm(df_na_r)
          rm(novel_junctions)
          rm(ref_junctions)
          gc()
        }
        
      } else {
        print(paste0("Error: sample '", sample, "' not found within the split-reads file!"))
        break;
      }
      
   } else {
     message(cluster, " sample '", sample, "' exists!")
   }
    
  }
  
}
