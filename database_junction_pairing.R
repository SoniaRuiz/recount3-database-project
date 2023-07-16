####################################################
## NOVEL JUNCTION PAIRING WITH ANNOTATED INTRONS ###
####################################################


## JUNCTION PAIRING AND DISTANCES--------------------

get_distances <- function(cluster, 
                          samples,
                          split_read_counts,
                          all_split_reads_details,
                          folder_name) {
  
  ## DECLARE SOME VARIABLES
  num_sample <- 1
  
  ## Per each sample from the current cluster, we obtain all junctions, counts and ratios
  for (sample in samples) { 
    
    # sample <- samples[1]
    # sample <- samples[2]
    # sample <- "54285"
    
    # if ( !file.exists(paste0(folder_name, "/", cluster, "_", sample, "_distances.rds")) ) {
    
    # Every sampleID is unique, so the result should be equal to 1
    if ( length(which((colnames(split_read_counts) == sample) == TRUE)) == 1 ) { 
      
      print(cluster)
      print(paste0(Sys.time(), " - sample '", num_sample, "'"))
      
      ############################################################
      ## OBTAINING ALL JUNCTIONS AND COUNTS FOR THE CURRENT SAMPLE
      ############################################################
      
      split_read_counts_sample <- split_read_counts %>%
        as_tibble() %>%
        dplyr::select(junID, all_of(sample %>% as.character())) %>%
        drop_na() 
      
      split_read_counts_sample <- split_read_counts_sample[(split_read_counts_sample[,2] >0),]
      
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
      data <- split_reads_details_sample[index, c(1, split_reads_details_sample %>% ncol())]
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
        GRanges()
      
      all_donor_forward <- all_donor[all_donor@strand == "+"]
      all_donor_reverse <- all_donor[all_donor@strand == "-"]
      
      
      ###################################################
      ##########    NOVEL ACCEPTOR      #################
      ###################################################
      all_acceptor <- split_reads_details_sample %>%
        dplyr::filter(type == "novel_acceptor") %>% 
        GRanges()
      
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
        ## the one with the highest number of counts. In case many ref junctions present the same number of counts,
        ## the closest one will be chosen.
        ## We do not consider the novel junction ambiguous at this stage, as the pairing has not been done yet.
        df_nd_f <- df_nd_f %>% 
          dplyr::group_by(novel_junID) %>% 
          dplyr::filter(ref_counts == max(ref_counts)) %>%
          dplyr::filter(abs(distance) == min(abs(distance))) %>%
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
        ## the one with the highest number of counts. In case many ref junctions present the same number of counts,
        ## the closest one will be chosen.
        ## We do not consider the novel junction ambiguous at this stage, as the pairing has not been done yet.
        df_nd_r <- df_nd_r %>% 
          dplyr::group_by(novel_junID) %>%
          dplyr::filter(ref_counts == max(ref_counts)) %>%
          dplyr::filter(abs(distance) == min(abs(distance))) %>%
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
        ## the one with the highest number of counts. In case many ref junctions present the same number of counts,
        ## the closest one will be chosen.
        ## We do not consider the novel junction ambiguous at this stage, as the pairing has not been done yet.
        df_na_f <- df_na_f %>% 
          dplyr::group_by(novel_junID) %>%
          dplyr::filter(ref_counts == max(ref_counts)) %>%
          dplyr::filter(abs(distance) == min(abs(distance))) %>%
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
        
        
        df_na_r <- df_na_r %>% 
          dplyr::group_by(novel_junID) %>%
          dplyr::filter(ref_counts == max(ref_counts)) %>%
          dplyr::filter(abs(distance) == min(abs(distance))) %>%
          ungroup() 
        
      } else {
        df_na_r <- NULL
      }
      
      #################################
      ########  SAVE RESULTS  #########
      #################################
      
      distances <- rbind(df_nd_f, df_nd_r,
                         df_na_f, df_na_r)
      
      distances <- distances %>% 
        dplyr::select(-novel_width, -ref_width, -sample)
      
      # if (any(distances$novel_junID %>% duplicated())) {
      #   print("ERROR! At least one of the novel junctions has been paired twice!")
      #   break;
      # }
      
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
    
    # } else {
    #   print(paste0("Sample '", sample, "' exists!"))
    # }
    
  }
}


extract_distances <- function(cluster,
                              samples,
                              folder_name) {
  
  ## Obtain the distances across all samples
  df_all <- map_df(samples, function(sample) { 
    
    # sample <- samples[1]
    # sample <- "GTEX-ZVT2-0426-SM-5E44S.1"
    
    file_name <- paste0(folder_name, "/", cluster, "_", sample, "_distances.rds")
    
    if (file.exists(file_name)) {
      
      print(paste0("sample: ", sample))
      
      df <- readRDS(file = file_name)
      
      return(df)
    } else {
      print(paste0("sample: ", sample, " doesn't exist!"))
      break;
    }
    
  })
  saveRDS(object = df_all %>%
            distinct(novel_junID, ref_junID, .keep_all = T) %>%
            mutate(tissue = cluster),
          file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds"))
  
  ## RELEASE SOME MEMORY
  rm(df_all)
  gc()
  
}





## NEVER MIS-SPLICED --------------------------------

get_never_misspliced <- function(cluster, 
                                 samples,
                                 split_read_counts,
                                 all_split_reads_details,
                                 folder_name) {
  
  #######################################
  ## GET ALL JUNCTIONS THAT WERE NOT 
  ## MIS-SPLICED IN THE 1ST JXN PAIRING
  #######################################
  
  print(paste0(Sys.time(), " - filtering junctions that are potentially not mis-spliced..."))
  
  if ( file.exists(paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds")) ) {
    
    # ## This should be zero
    if ( setdiff(all_split_reads_details$junID, 
                 split_read_counts$junID) %>% length() > 0) {
      print("ERROR")
      break;
    }
    
    all_split_reads_details$junID %>% length()
    
    split_read_counts_tidy <- split_read_counts %>%
      inner_join(y = all_split_reads_details%>% dplyr::select(seqnames, start, end, width, strand, junID),
                 by = c("junID" = "junID")) 
    
    ############################
    ## QC 
    ############################
    
    ## Split_Read_Counts file
    ind <- which(str_detect(string = split_read_counts_tidy$junID, pattern = "\\*"))
    if (ind %>% length() > 0) {
      split_read_counts_tidy[ind, "junID"] <- str_replace(string = split_read_counts_tidy[ind, "junID"]$junID, 
                                                          pattern = "\\*", 
                                                          replacement = split_read_counts_tidy[ind, "strand"]$strand %>% as.character() )
    }
    split_read_counts_tidy <- split_read_counts_tidy %>%
      dplyr::select(-seqnames, -start, -end, -width, -strand)
    any(str_detect(string = split_read_counts_tidy$junID, pattern = "\\*"))
    
    ## Annotated_split_reads file
    ind <- which(str_detect(string = all_split_reads_details$junID, pattern = "\\*"))
    if (ind %>% length() > 0) {
      all_split_reads_details[ind, "junID"] <- str_replace(string = all_split_reads_details[ind, "junID"]$junID, 
                                                           pattern = "\\*", 
                                                           replacement = all_split_reads_details[ind, "strand"]$strand %>% as.character() )
    }
    
    ## QC - 2: this should be zero
    ind <- which(str_detect(string = split_read_counts_tidy$junID, pattern = "\\*"))
    if (ind %>% length() > 0) {
      print("ERROR")
      break;
    }
    
    ## Get mis-spliced and junctions that were not mis-spliced in the first round
    df_all_misspliced <- readRDS(file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds")) %>%
      data.table::as.data.table() %>%
      dplyr::distinct(ref_junID, .keep_all = T)
    
    all_not_misspliced <- all_split_reads_details %>%
      data.table::as.data.table() %>%
      dplyr::filter( !(junID %in% df_all_misspliced$ref_junID) )
    
    if (intersect(all_not_misspliced$junID, df_all_misspliced$ref_junID) %>% length() > 0){
      print("Error: some of the never misspliced introns appear within the 'distances' file!")
      break;
    }
    if (any(df_all_misspliced$ref_junID %>% duplicated()) || any(all_not_misspliced$junID %>% duplicated())) {
      print("Error: none of the mis-spliced introns should within the never mis-spliced introns!")
      break;
    }
  } else {
    print(paste0("File: '", folder_name, "/", cluster, "_distances_raw.rds' doesn't exist."))
    break;
  }
  
  
  
  #######################################
  ## FROM THE JUNCTIONS NOT PAIRED IN THE FIRST ROUND,
  ## WE TRY TO PAIR THEM AGAIN.
  ## THE ONES THAT ARE NOT PAIRED AGAIN, ARE 
  ## CONSIDERED AS NEVER MIS-SPLICED
  #######################################
  
  
  num_sample <- 1
  junc_ignore <- NULL
  junc_not_misspliced <- NULL
  
  
  # if ( !file.exists(paste0(folder_name, "/not-misspliced/", cluster, "_all_notmisspliced.rds")) || 
  #      !file.exists(paste0(folder_name, "/not-misspliced/", cluster, "_all_misspliced_not_paired.rds")) ) {
  
  ## Per each sample from the current cluster, we obtain all junctions, counts and ratios
  for (sample in samples) { 
    
    # sample <- samples[1]
    # sample <- samples[16]
    
    if ( length(which((colnames(split_read_counts_tidy) == sample) == TRUE)) == 1 ) { # Every sampleID is unique, so the result should be equal to 1
      
      
      split_read_counts_sample <- split_read_counts_tidy %>%
        as_tibble()%>%
        dplyr::select(junID, all_of(sample %>% as.character())) 
      
      split_read_counts_sample[split_read_counts_sample == 0] <- NA
      split_read_counts_sample <- split_read_counts_sample %>%
        drop_na()
      
      split_reads_details_sample <- all_not_misspliced %>% 
        inner_join(y = split_read_counts_sample,
                   by = "junID") %>%
        dplyr::rename( counts = all_of(sample %>% as.character()) )
      
      split_reads_details_sample <- split_reads_details_sample %>% as_tibble()
      all_not_misspliced <- all_not_misspliced %>% as_tibble()
      
      
      ## Do some QC
      index <- runif(n = 1, 1, split_reads_details_sample %>% nrow()) %>% as.integer()
      data <- split_reads_details_sample[index, c(1, split_reads_details_sample %>% ncol())]
      data_col <- split_read_counts_sample %>%
        dplyr::filter(junID == data$junID)
      if (data_col[[2]] != data$counts) {
        print("Error: QC failed!")
        break;
      }
      
      ## Print an informative message
      print(cluster)
      
      ###################################################
      ##########    ANNOTATED JUNC   ####################
      ###################################################
      all_annotated <- split_reads_details_sample %>%
        dplyr::filter(type == "annotated") %>% 
        GenomicRanges::GRanges()
      
      print(paste0("sample '", num_sample, "': ", 
                   length(all_annotated$junID), " initial ref junctions."))
      
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
      all_acceptor <- split_reads_details_sample%>%
        dplyr::filter(type == "novel_acceptor") %>% 
        GRanges()
      
      all_acceptor_forward <- all_acceptor[all_acceptor@strand == "+"]
      all_acceptor_reverse <- all_acceptor[all_acceptor@strand == "-"]
      
      if (all_donor_forward %>% length() > 0 &&
          all_donor_reverse %>% length() > 0 &&
          all_annotated_forward %>% length() > 0 &&
          all_annotated_reverse %>% length() > 0) {
        
        ####################################################
        ########  NOVEL-DONOR & FORWARD STRAND  ############
        ####################################################
        
        overl_df <- GenomicRanges::findOverlaps(query = all_donor_forward,
                                                subject = all_annotated_forward,
                                                ignore.strand = FALSE,
                                                type = "end")
        ignore_df <- all_annotated_forward[subjectHits(overl_df),]$junID
        ref_annotated_d_forward <- all_annotated_forward[-subjectHits(overl_df),]$junID
        
        
        ####################################################
        ########  NOVEL-DONOR & REVERSE STRAND  ############
        ####################################################
        
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
      
      ####################################################
      ########  NOVEL-ACCEPTOR & FORWARD STRAND  #########
      ####################################################
      
      if (all_acceptor_forward %>% length() > 0 &&
          all_acceptor_reverse %>% length() > 0 &&
          all_annotated_forward %>% length() > 0 &&
          all_annotated_reverse %>% length() > 0) {
        
        overl_af <- GenomicRanges::findOverlaps(query = all_acceptor_forward,
                                                subject = all_annotated_forward,
                                                type = "start",
                                                ignore.strand = FALSE)
        ignore_af <- all_annotated_forward[subjectHits(overl_af),]$junID
        ref_annotated_a_forward <- all_annotated_forward[-subjectHits(overl_af),]$junID
        
        
        ####################################################
        ########  NOVEL-ACCEPTOR & REVERSE STRAND  #########
        ####################################################
        
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
                       ignore_af, ignore_ar,
                       ignore_df, ignore_dr) %>% unique()
      
      
      junc_not_misspliced <- c(junc_not_misspliced,
                               ref_annotated_d_forward, 
                               ref_annotated_d_reverse,
                               ref_annotated_a_forward, 
                               ref_annotated_a_reverse) %>% unique()
      
      
      
      
      print(paste0(Sys.time(), " --> Sample '", num_sample, "' processed // ", 
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
  
  ## It is possible that some introns not mis-spliced in a particular samples, have been mis-spliced in other sample
  ## We discard them as they have the potential to be mis-spliced and, therefore, they cannot be considered
  ## as never mis-spliced introns
  if ( (which(junc_not_misspliced %in% junc_ignore) %>% length()) > 0) {
    junc_not_misspliced <- junc_not_misspliced[-which(junc_not_misspliced %in% junc_ignore)]  
  }
  
  
  ###########################################
  ## SAVE THE RESULTS AND RELEASE MEMORY
  ###########################################
  
  if ( intersect(junc_not_misspliced, junc_ignore) %>% length() > 0) {
    print(paste0(Sys.time(), " - Error! There are overlapping not-misspliced and misspliced junctions!"))
    break;
    
  } else { 
    
    print(paste0(Sys.time(), " - ", junc_not_misspliced %>% length(), " not mis-spliced junctions!"))
    folder_name <- paste0(folder_name, "/not-misspliced/")
    dir.create(file.path(folder_name), showWarnings = F)
    
    saveRDS(object = junc_not_misspliced, 
            file = paste0(folder_name, "/", cluster, "_all_notmisspliced.rds"))
    
    saveRDS(object = junc_ignore, 
            file = paste0(folder_name, "/", cluster, "_all_misspliced_not_paired.rds"))
    
    print(paste0(Sys.time(), " - results saved!"))
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
  #}
  
}

