source("/home/dzhang/projects/constraint_conservation_wd/constraint_conservation/conservation/conservation_general/conservation_general_functions_from_bw.R")
source("/home/dzhang/projects/constraint_conservation_wd/constraint_conservation/constraint/constraint_general/constraint_general_functions_bw.R")

#' Title
#' Calculates the CDTS and PhastCons20 mean scores of the sequences overlapping the /+35bp sequence at the donor site of the intron 
#' (i.e. / = the exon-intron junction), and the /+35bp sequence at the acceptor splice 
#' site of the intron (i.e. the intron-exon junction)
#' @param cluster 
#' @param db_introns Dataframe of introns to calculate the scores
#' @param folder_name 
#'
#' @return
#' @export
#'
#' @examples
generate_cdts_phastcons_scores <- function(db_introns = NULL,
                                           cluster = NULL,
                                           folder_name = NULL,
                                           intron_size = c(100),
                                           phastcons_type = c(17)) {
  
  
  # load(file = "/data/constraint/CDTS_percentile_N7794_unrelated_all_chrs_gr.rda") 
  
  
  if (is.null(db_introns) && !is.null(cluster)) {
    ## Load the IDB 
    db_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds")) %>%
      distinct(ref_junID, .keep_all = T)
  }
  
  if ( !str_detect(string = db_introns[1, ]$seqnames, pattern = "chr") ) {
    db_introns <- db_introns %>%
      GRanges() %>% 
      diffloop::addchr()
  } else {
    db_introns <- db_introns %>%
      GRanges() 
  }
 
  
  ###########################################
  ## PHASTCONS SCORES
  ###########################################
  
  for (i_size in intron_size) {
    
    # i_size <-intron_size[1]
    
    for (p_type in phastcons_type) {
      
      # p_type <- phastcons_type[1]
      
      logger::log_info("PhastCons", p_type, "...")
      
      bw_path = paste0("/data/conservation/phastCons/hg38.phastCons", p_type, "way.bw")
      
      if ( !file.exists(bw_path) ) {
        paste0("PhastCons file '", bw_path, "' does not exist!") %>% logger::log_info()
        break;
      }
        
      # i_size <- 35
      
      ## Calculate donor scores
      logger::log_info(i_size, "bp - Calculating PhastCons" , p_type, " scores overlapping donor sequences...")
      gr <- GenomicRanges::GRanges(seqnames = db_introns %>% seqnames(),
                                   ranges = IRanges(start = db_introns %>% start(), 
                                                    end = db_introns %>% start() + i_size))
      values(gr) <- DataFrame(junID = (db_introns) %>% as.character() )
      phastCons_5ss <- get_conservation_score_for_regions_bw(bw_path = bw_path,
                                                             gr = gr, 
                                                             summaryFun = "mean") %>% 
        as_tibble() %>%
        dplyr::rename_with(.fn = ~paste0(., "5ss_", i_size), .cols = paste0("mean_phastCons", p_type, "way"))
      
      
      
      
      ## Calculate acceptor scores
      logger::log_info(i_size, "bp - Calculating PhastCons" , p_type, " scores overlapping acceptor sequences...")
      gr <- GenomicRanges::GRanges(seqnames = db_introns %>% seqnames(),
                                   ranges = IRanges(start = db_introns %>% end() - i_size, 
                                                    end = db_introns %>% end()))
      values(gr) <- DataFrame(junID = (db_introns) %>% as.character() )
      phastCons_3ss <- get_conservation_score_for_regions_bw(bw_path = bw_path,
                                                             gr = gr, 
                                                             summaryFun = "mean") %>% 
        as_tibble()  %>%
        dplyr::rename_with(.fn = ~paste0(., "3ss_", i_size), .cols = paste0("mean_phastCons", p_type, "way"))
      
      
      
      ## Add columns to master data
      db_introns <- db_introns %>%
        as_tibble() %>%
        left_join(phastCons_5ss %>% dplyr::select(junID, paste0("mean_phastCons", p_type, "way5ss_", i_size)),
                   by = "junID") %>%
        left_join(phastCons_3ss %>% dplyr::select(junID, paste0("mean_phastCons", p_type, "way3ss_", i_size)),
                   by = "junID") %>%
        GRanges()
      
      
    }
    
    ###########################################
    ## CDTS SCORES
    ###########################################
    
    
    # CDTS scores were calculated using a window size of 550bp, sliding every 10bp, 
    # attributting the calculated CDTS score across the 550-bp window to the middle 10-bp bin.
    
    # Hence, to calculate the CDTS scores overlapping with the proximal intronic regions,
    # we obtain the 10bp regions overlapping with the proximal intronic region and get the mean
    # CDTS score value
    

    bw_path = paste0(dependencies_folder, "/CDTS_percentile_N7794_unrelated_all_chrs.bw")
    if (!file.exists(bw_path)) {
      paste0("Constraint file '", bw_path, "' does not exist!") %>% logger::log_info()
      break;
    }


    ## Calculate donor scores
    logger::log_info(i_size, "bp - CDTS calculating donor sequences...")
    gr <- GenomicRanges::GRanges(seqnames = db_introns %>% seqnames(),
                                 ranges = IRanges(start = db_introns %>% start(),
                                                  end = db_introns %>% start() + i_size))
    values(gr) <- DataFrame(junID = (db_introns) %>% as.character() )
    CDTS_5ss <- get_constraint_score_for_regions_bw(bw_path = bw_path,
                                                    gr = gr,
                                                    summaryFun = "mean") %>%
      as_tibble() %>%
      dplyr::rename_with(.fn = ~paste0(., "5ss_", i_size), .cols = mean_CDTS)




    ## Calculate acceptor scores
    logger::log_info(i_size, "bp - CDTS calculating acceptor sequences....")
    gr <- GenomicRanges::GRanges(seqnames = db_introns %>% seqnames(),
                                 ranges = IRanges(start = db_introns %>% end() - i_size,
                                                  end = db_introns %>% end()))
    values(gr) <- DataFrame(junID = (db_introns) %>% as.character() )
    CDTS_3ss <- get_constraint_score_for_regions_bw(bw_path = bw_path,
                                                    gr = gr,
                                                    summaryFun = "mean") %>%
      as_tibble()  %>%
      dplyr::rename_with(.fn = ~paste0(., "3ss_", i_size), .cols = paste0("mean_CDTS"))

    
    
    ## Add columns to master data
    db_introns <- db_introns %>%
      as_tibble() %>%
      left_join(CDTS_5ss %>% dplyr::select(junID, paste0("mean_CDTS5ss_", i_size)),
                 by = "junID") %>%
      left_join(CDTS_3ss %>% dplyr::select(junID, paste0("mean_CDTS3ss_", i_size)),
                 by = "junID") %>%
      GRanges()
    
    # 
    # ## 5'ss
    # 
    # logger::log_info(i_size, "bp - Calculating CDTS scores overlapping donor sequences...")
    # overlaps_5ss <- GenomicRanges::findOverlaps(query = GRanges(seqnames = db_introns %>% seqnames,
    #                                                             ranges = IRanges(start = db_introns %>% start,
    #                                                                              end = db_introns %>% start + i_size),
    #                                                             strand = db_introns %>% strand),
    #                                             subject = CDTS_percentile_N7794_unrelated_all_chrs_gr,
    #                                             type = "any",
    #                                             select = "all")
    # overlaps_5ss
    # 
    # 
    # 
    # hits_overlaps_5ss <- data.frame(ref_junID = db_introns[queryHits(overlaps_5ss),]$ref_junID,
    #                                 CDTS_5ss_100 = CDTS_percentile_N7794_unrelated_all_chrs_gr[subjectHits(overlaps_5ss),]$CDTS) %>%
    #   as_tibble()
    # 
    # ## 3'ss
    # 
    # logger::log_info(i_size, "bp - Calculating CDTS scores overlapping acceptor sequences...")
    # overlaps_3ss <- GenomicRanges::findOverlaps(query = GRanges(seqnames = db_introns %>% seqnames,
    #                                                             ranges = IRanges(start = db_introns %>% end - i_size,
    #                                                                              end = db_introns %>% end),
    #                                                             strand = db_introns$strand),
    #                                             subject = CDTS_percentile_N7794_unrelated_all_chrs_gr,
    #                                             type = "any",
    #                                             select = "all")
    # overlaps_3ss
    # 
    # 
    # 
    # hits_overlaps_3ss <- data.frame(ref_junID = db_introns[queryHits(overlaps_3ss),]$ref_junID,
    #                                 CDTS_3ss_100 = CDTS_percentile_N7794_unrelated_all_chrs_gr[subjectHits(overlaps_3ss),]$CDTS) %>%
    #   as_tibble()
    # 
    # 
    # ## Merge 5' and 3'
    # logger::log_info(i_size, "bp - merging results ...")
    # 
    # db_introns <- db_introns %>%
    #   as_tibble() %>%
    #   left_join(y = hits_overlaps_5ss,
    #             by = "ref_junID") %>%
    #   left_join(y = hits_overlaps_3ss,
    #             by = "ref_junID") %>%
    #   as_tibble() %>%
    #   dplyr::group_by(ref_junID) %>%
    #   mutate(CDTS_5ss_100_mean = CDTS_5ss_100 %>% mean) %>%
    #   mutate(CDTS_3ss_100_mean = CDTS_3ss_100 %>% mean) %>%
    #   ungroup() %>%
    #   distinct(ref_junID, .keep_all = T) %>%
    #   dplyr::select(-c(CDTS_5ss_100, CDTS_3ss_100))
    # 
    
    # 
    # ## Add CDTS scores 
    # 
    # db_introns <- db_introns %>%
    #   as_tibble() %>%
    #   left_join(CDTS_scores,
    #             by = "ref_junID") 
    
    
    
    # # i_size <- 35
    # 
    # ## Calculate donor scores
    # logger::log_info(i_size, "bp - CDTS calculating donor sequences")
    # gr <- GenomicRanges::GRanges(seqnames = db_introns %>% seqnames(),
    #                              ranges = IRanges(start = db_introns %>% start(), 
    #                                               end = db_introns %>% start() + i_size))
    # values(gr) <- DataFrame(ref_junID = (db_introns)$ref_junID )
    # CDTS_5ss <- get_constraint_score_for_regions_bw(bw_path = bw_path,
    #                                                 gr = gr,
    #                                                 summaryFun = "mean") %>% 
    #   as_tibble() %>%
    #   dplyr::rename_with(.fn = ~paste0(., "5ss_", i_size), .cols = mean_CDTS)
    # 
    # 
    # 
    # 
    # ## Calculate acceptor scores
    # logger::log_info(i_size, "bp - CDTS calculating acceptor sequences")
    # gr <- GenomicRanges::GRanges(seqnames = db_introns %>% seqnames(),
    #                              ranges = IRanges(start = db_introns %>% end() - i_size, 
    #                                               end = db_introns %>% end()))
    # values(gr) <- DataFrame(ref_junID = (db_introns)$ref_junID )
    # CDTS_3ss <- get_constraint_score_for_regions_bw(bw_path = bw_path,
    #                                                 gr = gr, 
    #                                                 summaryFun = "mean") %>% 
    #   as_tibble()  %>%
    #   dplyr::rename_with(.fn = ~paste0(., "3ss_", i_size), .cols = paste0("mean_CDTS"))
    # 
    
    
    ## Add columns to master data
    # db_introns <- db_introns %>%
    #   # as_tibble() %>%
    #   # inner_join(CDTS_5ss %>% dplyr::select(ref_junID, paste0("mean_CDTS5ss_", i_size)),
    #   #            by = "ref_junID") %>%
    #   # inner_join(CDTS_3ss %>% dplyr::select(ref_junID, paste0("mean_CDTS3ss_", i_size)),
    #   #            by = "ref_junID") %>%
    #   GRanges()
    
  }
  
  
  ####################
  ## Save results
  ####################
  
  if (!is.null(cluster)) {

    saveRDS(object = db_introns,
            file = paste0(folder_name, "/", cluster, "_db_introns.rds"))

    logger::log_info(" - CDTS and PhastCons scores added! Database updated!")


    rm(gr)
    rm(db_introns)
    rm(phastCons_5ss)
    rm(phastCons_3ss)

  } else {
    return(db_introns)
  }
      
  
}
# 
# ## CONVERT TO BW FILE
# 
# CDTS_percentile_N7794_unrelated <- read.csv(file = "/data/constraint/coord_CDTS_percentile_N7794unrelated.txt", header = T, sep = "\t")
# 
# CDTS_percentile_N7794_unrelated %>% head()
# CDTS_percentile_N7794_unrelated %>% nrow()
# 
# CDTS_percentile_N7794_unrelated_all_chrs_gr
# 
# genome_build <- "hg38"
# 
# chrominfo <- fetchExtendedChromInfoFromUCSC(genome_build)
# 
# chromosome_lengths_df <-
#   chrominfo %>%
#   filter(UCSC_seqlevel %in% seqlevels(CDTS_percentile_N7794_unrelated_all_chrs_gr))
# 
# chromosome_lengths <-
#   chromosome_lengths_df[["UCSC_seqlength"]]
# 
# names(chromosome_lengths) <- chromosome_lengths_df[["UCSC_seqlevel"]]
# 
# CDTS_percentile_N7794_unrelated_all_chrs_gr <- CDTS_percentile_N7794_unrelated_all_chrs_gr %>% sortSeqlevels()
# 
# seqlengths(CDTS_percentile_N7794_unrelated_all_chrs_gr) <- chromosome_lengths
# genome(CDTS_percentile_N7794_unrelated_all_chrs_gr) <- genome_build
# 
# CDTS_percentile_N7794_unrelated_all_chrs_gr$percentile <- NULL
# CDTS_percentile_N7794_unrelated_all_chrs_gr$score <- CDTS_percentile_N7794_unrelated_all_chrs_gr$CDTS
# CDTS_percentile_N7794_unrelated_all_chrs_gr$CDTS <- NULL
# 
# end(CDTS_percentile_N7794_unrelated_all_chrs_gr) <- end(CDTS_percentile_N7794_unrelated_all_chrs_gr) - 1
# 
# export.bw(object = CDTS_percentile_N7794_unrelated_all_chrs_gr, con = "/data/constraint/CDTS_percentile_N7794_unrelated_all_chrs.bw")

