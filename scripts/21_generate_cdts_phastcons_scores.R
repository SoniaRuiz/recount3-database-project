library(rtracklayer)

#' Title
#' Calculates the CDTS and PhastCons20 mean scores of the sequences overlapping the /+35bp sequence at the donor site of the intron 
#' (i.e. / = the exon-intron junction), and the /+35bp sequence at the acceptor splice 
#' site of the intron (i.e. the intron-exon junction)
#' @param cluster 
#' @param db.introns Dataframe of introns to calculate the scores
#' @param folder.name 
#'
#' @return
#' @export
#'
#' @examples
GenerateCdtsPhastconsScores <- function(dependencies.folder,
                                        db.introns = NULL,
                                        cluster = NULL,
                                        folder.name = NULL,
                                        intron.size = 100,
                                        phastcons.type = 17) {
  
  if (is.null(db.introns) && !is.null(cluster)) {
    db.introns <- readRDS(file = paste0(folder.name, "/", cluster, "_db.introns.rds")) %>%
      distinct(ref_junID, .keep_all = T)
  }
  if ( !str_detect(string = db.introns[1, ]$seqnames, pattern = "chr") ) {
    db.introns <- db.introns %>% GRanges() %>% diffloop::addchr()
  } else {
    db.introns <- db.introns %>% GRanges() 
  }
 
  
  ###########################################
  ## PHASTCONS SCORES
  ###########################################
  
  for (i_size in intron.size) {
    
    # i_size <-intron.size[1]
    
    for (p_type in phastcons.type) {
      
      # p_type <- phastcons.type[1]
      
      logger::log_info("PhastCons", p_type, "...")
      
      bw_path = paste0(dependencies.folder, "/hg38.phastCons", p_type, "way.bw")
      
      if (!file.exists(bw_path)) {
        stop("PhastCons file '", bw_path, "' does not exist!") %>% logger::log_info()
      }
        
      # i_size <- 35
      
      ## Calculate donor scores
      logger::log_info(i_size, "bp - Calculating PhastCons" , p_type, " scores overlapping donor sequences...")
      gr <- GenomicRanges::GRanges(seqnames = db.introns %>% seqnames(),
                                   ranges = IRanges(start = db.introns %>% start(), 
                                                    end = db.introns %>% start() + i_size))
      values(gr) <- DataFrame(junID = (db.introns) %>% as.character() )
      phastCons_5ss <- get_conservation_score_for_regions_bw(bw_path = bw_path,
                                                             gr = gr, 
                                                             summaryFun = "mean") %>% 
        as_tibble() %>%
        dplyr::rename_with(.fn = ~paste0(., "5ss_", i_size), .cols = paste0("mean_phastCons", p_type, "way"))
      
      
      
      
      ## Calculate acceptor scores
      logger::log_info(i_size, "bp - Calculating PhastCons" , p_type, " scores overlapping acceptor sequences...")
      gr <- GenomicRanges::GRanges(seqnames = db.introns %>% seqnames(),
                                   ranges = IRanges(start = db.introns %>% end() - i_size, 
                                                    end = db.introns %>% end()))
      values(gr) <- DataFrame(junID = (db.introns) %>% as.character() )
      phastCons_3ss <- get_conservation_score_for_regions_bw(bw_path = bw_path,
                                                             gr = gr, 
                                                             summaryFun = "mean") %>% 
        as_tibble()  %>%
        dplyr::rename_with(.fn = ~paste0(., "3ss_", i_size), .cols = paste0("mean_phastCons", p_type, "way"))
      
      
      
      ## Add columns to master data
      db.introns <- db.introns %>%
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
    

    bw_path = file.path(dependencies.folder, "CDTS_percentile_N7794_unrelated_all_chrs.bw")
    if (!file.exists(bw_path)) {
      stop("Constraint file '", bw_path, "' does not exist!") %>% logger::log_info()
    }


    ## Calculate donor scores
    logger::log_info(i_size, "bp - CDTS calculating donor sequences...")
    gr <- GenomicRanges::GRanges(seqnames = db.introns %>% seqnames(),
                                 ranges = IRanges(start = db.introns %>% start(),
                                                  end = db.introns %>% start() + i_size))
    values(gr) <- DataFrame(junID = (db.introns) %>% as.character() )
    CDTS_5ss <- get_constraint_score_for_regions_bw(bw_path = bw_path,
                                                    gr = gr,
                                                    summaryFun = "mean") %>%
      as_tibble() %>%
      dplyr::rename_with(.fn = ~paste0(., "5ss_", i_size), .cols = mean_CDTS)




    ## Calculate acceptor scores
    logger::log_info(i_size, "bp - CDTS calculating acceptor sequences....")
    gr <- GenomicRanges::GRanges(seqnames = db.introns %>% seqnames(),
                                 ranges = IRanges(start = db.introns %>% end() - i_size,
                                                  end = db.introns %>% end()))
    values(gr) <- DataFrame(junID = (db.introns) %>% as.character() )
    CDTS_3ss <- get_constraint_score_for_regions_bw(bw_path = bw_path,
                                                    gr = gr,
                                                    summaryFun = "mean") %>%
      as_tibble()  %>%
      dplyr::rename_with(.fn = ~paste0(., "3ss_", i_size), .cols = paste0("mean_CDTS"))

    
    
    ## Add columns to master data
    db.introns <- db.introns %>%
      as_tibble() %>%
      left_join(CDTS_5ss %>% dplyr::select(junID, paste0("mean_CDTS5ss_", i_size)),
                 by = "junID") %>%
      left_join(CDTS_3ss %>% dplyr::select(junID, paste0("mean_CDTS3ss_", i_size)),
                 by = "junID") %>%
      GRanges()
    
  }
  
  
  ####################
  ## Save results
  ####################
  
  if (!is.null(cluster)) {

    saveRDS(object = db.introns, file = paste0(folder.name, "/", cluster, "_db.introns.rds"))
    logger::log_info(" - CDTS and PhastCons scores added! Database updated!")

    rm(gr)
    rm(db.introns)
    rm(phastCons_5ss)
    rm(phastCons_3ss)

  } else {
    return(db.introns)
  }
}




# Functions -------------------------------------------------------------------------------------------

get_conservation_score_for_regions_bw <- function(bw_path, gr, summaryFun  = "mean"){
  
  BigWigFile <- BigWigFile(bw_path)
  
  phast_cons_score <- bw_path %>% str_replace(".*/", "") %>% str_extract("phastCons.*way")
  
  gr_w_scores <- summary(BigWigFile, gr, size = 1L, type = summaryFun) %>% unlist()
  
  stopifnot((width(gr)  == width(gr_w_scores)))
  
  elementMetadata(gr)[[str_c(summaryFun, "_", phast_cons_score)]] <- gr_w_scores$score
  
  gr_w_scores <- gr
  
  return(gr_w_scores)
  
}

get_constraint_score_for_regions_bw <- function(bw_path, gr, summaryFun  = "mean"){
  
  BigWigFile <- BigWigFile(bw_path)
  
  gr_w_scores <- summary(BigWigFile, gr, size = 1L, type = summaryFun) %>% unlist()
  
  stopifnot((width(gr)  == width(gr_w_scores)))
  
  elementMetadata(gr)[[str_c(summaryFun, "_CDTS")]] <- gr_w_scores$score
  
  gr_w_scores <- gr
  
  return(gr_w_scores)
  
}


# Functions -------------------------------------------------------------------------------------------

# to check whether the two methods of getting conservation return comparable answers
# they do and also comparable to the UCSC browser, although slightly different for some which i believe if accounted for by rounding/aggregating errors
# gr <- GRanges(c("chr8:69105760-69105804", "chr8:69105760-69105804"))
#
# system.time(
#
#   bw_method <- get_conservation_score_for_regions_bw(bw_path = "/data/conservation/phastCons/hg38.phastCons7way.bw",
#                                                      gr = gr, summaryFun = "mean")
# )
#
# system.time(
#   old_method <- get_conservation_score_for_regions(conserv_score_to_load = "phast_cons_7", gr = gr, summaryFun = "mean")
# )
#
# bw_method$mean_phastCons7_old <- old_method$mean_phast_cons_7
#
# bw_method

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

