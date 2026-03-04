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
                                        phastcons.bw.path,
                                        cdts.bw.path,
                                        db.introns = NULL,
                                        cluster = NULL,
                                        folder.name = NULL,
                                        intron.size = 100,
                                        phastcons.type = 17) {
  
  if (is.null(db.introns) && !is.null(cluster)) {
    db.introns <- readRDS(file = paste0(folder.name, "/", cluster, "_db.introns.rds")) %>%
      distinct(ref_junID, .keep_all = T)
  }

 
  # Open the connection to the BigWigFile once 
  if (!file.exists(phastcons.bw.path)) {
    stop("Conservation file '", phastcons.bw.path, "' does not exist!") %>% logger::log_info()
  } else {
    BigWigFilePhastCons <- BigWigFile(phastcons.bw.path)
  }
  if (!file.exists(cdts.bw.path)) {
    stop("Constraint file '", cdts.bw.path, "' does not exist!") %>% logger::log_info()
  } else {
    BigWigFileCDTS <- BigWigFile(cdts.bw.path)
  }
  
  
  
  ###########################################
  ## PHASTCONS SCORES
  ###########################################
  
  for (i_size in intron.size) {
    
    # i_size <-intron.size[1]
    
    for (p_type in phastcons.type) {
      
      colname_5ss <- paste0("mean_phastCons", p_type, "way5ss_", i_size)
      colname_3ss <- paste0("mean_phastCons", p_type, "way3ss_", i_size)

      ## Calculate donor (5'ss) scores -----------------------------------------      
      logger::log_info(i_size, "bp - Calculating PhastCons" , p_type, " scores overlapping donor sequences...")
      gr <- GenomicRanges::GRanges(seqnames = db.introns %>% seqnames(),
                                   ranges = IRanges(start = db.introns %>% start(), end = db.introns %>% start() + i_size),
                                   strand = db.introns %>% strand())

      mcols(gr)[["junID"]] <- (db.introns %>% as.character())
      phastCons_5ss <- GetScoresFromBWRegion(bw_path = BigWigFilePhastCons, gr, summaryFun = "mean", col_name = colname_5ss) 
     
      
      ## Calculate acceptor (3'ss) scores --------------------------------------      
      logger::log_info(i_size, "bp - Calculating PhastCons" , p_type, " scores overlapping acceptor sequences...")
      gr <- GenomicRanges::GRanges(seqnames = db.introns %>% seqnames(),
                                   ranges = IRanges(start = db.introns %>% end() - i_size, db.introns %>% end()),
                                   strand = db.introns %>% strand)

      mcols(gr)[["junID"]] <- (db.introns %>% as.character())  
      phastCons_3ss <- GetScoresFromBWRegion(bw_path = BigWigFilePhastCons, gr = gr, summaryFun = "mean", col_name = colname_3ss) 

      
      ## Add phastcons columns to gr object   
      mcols(db.introns)[[colname_5ss]] <- mcols(phastCons_5ss)[[colname_5ss]]
      mcols(db.introns)[[colname_3ss]] <- mcols(phastCons_3ss)[[colname_3ss]]

   
    }
    
    ###########################################
    ## CDTS SCORES
    ###########################################
    
    
    # CDTS scores were calculated using a window size of 550bp, sliding every 10bp, 
    # attributting the calculated CDTS score across the 550-bp window to the middle 10-bp bin.
    
    # Hence, to calculate the CDTS scores overlapping with the proximal intronic regions,
    # we obtain the 10bp regions overlapping with the proximal intronic region and get the mean
    # CDTS score value
    
    ## Add columns to master data
    colname_5ss <- paste0("mean_CDTS5ss_", i_size)
    colname_3ss <- paste0("mean_CDTS3ss_", i_size)


    ## Calculate donor scores
    logger::log_info(i_size, "bp - CDTS calculating donor sequences...")
    gr <- GenomicRanges::GRanges(seqnames = db.introns %>% seqnames(),
                                 ranges = IRanges(start = db.introns %>% start(), end = db.introns %>% start() + i_size),
                                 strand = db.introns %>% strand)
    mcols(gr)[["junID"]] <- (db.introns %>% as.character())
    CDTS_5ss <- GetScoresFromBWRegion(bw_path = BigWigFileCDTS, gr = gr, summaryFun = "mean", col_name = colname_5ss)


    ## Calculate acceptor scores
    logger::log_info(i_size, "bp - CDTS calculating acceptor sequences....")
    gr <- GenomicRanges::GRanges(seqnames = db.introns %>% seqnames(),
                                 ranges = IRanges(start = db.introns %>% end() - i_size, end = db.introns %>% end()),
                                 strand = db.introns %>% strand)
    mcols(gr)[["junID"]] <- (db.introns %>% as.character())
    CDTS_3ss <- GetScoresFromBWRegion(bw_path = BigWigFileCDTS, gr = gr, summaryFun = "mean", col_name = colname_3ss) 
   
    
    ## Add CDTS columns to gr object
    mcols(db.introns)[[colname_5ss]] <- mcols(CDTS_5ss)[[colname_5ss]]
    mcols(db.introns)[[colname_3ss]] <- mcols(CDTS_3ss)[[colname_3ss]]
  }
  
  
  ####################
  ## Save results
  ####################
  
  if (!is.null(cluster)) {

    saveRDS(object = db.introns, file = paste0(folder.name, "/", cluster, "_db.introns.rds"))
    logger::log_info("CDTS and PhastCons scores added! Database updated!")

    rm(gr)
    rm(db.introns)
    rm(phastCons_5ss)
    rm(phastCons_3ss)

  } else {
    return(db.introns)
  }
}




# Functions -------------------------------------------------------------------------------------------

GetScoresFromBWRegion <- function(bw_path, gr, summaryFun  = "mean", col_name){
    
  gr_w_scores <- summary(object = bw_path, gr, size = 1L, type = summaryFun) %>% unlist()
  
  stopifnot((width(gr) == width(gr_w_scores)))

  mcols(gr)[[col_name]] <- gr_w_scores$score
  
  return(gr)
  
}




# Functions -------------------------------------------------------------------------------------------

# to check whether the two methods of getting conservation return comparable answers
# they do and also comparable to the UCSC browser, although slightly different for some which i believe if accounted for by rounding/aggregating errors
# gr <- GRanges(c("chr8:69105760-69105804", "chr8:69105760-69105804"))
#
# system.time(
#
#   bw_method <- GetConservationScoreForRegionsBW(bw_path = "/data/conservation/phastCons/hg38.phastCons7way.bw",
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

