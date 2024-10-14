
#' Title
#' Per junction, this function calculates the percentage of transcripts in which this junction may appear, that are protein-coding
#' @param gtf.version Version of the reference transcriptome to use. In this case it has been used '105' corresponding
#' to Ensembl v105
#' @param database.folder Path to the local folder that stores the database to be produced and the files needed to produce it 
#'
#' @return
#' @export
#'
#' @examples
GenerateTranscriptBiotypePercentage <- function(gtf.version,
                                                dependencies.folder,
                                                database.folder) {
  
  
  #######################################
  ## GET THE TRANSCRIPT BIOTYPE
  #######################################
  
  logger::log_info(paste0(Sys.time(), " - loading the human reference transcriptome ... "))
  
  ## Import HUMAN REFERENCE transcriptome
  homo_sapiens_hg <- rtracklayer::import(con = paste0(dependencies.folder, "/Homo_sapiens.GRCh38.",gtf.version,".chr.gtf")) %>% as_tibble()
  
  ## Get transcripts
  transcripts_hg <- homo_sapiens_hg %>%
    filter(type == "transcript") %>% 
    dplyr::select(transcript_id, transcript_biotype, gene_id)
  
  transcripts_hg %>% head()
  
  logger::log_info("starting protein-coding percentage calculation ...")
  
  #######################################
  ## LOAD ALL SPLIT READS
  #######################################

  logger::log_info("loading the split reads ... ")

  ## LOAD the all split reads from all recount3 GTEx projects
  all_split_reads_details_all_tissues <- readRDS(file = paste0(database.folder, "/all_split_reads_qc_level1.rds") )

  all_split_reads_details_all_tissues %>% head()
  all_split_reads_details_all_tissues %>% nrow()
  
  ##########################################
  ## Link SPLIT READS TO THEIR TRANSCRIPTS
  ##########################################
  
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = all_split_reads_details_all_tissues$junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    all_split_reads_details_all_tissues[ind, "junID"] <-
      str_replace(string = all_split_reads_details_all_tissues[ind, "junID"]$junID,
                  pattern = "\\*",
                  replacement = all_split_reads_details_all_tissues[ind, "strand"]$strand %>% as.character() )
  }
  
  ## Merge datasets to add transcript biotype
  logger::log_info(paste0(Sys.time(), " - adding human reference transcript info to the paired junctions..."))

  df_all_junctions <- all_split_reads_details_all_tissues %>% unnest(tx_id_junction) %>% data.table::as.data.table()
  transcripts_hg <- transcripts_hg %>%  data.table::as.data.table()

  df_all_junctions_tx_hg <- df_all_junctions %>%
    dplyr::select(-gene_id) %>%
    inner_join(y = transcripts_hg,
              by = c("tx_id_junction" = "transcript_id"))

  if ( !identical(df_all_junctions_tx_hg$gene_id.x %>% unlist, df_all_junctions_tx_hg$gene_id.y) ) {
    logger::log_info("ERROR! The gene IDs from the reference and the split reads in the database do not match!")
    stop("ERROR! The gene IDs from the reference and the split reads in the database do not match!");
  }

  
  ################################################
  ## CALCULATE THE TRANSCRIPT BIOTYPE PERCENTAGE
  ################################################
  
  logger::log_info(df_all_junctions$junID %>% unique() %>% length(), " total number of junctions.")
  logger::log_info("getting transcript biotype percentage per junction...")
  
  transcripts_hg_percentage <- df_all_junctions_tx_hg %>% 
    dplyr::group_by(junID, transcript_biotype) %>%
    distinct(tx_id_junction, .keep_all = T) %>% 
    summarise(n = n()) %>% 
    mutate(percent = (n / sum(n)) * 100) %>%
    ungroup()
  
  transcripts_hg_percentage_tidy <- transcripts_hg_percentage %>%
    dplyr::select(-n) %>%
    spread(transcript_biotype,percent) %>%
    replace(is.na(.), 0)
  
  saveRDS(object = transcripts_hg_percentage_tidy, 
          file = file.path(database.folder, "all_split_reads_qc_level1_PC_biotype.rds"))
  
  logger::log_info(paste0(Sys.time(), " - results saved!"))
  
  
}