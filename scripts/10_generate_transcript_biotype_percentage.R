# gtf.version <- 110
# database.folder <- "/home/sruiz/PROJECTS/splicing-accuracy-manuscript/database/splicing_2read/110/"
# results.folder <- "/home/sruiz/PROJECTS/splicing-accuracy-manuscript/results/splicing_2read/110/"

#' Title
#' Per junction, this function calculates the percentage of transcripts in which this junction may appear, that are protein-coding
#' @param gtf.version Version of the reference transcriptome to use. In this case it has been used '105' corresponding
#' to Ensembl v105
#' @param database.folder Path to the local folder that stores the database to be produced and the files needed to produce it 
#' @param results.folder Local path to the folder that contains the results of the analyses performed
#'
#' @return
#' @export
#'
#' @examples
generate_transcript_biotype_percentage <- function(gtf.version,
                                                   database.folder,
                                                   results.folder) {
  
  
  #######################################
  ## GET THE TRANSCRIPT BIOTYPE
  #######################################
  
  logger::log_info(paste0(Sys.time(), " - loading the human reference transcriptome ... "))
  
  ## Import HUMAN REFERENCE transcriptome
  homo_sapiens_hg <- rtracklayer::import(con = paste0(dependencies_folder, 
                                                      "/Homo_sapiens.GRCh38.",gtf.version,".chr.gtf")) %>% as_tibble()
  
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
  #logger::log_info(object.size(all_split_reads_details_all_tissues), units = "Gb")

  
  ## Merge datasets to add transcript biotype
  logger::log_info(paste0(Sys.time(), " - adding human reference transcript info to the paired junctions..."))

  df_all_junctions <- all_split_reads_details_all_tissues %>% 
    unnest(tx_id_junction) %>% 
    data.table::as.data.table()
 
  transcripts_hg <- transcripts_hg %>% 
    data.table::as.data.table()

  df_all_junctions %>% head()
  df_all_junctions %>% nrow()
  transcripts_hg %>% head()
  
  df_all_junctions_tx_hg <- df_all_junctions %>%
    dplyr::select(-gene_id) %>%
    inner_join(y = transcripts_hg,
              by = c("tx_id_junction" = "transcript_id"))

  df_all_junctions_tx_hg %>% nrow
  df_all_junctions_tx_hg %>% head
  
  # logger::log_info(object.size(df_all_junctions_tx_hg), units = "Gb")

  if ( !identical(df_all_junctions_tx_hg$gene_id.x %>% unlist, 
                  df_all_junctions_tx_hg$gene_id.y) ) {
    logger::log_info("ERROR! The gene IDs from the reference and the split reads in the database do not match!")
    break;
  }

  
  ################################################
  ## CALCULATE THE TRANSCRIPT BIOTYPE PERCENTAGE
  ################################################
  
  logger::log_info(df_all_junctions$junID %>% unique() %>% length(), " total number of junctions.")
  
  
  ## Only keep non ambiguous jxn (i.e. junctions belonging to multiple genes)
  
  # junID_OK <- df_all_junctions %>%
  #   dplyr::group_by(junID) %>%
  #   distinct(gene_id, .keep_all = T) %>%
  #   dplyr::count() %>%
  #   filter(n == 1) %>%
  #   pull(junID)
  # 
  # 
  #   if (!identical(junID_OK %>% unique() %>% length(),
  #                  df_all_junctions$junID %>% unique() %>% length())) {
  #     logger::log_info("ERROR! There are still some split reads ambiguously linked to multiple genes.")
  #     break;
  #   }
  
  # df_all_junctions_tx_hg%>%
  #   filter(junID == "chr1:1290695-1290894:+")
  # junID_OK <- c("chr1:1290695-1290894:+","chr1:2583496-2584124:+")
  
  logger::log_info("getting transcript biotype percentage per junction...")
  
  transcripts_hg_percentage <- df_all_junctions_tx_hg %>% 
    #filter(junID %in% junID_OK) %>%
    dplyr::group_by(junID, transcript_biotype) %>%
    distinct(tx_id_junction, .keep_all = T) %>% 
    summarise(n = n()) %>% 
    mutate(percent = (n / sum(n)) * 100) %>%
    ungroup()
  
  # transcripts_hg_percentage %>% filter(junID == "chr1:1290695-1290894:+")
  # transcripts_hg_percentage %>% filter(junID == "chr1:2583496-2584124:+")
  
  # logger::log_info(Sys.time(), " - getting protein-coding percentage per junction...")
  # 
  # transcripts_hg_percentage_PC <- transcripts_hg_percentage %>% 
  #   dplyr::group_by(junID) %>%
  #   rowwise() %>%
  #   mutate(percent = ifelse (transcript_biotype == "protein_coding", percent, 0)) %>%
  #   filter(percent == max(percent)) %>%
  #   dplyr::group_by(junID) %>%
  #   filter(transcript_biotype == "protein_coding") %>%
  #   ungroup() %>%
  #   dplyr::select(-transcript_biotype, -n) %>%
  #   distinct(junID, .keep_all = T) %>%
  #   ungroup()
  # 
  # transcripts_hg_percentage_PC %>% 
  #   dplyr::count(percent)
  # 
  # # transcripts_hg_percentage_PC %>% filter(junID == "chr1:1290695-1290894:+")
  # # transcripts_hg_percentage_PC %>% filter(junID == "chr1:2583496-2584124:+")
  # 
  # 
  # logger::log_info(Sys.time(), " - getting lncRNA percentage per junction...")
  # 
  # transcripts_hg_percentage_lncRNA <- transcripts_hg_percentage %>% 
  #   dplyr::group_by(junID) %>%
  #   rowwise() %>%
  #   mutate(percent = ifelse (transcript_biotype != "protein_coding", percent, 0)) %>%
  #   filter(percent == max(percent)) %>%
  #   dplyr::group_by(junID) %>%
  #   filter(transcript_biotype == "lncRNA") %>%
  #   ungroup() %>%
  #   dplyr::select(-transcript_biotype, -n) %>%
  #   distinct(junID, .keep_all = T) %>%
  #   ungroup()
  # 
  # # transcripts_hg_percentage_lncRNA %>% filter(junID == "chr1:1290695-1290894:+")
  # # transcripts_hg_percentage_lncRNA %>% filter(junID == "chr1:2583496-2584124:+")
  # 
  # logger::log_info(Sys.time(), " - saving results...")
  # 
  # 
  # transcripts_hg_percentage_all <- plyr::rbind.fill(transcripts_hg_percentage_PC %>%
  #                                                     dplyr::rename(protein_coding = percent),
  #                                                   transcripts_hg_percentage_lncRNA %>%
  #                                                     dplyr::rename(lncRNA = percent)) %>% 
  #   replace(is.na(.), 0)
    

  # transcripts_hg_percentage_all %>% filter(junID == "chr1:1290695-1290894:+")
  # transcripts_hg_percentage_all %>% filter(junID == "chr1:2583496-2584124:+")
  
  # if (transcripts_hg_percentage_all %>% filter(protein_coding == 100) %>% distinct(lncRNA) %>% pull() != 0) {
  #   logger::log_info("ERROR! some only protein-coding introns have been also classified as lncRNAs!")
  #   break;
  # }
  # if (transcripts_hg_percentage_all %>% filter(lncRNA == 100) %>% distinct(protein_coding) %>% pull() != 0) {
  #   logger::log_info("ERROR! some only lncRNA introns have been also classified as protein-coding!")
  #   break;
  # }
  # 
  # logger::log_info(object.size(transcripts_hg_percentage_all), units = "Gb")
  
  transcripts_hg_percentage_tidy <- transcripts_hg_percentage %>%
    dplyr::select(-n) %>%
    #filter(junID == "chr1:10000015-10007569:+") %>%
    spread(transcript_biotype,percent) %>%
    replace(is.na(.), 0)
  
  saveRDS(object = transcripts_hg_percentage_tidy, #transcripts_hg_percentage_all,
          file = paste0(results.folder, "/all_split_reads_qc_level1_PC_biotype.rds"))
  
  logger::log_info(paste0(Sys.time(), " - results saved!"))
  
  
  # ## Calculate the biotype percentage
  # df_all_percentage <- df_all_junctions_tx_hg %>% #filter(junID == "chr1:9997332-10004092:+") %>%
  #   #filter(junID %in% junID_OK) %>%
  #   dplyr::group_by(junID, transcript_biotype) %>%
  #   distinct(tx_id_junction, .keep_all = T) %>% 
  #   summarise(n = n()) %>% 
  #   mutate(percent = (n / sum(n)) * 100) %>%
  #   ungroup()
  # 
  # 
  # df_all_percentage %>% filter(junID == "chr1:9997332-10007569:+")
  # df_all_junctions_tx_hg %>% filter(junID == "chr1:9997332-10007569:+")
  # 
  # # df_all_junctions_tx_hg %>% filter(junID == "chr1:100022245-100022385:+") 
  # # df_all_percentage %>% head()
  # # df_all_percentage %>%
  # #   filter(percent < 100, percent > 50)
  # 
  # saveRDS(object = df_all_percentage,
  #         file = paste0(results.folder, "/all_split_reads_qc_level2_biotype.rds"))
  
  ## Only filter by the protein-coding biotype
  # df_all_percentage_tidy_PC <- df_all_percentage %>% 
  #   dplyr::group_by(junID) %>%
  #   rowwise() %>%
  #   mutate(percent = ifelse (transcript_biotype == "protein_coding", percent, 0)) %>%
  #   ungroup() %>%
  #   dplyr::group_by(junID) %>%
  #   filter(percent == max(percent)) %>%
  #   dplyr::select(-transcript_biotype, -n) %>%
  #   distinct(junID, .keep_all = T) %>%
  #   ungroup()
  
  ## Only filter by the lncRNA biotype
  # df_all_percentage_tidy_lncRNA <- df_all_percentage %>% 
  #   dplyr::group_by(junID) %>%
  #   rowwise() %>%
  #   mutate(percent = ifelse (transcript_biotype == "lncRNA", percent, 0)) %>%
  #   ungroup() %>%
  #   dplyr::group_by(junID) %>%
  #   filter(percent == max(percent)) %>%
  #   dplyr::select(-transcript_biotype, -n) %>%
  #   distinct(junID, .keep_all = T) %>%
  #   ungroup()
  
  
  # df_all_percentage_tidy_merged <- merge(x = df_all_percentage_tidy_PC %>% dplyr::rename(percent_PC = percent),
  #                                        y = df_all_percentage_tidy_lncRNA %>% dplyr::rename(percent_lncRNA = percent),
  #                                        by = "junID")
  # df_all_percentage_tidy_merged %>% nrow()
  
  # if (df_all_percentage_tidy_merged %>% filter(percent_PC == 100) %>% distinct(percent_lncRNA) %>% pull() != 0) {
  #   logger::log_info("ERROR! some only protein-coding introns have been also classified as lncRNAs!")
  # }
  # if (df_all_percentage_tidy_merged %>% filter(percent_lncRNA == 100) %>% distinct(percent_PC) %>% pull() != 0) {
  #   logger::log_info("ERROR! some only lncRNA introns have been also classified as protein-coding!")
  # }
  # 
  # logger::log_info(object.size(df_all_percentage_tidy_merged), units = "Gb")
  # 
  # saveRDS(object = df_all_percentage_tidy_merged,
  #         file = paste0(results.folder, "/all_split_reads_qc_level2_PC_biotype.rds"))
  # logger::log_info(paste0(Sys.time(), " - results saved!"))
  
  ##########################################
  ## FREE UP SOME MEMORY 
  ##########################################
  
  # rm(all_split_reads_details_all_tissues)
  # rm(df_all_junctions)
  # rm(df_all_junctions_tx_hg)
  # rm(homo_sapiens_hg)
  # rm(transcripts_hg)
  # gc()
  
}