main_path <- normalizePath(path = "./")
source(paste0(main_path, "/testthat/helper_Files/helper_Global.R"))
skip_if(!test_TranscriptTable, "Transcript table tests not executed. Variable test_TranscriptTable set to FALSE in global options.")
source(paste0(main_path, "/testthat/helper_Files/helper_TranscriptTable.R"))

context("\tTest that transcript exists in the reference transcriptome")
test_that("Test that transcript exists in the reference transcriptome", {
  reference_transcripts <- hg38 %>% as_tibble() %>% select(transcript_id) %>% distinct()
  df_combined <- df_transcript %>% left_join(reference_transcripts, by = "transcript_id")
  
  ## All transcript IDs are found in the reference transcriptome
  expect_true(all(df_transcript$transcript_id %in% reference_transcripts$transcript_id))
  
})

context("\tTest reference between transcript and gene table")
test_that("Test reference between transcript and gene table", {
  gene_id <- df_transcript %>% pull(gene_id)
  
  ## All `transcript_id`` are properly referenced to a row in the gene table
  expect_true(all(gene_id %in% df_gene$id))
})


context("\tTest that transcript maximum TSL number is correct")
test_that("Test that transcript maximum TSL number is correct", {
  
  hg38_TSL <- hg38 %>% 
    as_tibble() %>%
    filter(type == "transcript") %>%
    dplyr::select(transcript_id, "transcript_support_level") %>%
    mutate(transcript_support_level = str_sub(string = transcript_support_level,
                                              start = 1,
                                              end = 2))
  
  df_combined <- df_transcript %>% left_join(hg38_TSL, by = "transcript_id") %>%
    mutate(transcript_support_level = ifelse(is.na(transcript_support_level),10,transcript_support_level)) %>%
    mutate(transcript_support_level = ifelse(transcript_support_level == "NA",10,transcript_support_level)) %>%
    mutate(transcript_support_level = transcript_support_level %>% as.integer())
  
  ## All gene widths are consistent with the reference transcriptome
  expect_equal(df_combined$TSL, df_combined$transcript_support_level)
})

context("\tTest that transcript MANE is correct")
test_that("Test that transcript MANE is correct", {
  
  hg_mane_tidy <- hg_mane %>% 
    as_tibble() %>%
    filter(type == "transcript") %>%
    mutate(transcript_id = str_sub(string = transcript_id,
                                   start = 1,
                                   end = 15)) %>%
    dplyr::select(transcript_id, tag) %>%
    mutate(tag = ifelse(is.na(tag), 0,1))
  
  
  df_combined <- df_transcript %>% left_join(hg_mane_tidy, by = "transcript_id")%>%
    mutate(tag = ifelse(is.na(tag), 0,1))
  
  ## All transcripts that are MANE Select have are correctly set
  expect_equal(df_combined$MANE, df_combined$tag)
})

clearAllVariables()

