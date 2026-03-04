# module load R/4.3.1-icelake
# module load gcc/11
# module load libiconv/1.17/gcc/rlhw7p3w
# module load curl/8.7.1/gcc/d76pzrod  

# options(width=150)

library(readr)
library(dplyr)
library(stringr)
library(DBI)
library(RSQLite)
library(tidyr)
library(data.table)


run = FALSE

if (run) {
  options(width = 200)
  
  # BiocManager::install("SGSeq")
  # source("/home/sruiz/POST_DOC/recount3-database-project/scripts/34_database_isopedia.R")
  
  
  isopedia_results <- "/home/sruiz/POST_DOC/recount3-database-project/database/bed/isopedia_introverse_splice_junctions.tsv.gz"
  database_path <- "/home/sruiz/POST_DOC/recount3-database-project/database/isopedia"
  
  # ------------------------------------------------------------------------------
  # Load sample metadata from Isopedia
  # ------------------------------------------------------------------------------
  
  sample_metadata_isopedia <- read_tsv(isopedia_results, n_max = 1008, show_col_types = FALSE)
  
  # Remove '##[SAMPLE]' string on the first column
  sample_metadata_isopedia <- sample_metadata_isopedia %>%
    dplyr::rename(Sample = `##[SAMPLE]Sample`) %>%
    mutate(Sample = str_remove(Sample, fixed("##[SAMPLE]")))
  
  # ------------------------------------------------------------------------------
  # Load junction data from Isopedia
  # ------------------------------------------------------------------------------
  
  ## Load junction data from Isopedia that matches information on IntroVerse
  df <- fread(cmd = paste("zcat", isopedia_results), skip = 1009, sep = "\t")
  df %>% distinct(dist_to_matched_sj)
  
  #df = df %>% filter(dist_to_matched_sj == "0,0" | dist_to_matched_sj == "-1,0")
  
  # Identify sample columns
  sample_cols <- colnames(df)[17:ncol(df)]
  
  # grep("^sample_", colnames(df), value = TRUE)
  
  
  # Step 2 - Build transcript table ----------------------------------------------
  
  transcript_tbl <- df %>%
    dplyr::select(chr = chr1, 
                  transcript_start = start_pos_left,
                  transcript_end = end_pos_right) %>%
    distinct() %>%
    mutate(txID = row_number()) %>%
    dplyr::select(txID, chr, transcript_start, transcript_end)
  transcript_tbl
  
  
  # Step 3 - Extract junctions and keep transcript identity ----------------------
  
  junction_long <- df %>%
    dplyr::select(chr = chr1,
                  transcript_start = start_pos_left,
                  transcript_end = end_pos_right,
                  splice_junctions) %>%
    separate_rows(splice_junctions, sep = ",") %>%
    separate(splice_junctions, into = c("start", "end"), sep = "-") %>%
    mutate(start = as.integer(start),
           end = as.integer(end))
  # rm(df)
  gc()
  
  # Step 4: Link junctions to transcripts ----------------------------------------
  
  junction_long <- junction_long %>% left_join(transcript_tbl, by = c("chr", "transcript_start", "transcript_end"))
  junction_long
  
  
  # Step 5: Create unique junction table -----------------------------------------
  
  junction_tbl <- junction_long %>%
    distinct(chr, start, end) %>%
    mutate(junctionID = row_number())
  
  
  # Step 6: Create N:N table bridge between transcripts and junctions ------------
  
  transcript_junction_tbl <- junction_long %>%
    left_join(junction_tbl, by = c("chr", "start", "end")) %>%
    dplyr::select(txID, junctionID)
  
  
  # Step 7: Create sample_longread table -----------------------------------------
  # sample_cols <- grep("^sample_", colnames(df), value = TRUE)
  
  sample_tbl <- tibble(
    sampleID = seq_along(sample_cols),
    sampleName = sample_cols
  ) %>%
    inner_join(y = sample_metadata_isopedia, by = c("sampleName"="Sample")) %>%
    dplyr::select_if(~ !all(is.na(.)))
  
  
  
  
  # ------------------------------------------------------------------------------
  # Step 8: bridge table N:N between sampleIDs and TxIDs to provide the read counts
  # ------------------------------------------------------------------------------
  
  # Work on the sample columns as a matrix, not a pivoted tibble
  sample_mat <- as.matrix(df[, ..sample_cols]) # Equivalent to: sample_mat <- as.matrix(df[, sample_cols, with = FALSE])
  
  # Find only the non-NA cells directly (row/col indices), skipping the huge pivot
  idx <- which(!is.na(sample_mat), arr.ind = TRUE)
  nrow(idx)
  
  rm(df)
  gc()
  
  conn <- dbConnect(SQLite(), file.path(database_path, "isopedia.sqlite"))
  
  transcript_dt <- as.data.table(transcript_tbl)
  transcript_dt[, row := .I]
  
  sample_dt <- as.data.table(sample_tbl)
  
  chunk_size <- 150   # number of sample columns per batch — adjust based on RAM
  n_samples  <- ncol(sample_mat)
  
  dbExecute(conn, "DROP TABLE IF EXISTS transcript_sample")
  
  for (start_col in seq(1, n_samples, by = chunk_size)) {
    
    end_col  <- min(start_col + chunk_size - 1, n_samples)
    cols_idx <- start_col:end_col
    
    sub_mat <- sample_mat[, cols_idx, drop = FALSE]
    idx_sub <- which(!is.na(sub_mat), arr.ind = TRUE)
    
    if (nrow(idx_sub) == 0) next
    
    dt <- data.table(
      row        = idx_sub[, "row"],
      sampleName = sample_cols[cols_idx][idx_sub[, "col"]],
      count       = tstrsplit(sub_mat[idx_sub], ":", fixed = TRUE, keep = 1L)[[1]]
    )
    
    dt[, count := as.integer(count)]
    dt <- dt[count > 0]
    
    dt <- transcript_dt[dt, on = "row"]
    dt <- sample_dt[dt, on = "sampleName"]
    
    result <- dt[, .(txID, sampleID, totalEvidence = count)]
    
    dbWriteTable(conn, "transcript_sample", result, append = TRUE)
    
    rm(sub_mat, idx_sub, dt, result)
    gc()
  }
  
  dbDisconnect(conn)
  
  tx_sample_tbl <- df %>%
    dplyr::select(chr = chr1,
                  transcript_start = start_pos_left,
                  transcript_end = end_pos_right,
                  all_of(sample_cols)) %>%
    left_join(transcript_tbl,
              by = c("chr", "transcript_start", "transcript_end")) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sampleName",
      values_to = "value"
    ) %>%
    separate(value, into = c("count", "cpm", "coords"), sep = ":", fill = "right") %>%
    mutate(count = as.integer(count)) %>%
    left_join(sample_tbl, by = "sampleName") %>%
    filter(!is.na(count) & count > 0) %>%
    dplyr::select(txID, sampleID, totalEvidence = count)
  
  # ------------------------------------------------------------------------------
  # Step 9: write to SQL
  # ------------------------------------------------------------------------------
  
  conn <- dbConnect(SQLite(), paste0(database_path, "isopedia.sqlite"))
  dbWriteTable(conn, "transcript", transcript_tbl, overwrite = TRUE)
  dbWriteTable(conn, "junction", junction_tbl, overwrite = TRUE)
  dbWriteTable(conn, "transcript_junction", transcript_junction_tbl, overwrite = TRUE)
  dbWriteTable(conn, "sample", sample_tbl, overwrite = TRUE)
  dbWriteTable(conn, "transcript_sample", tx_sample_tbl, overwrite = TRUE)
  
  
}
