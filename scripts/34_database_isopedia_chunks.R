# module load R/4.3.1-icelake
# module load gcc/11
# module load libiconv/1.17/gcc/rlhw7p3w
# module load curl/8.7.1/gcc/d76pzrod

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(DBI)
library(RSQLite)

options(width = 200)

# ------------------------------------------------------------------------------
# Paths / settings
# ------------------------------------------------------------------------------

run = FALSE

if (run) {
  
  chunk_dir     <- "/rds/user/sg2173/hpc-work/moved_from_home_dir/PROJECTS/srRNAseq/introverse_v2/database/bed/chunks"
  database_file <- "/rds/user/sg2173/hpc-work/moved_from_home_dir/PROJECTS/srRNAseq/introverse_v2/database/isopedia_large.sqlite"  # <-- set where you want the database
  batch_size    <- 2000   # rows read/pivoted at a time
  
  # Layout confirmed from chunk_000.tsv.gz:
  #   line 1            : ##[SAMPLE] column header
  #   lines 2-1008      : 1007 ##[SAMPLE] rows
  #   line 1009         : ##[DBINFO]
  #   line 1010         : data header (#id, chr1, ..., FORMAT, <sample columns>)
  #   sample columns start at column 17
  n_samples  <- 1007
  skip_lines <- 1009
  first_sample_col <- 17
  
  chunk_files <- sort(list.files(chunk_dir, pattern = "^chunk_[0-9]+\\.tsv\\.gz$", full.names = TRUE))
  stopifnot(length(chunk_files) > 0)
  
  # ------------------------------------------------------------------------------
  # Sample table (read once from the first chunk; assumes all chunks share the same
  # sample metadata block, as they were produced from the same index)
  # ------------------------------------------------------------------------------
  
  sample_metadata_isopedia <- read_tsv(
    chunk_files[1],
    n_max = n_samples,
    col_types = cols(.default = col_character())
  ) %>%
    rename(Sample = `##[SAMPLE]Sample`) %>%
    mutate(Sample = str_remove(Sample, fixed("##[SAMPLE]")))
  
  header <- read_tsv(
    chunk_files[1],
    skip = skip_lines,
    n_max = 0,
    col_types = cols(.default = col_character())
  )
  sample_cols <- colnames(header)[first_sample_col:ncol(header)]
  
  sample_tbl <- tibble(
    sampleID = seq_along(sample_cols),
    sampleName = sample_cols
  ) %>%
    inner_join(y = sample_metadata_isopedia, by = c("sampleName" = "Sample")) %>%
    dplyr::select_if(~ !all(is.na(.)))
  
  # ------------------------------------------------------------------------------
  # Database + global ID state (shared across batches and chunks)
  # ------------------------------------------------------------------------------
  
  conn <- dbConnect(SQLite(), database_file)
  
  for (tbl in c("transcript", "junction", "transcript_junction", "transcript_sample")) {
    if (dbExistsTable(conn, tbl)) dbRemoveTable(conn, tbl)
  }
  dbWriteTable(conn, "sample", sample_tbl, overwrite = TRUE)
  
  state <- new.env()
  state$next_tx <- 1L
  state$next_junction <- 1L
  
  tx_ids       <- new.env(hash = TRUE)  # "chr:start:end"             -> txID
  junction_ids <- new.env(hash = TRUE)  # "chr:start:end"             -> junctionID
  seen_rows    <- new.env(hash = TRUE)  # "chr|start|end|junctions"   -> 1
  
  lookup_ids <- function(keys, env) {
    vapply(keys, function(k) {
      v <- env[[k]]
      if (is.null(v)) NA_integer_ else v
    }, integer(1), USE.NAMES = FALSE)
  }
  
  # ------------------------------------------------------------------------------
  # Process one batch of rows
  # ------------------------------------------------------------------------------
  
  process_batch <- function(df, pos) {
    
    df <- df %>%
      rename(chr = chr1,
             transcript_start = start_pos_left,
             transcript_end = end_pos_right) %>%
      mutate(transcript_start = as.integer(transcript_start),
             transcript_end = as.integer(transcript_end))
    
    # Drop rows already seen (same isoform matched by several query junctions,
    # within this batch or in an earlier batch/chunk)
    row_key <- paste(df$chr, df$transcript_start, df$transcript_end, df$splice_junctions, sep = "|")
    keep <- !duplicated(row_key) & is.na(lookup_ids(row_key, seen_rows))
    if (!any(keep)) return(invisible(NULL))
    
    df <- df[keep, ]
    row_key <- row_key[keep]
    list2env(setNames(as.list(rep(1L, length(row_key))), row_key), envir = seen_rows)
    
    dbBegin(conn)
    
    # --- transcript table ---
    tx_key <- paste(df$chr, df$transcript_start, df$transcript_end, sep = ":")
    tx_id <- lookup_ids(tx_key, tx_ids)
    new_tx <- which(is.na(tx_id) & !duplicated(tx_key))
    if (length(new_tx) > 0) {
      new_ids <- seq.int(state$next_tx, length.out = length(new_tx))
      state$next_tx <- state$next_tx + length(new_tx)
      list2env(setNames(as.list(new_ids), tx_key[new_tx]), envir = tx_ids)
      dbWriteTable(conn, "transcript",
                   tibble(txID = new_ids,
                          chr = df$chr[new_tx],
                          transcript_start = df$transcript_start[new_tx],
                          transcript_end = df$transcript_end[new_tx]),
                   append = TRUE)
      tx_id <- lookup_ids(tx_key, tx_ids)
    }
    df$txID <- tx_id
    
    # --- junction + transcript_junction tables ---
    junction_long <- df %>%
      select(txID, chr, splice_junctions) %>%
      separate_rows(splice_junctions, sep = ",") %>%
      separate(splice_junctions, into = c("start", "end"), sep = "-") %>%
      mutate(start = as.integer(start),
             end = as.integer(end))
    
    j_key <- paste(junction_long$chr, junction_long$start, junction_long$end, sep = ":")
    j_id <- lookup_ids(j_key, junction_ids)
    new_j <- which(is.na(j_id) & !duplicated(j_key))
    if (length(new_j) > 0) {
      new_jids <- seq.int(state$next_junction, length.out = length(new_j))
      state$next_junction <- state$next_junction + length(new_j)
      list2env(setNames(as.list(new_jids), j_key[new_j]), envir = junction_ids)
      dbWriteTable(conn, "junction",
                   tibble(chr = junction_long$chr[new_j],
                          start = junction_long$start[new_j],
                          end = junction_long$end[new_j],
                          junctionID = new_jids),
                   append = TRUE)
      j_id <- lookup_ids(j_key, junction_ids)
    }
    dbWriteTable(conn, "transcript_junction",
                 tibble(txID = junction_long$txID, junctionID = j_id),
                 append = TRUE)
    
    # --- transcript_sample table (read counts) ---
    tx_sample <- df %>%
      select(txID, all_of(sample_cols)) %>%
      pivot_longer(cols = all_of(sample_cols),
                   names_to = "sampleName",
                   values_to = "value") %>%
      filter(!is.na(value)) %>%
      mutate(count = as.integer(str_extract(value, "^[0-9]+"))) %>%
      filter(!is.na(count), count > 0) %>%
      left_join(sample_tbl %>% select(sampleID, sampleName), by = "sampleName") %>%
      select(txID, sampleID, totalEvidence = count)
    
    dbWriteTable(conn, "transcript_sample", tx_sample, append = TRUE)
    
    dbCommit(conn)
    invisible(NULL)
  }
  
  # ------------------------------------------------------------------------------
  # Loop over chunks
  # ------------------------------------------------------------------------------
  
  for (f in chunk_files) {
    message(Sys.time(), "  reading ", basename(f))
    read_tsv_chunked(
      f,
      callback = SideEffectChunkCallback$new(process_batch),
      chunk_size = batch_size,
      skip = skip_lines,
      col_types = cols(.default = col_character()),
      progress = FALSE
    )
    message(Sys.time(), "  done ", basename(f),
            " | transcripts so far: ", state$next_tx - 1L,
            " | junctions so far: ", state$next_junction - 1L)
  }
  
  dbDisconnect(conn)
}
