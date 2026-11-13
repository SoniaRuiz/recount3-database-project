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

chunk_dir     <- "/rds/user/sg2173/hpc-work/moved_from_home_dir/PROJECTS/srRNAseq/introverse_v2/database/bed/chunks"
database_file <- "/rds/user/sg2173/hpc-work/moved_from_home_dir/PROJECTS/srRNAseq/introverse_v2/database/isopedia_test.sqlite"
batch_size    <- 2000   # rows read/pivoted at a time

# For a timing test: point this at a small file and use a different database_file
chunk_files_override <- "/rds/user/sg2173/hpc-work/moved_from_home_dir/PROJECTS/srRNAseq/introverse_v2/database/bed/chunks/chunk_test.f1.tsv.gz"

# Layout confirmed from chunk_000.tsv.gz:
#   lines 1-1008 : ##[SAMPLE] header + 1007 samples
#   line 1009    : ##[DBINFO]
#   line 1010    : data header; sample columns start at column 17
n_samples  <- 1007
skip_lines <- 1009
first_sample_col <- 17

if (is.null(chunk_files_override)) {
  chunk_files <- sort(list.files(chunk_dir, pattern = "^chunk_[0-9]+\\.f1\\.tsv\\.gz$", full.names = TRUE))
} else {
  chunk_files <- chunk_files_override
}
stopifnot(length(chunk_files) > 0)

# ------------------------------------------------------------------------------
# Sample table (read once from the first file; assumes identical across chunks)
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
# Database schema (UNIQUE / PRIMARY KEY constraints do the deduplication)
# ------------------------------------------------------------------------------

conn <- dbConnect(SQLite(), database_file)
dbExecute(conn, "PRAGMA synchronous = OFF")
dbExecute(conn, "PRAGMA temp_store = MEMORY")

for (tbl in c("transcript", "junction", "transcript_junction", "transcript_sample", "sample")) {
  if (dbExistsTable(conn, tbl)) dbRemoveTable(conn, tbl)
}
dbWriteTable(conn, "sample", sample_tbl)

dbExecute(conn, '
  CREATE TABLE transcript (
    txID INTEGER PRIMARY KEY,
    chr TEXT,
    transcript_start INTEGER,
    transcript_end INTEGER,
    UNIQUE (chr, transcript_start, transcript_end)
  )')
dbExecute(conn, '
  CREATE TABLE junction (
    junctionID INTEGER PRIMARY KEY,
    chr TEXT,
    "start" INTEGER,
    "end" INTEGER,
    UNIQUE (chr, "start", "end")
  )')
dbExecute(conn, '
  CREATE TABLE transcript_junction (
    txID INTEGER,
    junctionID INTEGER,
    PRIMARY KEY (txID, junctionID)
  )')
dbExecute(conn, '
  CREATE TABLE transcript_sample (
    txID INTEGER,
    sampleID INTEGER,
    totalEvidence INTEGER,
    PRIMARY KEY (txID, sampleID)
  )')

# ------------------------------------------------------------------------------
# Process one batch of rows
# ------------------------------------------------------------------------------

state <- new.env()
state$batches <- 0L
state$start_time <- Sys.time()

process_batch <- function(df, pos) {

  df <- df %>%
    rename(chr = chr1,
           transcript_start = start_pos_left,
           transcript_end = end_pos_right) %>%
    mutate(transcript_start = as.integer(transcript_start),
           transcript_end = as.integer(transcript_end))

  stage_tx <- df %>% distinct(chr, transcript_start, transcript_end)

  stage_junc <- df %>%
    select(chr, transcript_start, transcript_end, splice_junctions) %>%
    separate_rows(splice_junctions, sep = ",") %>%
    separate(splice_junctions, into = c("jstart", "jend"), sep = "-") %>%
    mutate(jstart = as.integer(jstart), jend = as.integer(jend)) %>%
    distinct()

  stage_ts <- df %>%
    select(chr, transcript_start, transcript_end, all_of(sample_cols)) %>%
    pivot_longer(cols = all_of(sample_cols),
                 names_to = "sampleName",
                 values_to = "value") %>%
    filter(!is.na(value), value != ".") %>%
    mutate(count = as.integer(str_extract(value, "^[0-9]+"))) %>%
    filter(!is.na(count), count > 0) %>%
    inner_join(sample_tbl %>% select(sampleID, sampleName), by = "sampleName") %>%
    select(chr, transcript_start, transcript_end, sampleID, count)

  dbBegin(conn)

  dbWriteTable(conn, "stage_tx",   stage_tx,   temporary = TRUE, overwrite = TRUE)
  dbWriteTable(conn, "stage_junc", stage_junc, temporary = TRUE, overwrite = TRUE)
  dbWriteTable(conn, "stage_ts",   stage_ts,   temporary = TRUE, overwrite = TRUE)

  dbExecute(conn, '
    INSERT OR IGNORE INTO transcript (chr, transcript_start, transcript_end)
    SELECT chr, transcript_start, transcript_end FROM stage_tx')

  dbExecute(conn, '
    INSERT OR IGNORE INTO junction (chr, "start", "end")
    SELECT DISTINCT chr, jstart, jend FROM stage_junc')

  dbExecute(conn, '
    INSERT OR IGNORE INTO transcript_junction (txID, junctionID)
    SELECT t.txID, j.junctionID
    FROM stage_junc s
    JOIN transcript t ON t.chr = s.chr
                     AND t.transcript_start = s.transcript_start
                     AND t.transcript_end = s.transcript_end
    JOIN junction j ON j.chr = s.chr AND j."start" = s.jstart AND j."end" = s.jend')

  dbExecute(conn, '
    INSERT OR IGNORE INTO transcript_sample (txID, sampleID, totalEvidence)
    SELECT t.txID, s.sampleID, s.count
    FROM stage_ts s
    JOIN transcript t ON t.chr = s.chr
                     AND t.transcript_start = s.transcript_start
                     AND t.transcript_end = s.transcript_end')

  dbCommit(conn)

  state$batches <- state$batches + 1L
  if (state$batches %% 50 == 0) {
    message(format(Sys.time()), "  batches: ", state$batches,
            " | elapsed (min): ", round(as.numeric(difftime(Sys.time(), state$start_time, units = "mins")), 1))
  }
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# Loop over chunk files
# ------------------------------------------------------------------------------

for (f in chunk_files) {
  message(format(Sys.time()), "  reading ", basename(f))
  read_tsv_chunked(
    f,
    callback = SideEffectChunkCallback$new(process_batch),
    chunk_size = batch_size,
    skip = skip_lines,
    col_types = cols(.default = col_character()),
    progress = FALSE
  )
  message(format(Sys.time()), "  done ", basename(f))
}

message("transcripts: ", dbGetQuery(conn, "SELECT COUNT(*) FROM transcript")[[1]],
        " | junctions: ", dbGetQuery(conn, "SELECT COUNT(*) FROM junction")[[1]],
        " | transcript_sample rows: ", dbGetQuery(conn, "SELECT COUNT(*) FROM transcript_sample")[[1]])

dbDisconnect(conn)