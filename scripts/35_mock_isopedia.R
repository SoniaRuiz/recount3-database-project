# =================================================================
# Mock Isopedia SQLite database for RAB32 (gene_id ENSG00000118508),
# matching the exact table schema produced by 34_database_isopedia.R:
#   transcript, junction, transcript_junction, sample, transcript_sample
#
# Coordinates are taken directly from the real IntroVerse short-read
# data for RAB32 you pasted (chr6, + strand):
#   - annotated intron (ref_junID 230435): 146544122-146549463
#   - downstream annotated intron (ref_junID 25623, "never" mis-spliced):
#     146549742-146554455
#   - 4 novel donor variants of the first intron (novel_junID 40950,
#     276451, 256784, 17135)
#   - 1 novel acceptor variant of the first intron (novel_junID 156581)
#   - 2 "novel_combo" variants from the 'other' table (ref_junID 49131, 49132)
#
# SIMPLIFICATION: the real generation script derives txID from
# distinct(chr, transcript_start, transcript_end) BEFORE joining
# junctions, which would collapse multiple different splice patterns
# into one txID if they happen to share the same outer transcript span.
# For this mock, each isoform is instead given its own txID directly
# (all sharing the same placeholder transcript_start/end), since the
# point here is to test isoform-structure rendering, not reproduce
# that span-collapsing behaviour.
# =================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(DBI)
library(RSQLite)


run = FALSE

if (run) {
  ## ---- 1. Isoform definitions -----------------------------------------------
  ## One row per isoform: its "first intron" (the variable alternative
  ## splice-site choice) plus the label matching the real event it represents.
  
  gene_chr <- 6
  transcript_start_default <- 146525000L   # before the furthest upstream novel donor
  transcript_end_default   <- 146556000L   # after the downstream annotated intron
  
  downstream_intron <- c(start = 146549742L, end = 146554455L)  # ref_junID 25623, constitutive
  
  isoform_defs <- tribble(
    ~txID, ~label,                                    ~intron1_start, ~intron1_end,
    1L,    "canonical (annotated intron, ref 230435)", 146544122L,    146549463L,
    2L,    "novel donor (novel 40950, -16177bp)",       146527945L,    146549463L,
    3L,    "novel donor (novel 276451)",                146544291L,    146549463L,
    4L,    "novel donor (novel 256784)",                146544566L,    146549463L,
    5L,    "novel donor (novel 17135)",                 146545057L,    146549463L,
    6L,    "novel acceptor (novel 156581, short)",       146544122L,    146544577L,
    7L,    "novel combo A ('other' ref 49131)",          146544094L,    146549463L,
    8L,    "novel combo B ('other' ref 49132)",          146544122L,    146549468L
  )
  
  ## ---- 2. transcript table ---------------------------------------------------
  
  transcript_tbl <- isoform_defs %>%
    transmute(
      txID,
      chr = gene_chr,
      transcript_start = transcript_start_default,
      transcript_end   = transcript_end_default
    )
  
  ## ---- 3. junction_long: two junctions per isoform (variable first intron +
  ## constitutive downstream intron) ------------------------------------------
  
  junction_long <- bind_rows(
    isoform_defs %>% transmute(txID, chr = gene_chr, start = intron1_start, end = intron1_end),
    isoform_defs %>% transmute(txID, chr = gene_chr,
                               start = downstream_intron["start"],
                               end   = downstream_intron["end"])
  ) %>%
    arrange(txID, start) %>%
    ## IMPORTANT: query_isopedia() looks up WHERE start = (typed_start - 1),
    ## i.e. it expects the junction table to store a 0-based start. The
    ## coordinates in isoform_defs above are the real, 1-based IntroVerse
    ## coordinates, so shift start down by 1 here before writing to the DB.
    mutate(start = start - 1L)
  
  ## ---- 4. junction table (unique junctions across all isoforms) -------------
  
  junction_tbl <- junction_long %>%
    distinct(chr, start, end) %>%
    mutate(junctionID = row_number())
  
  ## ---- 5. transcript_junction bridge table -----------------------------------
  
  transcript_junction_tbl <- junction_long %>%
    left_join(junction_tbl, by = c("chr", "start", "end")) %>%
    dplyr::select(txID, junctionID)
  
  ## ---- 6. sample table -------------------------------------------------------
  ## Mock samples spanning a few tissues, projects, and platforms for variety.
  
  ## NOTE: app.R's downstream code (output$queryIsopedia) does
  ## dplyr::select(-descrption) and renames sampleName/sample/date/source/
  ## tissue_subtype/cell_line -- so every one of those columns must exist here,
  ## even if only NA, or that select()/rename() chain will fail exactly like
  ## the error you just hit.
  sample_tbl <- tribble(
    ~sampleName,       ~sample,        ~tissue,              ~tissue_subtype, ~project,               ~platform,          ~disease,  ~cell_line,     ~date,        ~source,    ~descrption,
    "RAB32_ISO_S01",   "SAMN_MOCK01",  "Brain - Cortex",     NA_character_,   "IsopediaBrainAtlas",   "ONT PromethION",   "Control", NA_character_,  "2024-03-01", "Isopedia", "mock data for testing",
    "RAB32_ISO_S02",   "SAMN_MOCK02",  "Brain - Cortex",     NA_character_,   "IsopediaBrainAtlas",   "ONT PromethION",   "Control", NA_character_,  "2024-03-01", "Isopedia", "mock data for testing",
    "RAB32_ISO_S03",   "SAMN_MOCK03",  "Brain - Cerebellum", NA_character_,   "IsopediaBrainAtlas",   "ONT PromethION",   "Control", NA_character_,  "2024-03-01", "Isopedia", "mock data for testing",
    "RAB32_ISO_S04",   "SAMN_MOCK04",  "Liver",              NA_character_,   "IsopediaTissueSurvey", "PacBio Sequel II", "Control", NA_character_,  "2024-03-01", "Isopedia", "mock data for testing",
    "RAB32_ISO_S05",   "SAMN_MOCK05",  "Liver",              NA_character_,   "IsopediaTissueSurvey", "PacBio Sequel II", "Control", NA_character_,  "2024-03-01", "Isopedia", "mock data for testing",
    "RAB32_ISO_S06",   "SAMN_MOCK06",  "Muscle",             NA_character_,   "IsopediaTissueSurvey", "ONT PromethION",   "Control", NA_character_,  "2024-03-01", "Isopedia", "mock data for testing",
    "RAB32_ISO_S07",   "SAMN_MOCK07",  "K562",               NA_character_,   "IsopediaCellLines",    "ONT PromethION",   NA_character_, "K562",     "2024-03-01", "Isopedia", "mock data for testing",
    "RAB32_ISO_S08",   "SAMN_MOCK08",  "HepG2",              NA_character_,   "IsopediaCellLines",    "PacBio Sequel II", NA_character_, "HepG2",    "2024-03-01", "Isopedia", "mock data for testing"
  ) %>%
    mutate(sampleID = row_number(), .before = 1)
  
  ## ---- 7. transcript_sample: which samples support which isoform, with counts
  
  tx_sample_defs <- tribble(
    ~txID, ~sampleName,      ~totalEvidence,
    1L,    "RAB32_ISO_S01",  4L,
    1L,    "RAB32_ISO_S02",  3L,
    2L,    "RAB32_ISO_S03",  2L,
    3L,    "RAB32_ISO_S04",  1L,
    4L,    "RAB32_ISO_S05",  2L,
    5L,    "RAB32_ISO_S06",  3L,
    6L,    "RAB32_ISO_S07",  1L,
    7L,    "RAB32_ISO_S01",  2L,   # same sample can support >1 isoform, as in real data
    8L,    "RAB32_ISO_S08",  1L
  )
  
  transcript_sample_tbl <- tx_sample_defs %>%
    left_join(sample_tbl %>% dplyr::select(sampleID, sampleName), by = "sampleName") %>%
    dplyr::select(txID, sampleID, totalEvidence)
  
  ## ---- 8. Write to a mock SQLite file, matching the real table names --------
  
  mock_db_path <- "database/mock_isopedia_RAB32.sqlite"
  
  conn <- dbConnect(SQLite(), mock_db_path)
  dbWriteTable(conn, "transcript",          transcript_tbl,          overwrite = TRUE)
  dbWriteTable(conn, "junction",            junction_tbl,            overwrite = TRUE)
  dbWriteTable(conn, "transcript_junction", transcript_junction_tbl, overwrite = TRUE)
  dbWriteTable(conn, "sample",              sample_tbl,              overwrite = TRUE)
  dbWriteTable(conn, "transcript_sample",   transcript_sample_tbl,   overwrite = TRUE)
  dbDisconnect(conn)
  
  message("Mock Isopedia database written to: ", normalizePath(mock_db_path))
  
  ## ---- 9. How to test against query_isopedia() -------------------------------
  ## query_isopedia() has the real database path hardcoded inside the function
  ## body (see 01_backend_tab1.R), so for testing you have two options:
  ##
  ## Option A -- temporarily edit query_isopedia()'s dbConnect() call (or the
  ## ISOPEDIA_DATABASE constant in config.R, if your working copy uses that
  ## instead) to point at mock_isopedia_RAB32.sqlite, then call it normally:
  #
  # query_isopedia(chr = 6, start = 146544122, end = 146549463)
  # ## note: query_isopedia() filters on start = (start - 1). Since junction_tbl
  # ## now stores (true_start - 1) after the fix above, pass the TRUE/actual
  # ## junction start directly here (146544122, not +1 or -1) -- query_isopedia()
  # ## does the -1 itself internally, matching how jx_start is used elsewhere
  # ## in the app (e.g. the search box's Start field).
  #
  ## Option B -- if you'd rather not touch query_isopedia() at all, copy its
  ## body into a throwaway query_isopedia_mock(chr, start, end, db_path) that
  ## takes the db path as a parameter, for isolated testing:
  #
  # query_isopedia_mock <- function(chr, start, end, db_path) {
  #   conn <- dbConnect(RSQLite::SQLite(), db_path)
  #   junction_query <- paste0("SELECT junctionID FROM junction WHERE chr=", chr, " AND start=", (start - 1), " AND end=", end)
  #   db_junction <- DBI::dbGetQuery(conn, junction_query)
  #   if (nrow(db_junction) == 0) return(NULL)
  #   junction_bridge_tx_query <- paste0("SELECT * FROM transcript_junction WHERE junctionID IN (", paste(db_junction$junctionID, collapse=","), ")")
  #   db_junction_bridge_tx <- DBI::dbGetQuery(conn, junction_bridge_tx_query)
  #   tx_junction_query <- paste0("SELECT * FROM transcript_junction WHERE txID IN (", paste(db_junction_bridge_tx$txID, collapse=","), ")")
  #   db_tx_junction <- DBI::dbGetQuery(conn, tx_junction_query)
  #   tx_junction_coord_query <- paste0("SELECT * FROM junction WHERE junctionID IN (", paste(db_tx_junction$junctionID, collapse=","), ")")
  #   db_tx_junction_coord <- DBI::dbGetQuery(conn, tx_junction_coord_query)
  #   db_tx_junx_coord <- db_tx_junction %>% left_join(db_tx_junction_coord, by = "junctionID")
  #   tx_sample_query <- paste0("SELECT * FROM transcript_sample WHERE txID IN (", paste(unique(db_tx_junx_coord$txID), collapse=","), ")")
  #   db_tx_sample <- DBI::dbGetQuery(conn, tx_sample_query)
  #   db_tx_junx_coord_w_sample <- db_tx_junx_coord %>% left_join(db_tx_sample, by = "txID")
  #   tx_query <- paste0("SELECT txID, transcript_start, transcript_end FROM transcript WHERE txID IN (", paste(unique(db_tx_junx_coord_w_sample$txID), collapse=","), ")")
  #   db_tx <- DBI::dbGetQuery(conn, tx_query)
  #   db_tx_junx_coord_w_sample_w_tx <- db_tx_junx_coord_w_sample %>% left_join(db_tx, by = "txID")
  #   sample_query <- paste0("SELECT * FROM sample WHERE sampleID IN (", paste(unique(db_tx_junx_coord_w_sample_w_tx$sampleID), collapse=","), ")")
  #   db_sample <- DBI::dbGetQuery(conn, sample_query)
  #   db_tx_junx_coord_w_sample_w_tx %>%
  #     left_join(db_sample, by = "sampleID") %>%
  #     dplyr::select(-c(sampleID, junctionID))
  # }
  #
  # query_isopedia_mock(chr = 6, start = 146544122, end = 146549463, db_path = "mock_isopedia_RAB32.sqlite")
}
