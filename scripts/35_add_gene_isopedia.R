## ------------------------------------------------------------------------------
## Build a GTF file from the Isopedia SQLite database
##
## Strand, gene_id, gene_name and transcript_biotype are NOT available in the
## Isopedia schema (transcript / junction / transcript_junction / sample /
## transcript_sample), so they are written as placeholders:
##   - strand, score, frame  -> "."   (standard GTF "unknown" placeholder)
##   - gene_id                -> the transcript's txID (as decided)
##   - transcript_id          -> "TX<txID>" (as decided)
## These are expected to be resolved downstream by ORFannotate.
## ------------------------------------------------------------------------------

library(DBI)
library(RSQLite)
library(dplyr)
library(tibble)

run = FALSE

if (run) {
  # ---- CONFIG -------------------------------------------------------------------
  
  ## Update this path if your local Isopedia database lives elsewhere
  database_path <- "/home/sruiz/POST_DOC/recount3-database-project/database/isopedia/isopedia.sqlite"
  output_gtf    <- "/home/sruiz/POST_DOC/recount3-database-project/database/isopedia/isopedia.gtf"
  gtf_source    <- "Isopedia"
  
  placeholder_score  <- "."
  placeholder_strand <- "."
  placeholder_frame  <- "."
  
  # ---- LOAD DATA FROM DATABASE ---------------------------------------------------
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  DBI::dbListTables(conn = con)
  
  transcript_tbl          <- dbGetQuery(con, "SELECT * FROM 'transcript'")
  transcript_junction_tbl <- dbGetQuery(con, "SELECT * FROM 'transcript_junction'")
  junction_tbl            <- dbGetQuery(con, "SELECT * FROM 'junction'")
  
  dbDisconnect(con)
  
  # ---- JOIN INTRONS (JUNCTIONS) TO THEIR TRANSCRIPTS -----------------------------
  
  introns_per_tx <- transcript_junction_tbl %>%
    inner_join(junction_tbl, by = "junctionID") %>%
    dplyr::select(txID, intron_start = start, intron_end = end)
  
  # ---- BUILD EXON COORDINATES PER TRANSCRIPT -------------------------------------
  # Exons are the genomic gaps between transcript_start/transcript_end and the
  # transcript's introns (junctions), ordered by ascending genomic position.
  # (Strand is unknown, so exon_number below reflects genomic, not necessarily
  # transcriptional, order.)
  
  build_exons <- function(tx_start, tx_end, introns) {
    if (nrow(introns) == 0) {
      return(tibble(exon_start = tx_start, exon_end = tx_end))
    }
    introns <- introns %>% arrange(intron_start)
    exon_start <- c(tx_start, introns$intron_end + 1)
    exon_end   <- c(introns$intron_start - 1, tx_end)
    tibble(exon_start = exon_start, exon_end = exon_end) %>%
      filter(exon_end >= exon_start)
  }
  
  # ---- WRITE GTF ------------------------------------------------------------------
  
  gtf_lines <- character()
  
  for (tx in transcript_tbl$txID) {
    print(tx)
    tx_row <- transcript_tbl %>% filter(txID == tx)
    
    gene_id       <- as.character(tx)
    transcript_id <- paste0("TX", tx)
    
    tx_introns <- introns_per_tx %>% filter(txID == tx)
    tx_exons   <- build_exons(tx_row$transcript_start, tx_row$transcript_end, tx_introns) %>%
      arrange(exon_start)
    
    ## Transcript feature line
    gtf_lines <- c(gtf_lines, paste(
      tx_row$chr, gtf_source, "transcript",
      tx_row$transcript_start, tx_row$transcript_end,
      placeholder_score, placeholder_strand, placeholder_frame,
      paste0('gene_id "', gene_id, '"; transcript_id "', transcript_id, '";'),
      sep = "\t"
    ))
    
    ## Exon feature lines
    for (i in seq_len(nrow(tx_exons))) {
      gtf_lines <- c(gtf_lines, paste(
        tx_row$chr, gtf_source, "exon",
        tx_exons$exon_start[i], tx_exons$exon_end[i],
        placeholder_score, placeholder_strand, placeholder_frame,
        paste0('gene_id "', gene_id, '"; transcript_id "', transcript_id,
               '"; exon_number "', i, '";'),
        sep = "\t"
      ))
    }
  }
  
  writeLines(gtf_lines, output_gtf)
  
  message("GTF written to: ", normalizePath(output_gtf))
  
  
  
  
  # # module load R/4.3.1-icelake
  # # module load gcc/11
  # # module load libiconv/1.17/gcc/rlhw7p3w
  # # module load curl/8.7.1/gcc/d76pzrod  
  # 
  # # options(width=150)
  # 
  # library(readr)
  # library(dplyr)
  # library(stringr)
  # library(DBI)
  # library(RSQLite)
  # library(tidyr)
  # library(data.table)
  # 
  # options(width = 200)
  # # BiocManager::install("SGSeq")
  # 
  # isopedia_database_path <- "/home/sruiz/POST_DOC/recount3-database-project/database/isopedia"
  # conn_isopedia <- dbConnect(SQLite(), file.path(isopedia_database_path, "isopedia.sqlite"))
  # DBI::dbListTables(conn_isopedia)
  # 
  # introverse_database_path <- "/home/sruiz/POST_DOC/recount3-database-project/database/SRP058181_1read_subsampleFALSE/115/"
  # conn_introverse <- dbConnect(SQLite(), file.path(introverse_database_path, "SRP058181_1read_subsampleFALSE.sqlite"))
  # 
  # 
  # # ------------------------------------------------------------------------------
  # # INTROVERSE: Get gene and junction data from introverse
  # # ------------------------------------------------------------------------------
  # 
  # introverse_introns <- DBI::dbGetQuery(conn_introverse, "SELECT intron.seqnames, intron.start, intron.end, gene.gene_name FROM intron 
  #                 INNER JOIN gene ON gene.gene_id = intron.gene_id") %>% as_tibble() %>%
  #   mutate(seqnames = str_remove(pattern = "chr", seqnames))
  # 
  # introverse_novel <- DBI::dbGetQuery(conn_introverse, "SELECT novel.seqnames, novel.start, novel.end, gene.gene_name FROM novel
  #                 INNER JOIN intron ON intron.ref_junID = novel.ref_junID
  #                 INNER JOIN gene ON gene.gene_id = intron.gene_id") %>% as_tibble() %>%
  #   mutate(seqnames = str_remove(pattern = "chr", seqnames))
  # 
  # introverse_combo <- DBI::dbGetQuery(conn_introverse, "SELECT other.seqnames, other.start, other.end, gene.gene_name FROM other
  #                 INNER JOIN gene ON gene.gene_id = other.gene_id") %>% as_tibble() %>%
  #   mutate(seqnames = str_remove(pattern = "chr", seqnames))
  # 
  # introverse_jxn = rbind(introverse_introns, introverse_combo, introverse_novel)
  # 
  # 
  # # ------------------------------------------------------------------------------
  # # ISOPEDIA: Load junction data from Isopedia
  # # ------------------------------------------------------------------------------
  # 
  # isopedia_junctions <- DBI::dbGetQuery(conn_isopedia, "SELECT * FROM junction") %>% as_tibble()
  # transcript_junction <- DBI::dbGetQuery(conn_isopedia, "SELECT * FROM transcript_junction") %>% as_tibble()
  # isopedia_transcript <- DBI::dbGetQuery(conn_isopedia, "SELECT * FROM transcript") %>% as_tibble()
  # 
  # txID_gene <- isopedia_junctions %>%
  #   mutate(start = start + 1) %>%
  #   left_join(y = introverse_jxn,
  #              by = c("chr" = "seqnames","start", "end")) %>%
  #   inner_join(y = transcript_junction,
  #              by = "junctionID") %>%
  #   inner_join(y = isopedia_transcript %>% dplyr::select(txID, transcript_start, transcript_end),
  #              by = "txID") %>%
  #   dplyr::select(gene_name, txID)
  # 
  # txID_gene %>% dplyr::count(gene_name) %>% arrange(gene_name) %>% 
  # filter(gene_name == "ADAR")  
  # 
  # dbWriteTable(conn, "gene", txID_gene %>% distinct(txID, .keep_all = T), overwrite = TRUE)
  # 
  # 
  
}
