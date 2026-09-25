# module load R/4.3.1-icelake
# module load gcc/11
# module load libiconv/1.17/gcc/rlhw7p3w
# module load curl/8.7.1/gcc/d76pzrod  

# options(width=150)

library(DBI)
library(RSQLite)
library(dplyr)
library(tidyverse)


database_path <- "/rds/user/sg2173/hpc-work/moved_from_home_dir/PROJECTS/srRNAseq/introverse_v2/database/"
out_dir <- file.path(database_path, "bed")

dir.create(path = out_dir, recursive = T, showWarnings = F)
gtf_version <- 115

run = FALSE

if (run) {
  
  # options(width = 150)
   
  args <-
    list(
      "SRP058181" = file.path(database_path, "SRP058181_1read_subsampleFALSE", gtf_version, "SRP058181_1read_subsampleFALSE.sqlite"),
      "SRP100948" = file.path(database_path, "SRP100948_1read_subsampleFALSE", gtf_version, "SRP100948_1read_subsampleFALSE.sqlite"),
      "TCGA" = file.path(database_path, "TCGA_1read_subsampleFALSE", gtf_version, "TCGA_1read_subsampleFALSE.sqlite")
    )
  
  
  splicing_events <- map_df(names(args), function(database){
    # database = names(args)[1] # Debugging purposes
    
    # Connect to database
    con <- dbConnect(SQLite(), args[[database]])  
    # Check tables
    dbListTables(con)
    
    message(database)
    
    ## Get introns
    introns <- DBI::dbGetQuery(con, paste0("SELECT seqnames, start, end, strand, gene_id FROM 'intron'")) %>% as_tibble()
    
    # Get novel donor and novel acceptor junctions
    alternative5ss_3ss <- DBI::dbGetQuery(con, paste0("SELECT novel.seqnames, novel.start, novel.end, novel.strand, intron.gene_id 
                                                    FROM 'novel' as novel
                                                    INNER JOIN 'intron' as intron ON novel.ref_junID = intron.ref_junID")) %>% as_tibble()
    
    # Get novel combo and completely unannotated junctions
    novel_combo_unannotated <- DBI::dbGetQuery(con, paste0("SELECT seqnames, start, end, strand, gene_id FROM 'other'")) %>% as_tibble()
    
    # Get gene names
    gene_names <- DBI::dbGetQuery(con, paste0("SELECT gene_id, gene_name FROM 'gene'")) %>% as_tibble()
    
    rbind(introns, alternative5ss_3ss, novel_combo_unannotated) %>%
      left_join(y = gene_names,
                by = "gene_id")
  })
  
  
  
  splicing_events_bed <- splicing_events %>% 
    distinct(seqnames, start, end, strand, .keep_all = T) %>% 
    dplyr::select(seqnames, start, end, strand, gene_name) %>%
    mutate(
      chromStart = start,   # BED uses 0-based start
      chromEnd = end,           # end stays as-is
      name = paste0("junction_", row_number()),  # unique name per junction
      score = 0
    ) %>%
    dplyr::select(chrom = seqnames, chromStart,chromEnd, name, score, strand) %>%
    arrange(chrom, chromStart, chromEnd, strand) |>
    mutate(chrom2 = chrom) %>%
    relocate(chrom2, .after = chromStart)
  
  splicing_events_bed |> distinct(chrom, chromStart, chrom2, chromEnd, .keep_all = T)
  splicing_events_bed |> filter(name == "junction_13647393")
  
  write.table(
    splicing_events_bed,
    file = file.path(out_dir, "junctions.bed"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  
  # ----------------------------------------------------------------------------
  # TEST SUBSET
  # ----------------------------------------------------------------------------
  
  # Assumes `splicing_events_bed` (already sorted by chrom, chromStart, chromEnd, strand)
  # is already in memory from the previous step.
  
  # 27,049,404
  # splicing_events_bed <- read.table(file.path(out_dir, "junctions.bed"))


  chunk_out_dir <- file.path(out_dir, "chunks")
  dir.create(chunk_out_dir, showWarnings = FALSE, recursive = TRUE)


  chunk_size <- 1000000

  
  
  n_total <- nrow(splicing_events_bed)
  n_chunks <- ceiling(n_total / chunk_size)
  
  for (i in seq_len(n_chunks)) {

    # i = 1

    start_row <- (i - 1) * chunk_size + 1
    end_row   <- min(i * chunk_size, n_total)
  
    chunk_bed <- splicing_events_bed %>%
      slice(start_row:end_row) %>%
      mutate(
        chromStart = as.integer(chromStart - 1),
        chromEnd = as.integer(chromEnd)
      )
  
    write.table(
      chunk_bed,
      file      = file.path(chunk_out_dir, sprintf("chunk_%03d.bed", i - 1)),
      quote     = FALSE,
      sep       = "\t",
      row.names = FALSE,
      col.names = FALSE
    )
  }
  
  cat("Wrote", n_chunks, "chunks to", out_dir, "\n") 
  # test_n <- 1000000  # <-- adjust subset size here if needed
  
  # test_bed <- read.table(file.path(out_dir, "junctions.bed")) %>%
  #   slice(1:test_n) 

  
  # test_bed_tidy <- test_bed |> 
  #   as_tibble() |> 
  #   dplyr::rename(chrom = V1, chromStart =V2, chrom2 = V3, chromEnd = V4, name = V5, score = V6, strand = V7)  %>%
  #   mutate(
  #     chromStart = as.integer(chromStart - 1),
  #     chromEnd = as.integer(chromEnd)
  #   )
  
  # write.table(
  #   test_bed_tidy,
  #   file      = file.path(out_dir, "test_junctions_1m.bed"),
  #   quote     = FALSE,
  #   sep       = "\t",
  #   row.names = FALSE,
  #   col.names = FALSE
  # )
   
  # variants <- list(
  #   pos1_minus1 = test_bed %>% mutate(chromStart = chromStart - 1),
  #   pos2_minus1 = test_bed %>% mutate(chromEnd = chromEnd - 1),
  #   both_minus1 = test_bed %>% mutate(chromStart = chromStart - 1, chromEnd = chromEnd - 1),
  #   pos1_plus1  = test_bed %>% mutate(chromStart = chromStart + 1),
  #   pos2_plus1  = test_bed %>% mutate(chromEnd = chromEnd + 1),
  #   both_plus1  = test_bed %>% mutate(chromStart = chromStart + 1, chromEnd = chromEnd + 1)
  # )
  # 
  # out_dir <- file.path(out_dir, "offset_test")
  # dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  # 
  # for (vname in names(variants)) {
  #   write.table(
  #     variants[[vname]],
  #     file      = file.path(out_dir, paste0("test_junctions_10k_", vname, ".bed")),
  #     quote     = FALSE,
  #     sep       = "\t",
  #     row.names = FALSE,
  #     col.names = FALSE
  #   )
  # }
  # 
  # ------------------------------------------------------------------------------
  # QUERY ISOPEDIA INDEX
  # 
  # Run using screen
  # ------------------------------------------------------------------------------
  
  # cd /home/sruiz/POST_DOC/isopedia/
  # ./isopedia splice -i /home/sruiz/POST_DOC/isopedia_index/ \
  #     -S /home/sruiz/POST_DOC/recount3-database-project/database/bed/junctions.bed \
  #     -o /home/sruiz/POST_DOC/recount3-database-project/database/bed/isopedia_introverse_splice_junctions.tsv.gz \
  #     -f 1
}

