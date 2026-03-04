# module load R/4.3.1-icelake
# module load gcc/11
# module load libiconv/1.17/gcc/rlhw7p3w
# module load curl/8.7.1/gcc/d76pzrod  

# options(width=150)

library(DBI)
library(RSQLite)
library(dplyr)

run = FALSE

if (run) {
  # options(width = 150)
  database_path <- "/home/sruiz/POST_DOC/recount3-database-project/database/"
  gtf_version <- 115
  
  args <-
    list(
      SRP058181 = file.path(database_path, "SRP058181_1read_subsampleFALSE", gtf_version, "SRP058181_1read_subsampleFALSE.sqlite"),
      SRP100948 = file.path(database_path, "SRP100948_1read_subsampleFALSE", gtf_version, "SRP100948_1read_subsampleFALSE.sqlite")
    )
  
  
  splicing_events <- map_df(names(args), function(database){
    # database = names(args)[1] # Debugging purposes
    
    # Connect to database
    con <- dbConnect(SQLite(), args[[database]])  
    # Check tables
    dbListTables(con)
    
    ## Get introns
    introns <- DBI::dbGetQuery(con, paste0("SELECT seqnames, start, end, strand, gene_id FROM 'intron'")) %>% as_tibble()
    alternative5ss_3ss <- DBI::dbGetQuery(con, paste0("SELECT novel.seqnames, novel.start, novel.end, novel.strand, intron.gene_id 
                                                    FROM 'novel' as novel
                                                    INNER JOIN 'intron' as intron ON novel.ref_junID = intron.ref_junID")) %>% as_tibble()
    
    novel_combo_unannotated <- DBI::dbGetQuery(con, paste0("SELECT seqnames, start, end, strand, gene_id FROM 'other'")) %>% as_tibble()
    
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
    mutate(chrom2 = chrom) %>%
    relocate(chrom2, .after = chromStart)
  
  
  
  write.table(
    splicing_events_bed,
    file = file.path(database_path, "bed/junctions.bed"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  
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

