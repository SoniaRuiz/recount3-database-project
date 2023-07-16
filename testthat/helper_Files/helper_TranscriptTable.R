print("Loading reference transcriptome...")

if ( !exists("hg38") ) {
  hg38 <- rtracklayer::import(reference_transcriptome_path) 
}

if ( !exists("hg_mane") ) {
  hg_mane <- rtracklayer::import(hg_mane_transcripts_path)
}



## Load the necessary tables from the database
con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
df_transcript <- dplyr::tbl(con, "transcript") %>% dplyr::collect()
df_gene <- dplyr::tbl(con, "gene") %>% dplyr::collect()
DBI::dbDisconnect(con)
