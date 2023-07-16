## Load the necessary tables from the database
con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
df_novel <- dplyr::tbl(con, "novel") %>% dplyr::collect()
df_intron <- dplyr::tbl(con, "intron") %>% dplyr::collect()
df_gene <- dplyr::tbl(con, "gene") %>% dplyr::collect()
df_transcript <- dplyr::tbl(con, "transcript") %>% dplyr::collect()
if ( !exists("hg_mane") ) {
  hg_mane <- rtracklayer::import(hg_mane_transcripts_path)
}
