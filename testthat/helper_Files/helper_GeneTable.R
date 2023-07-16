print("Loading reference transcriptome...")
hg38 <- rtracklayer::import(reference_transcriptome_path)

## Load the necessary tables from the database
con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
df_gene <- dplyr::tbl(con, "gene") %>% dplyr::collect()
DBI::dbDisconnect(con)
