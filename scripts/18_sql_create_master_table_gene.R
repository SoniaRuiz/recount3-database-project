#' Title
#' Creates the master table 'gene'
#' @param database.path Local path to the .sqlite databse
#' @param hg38 Human genome. Needed to extract info about the transcripts to be stored.
#' @param gene.ids List of Ensembl gene IDs
#'
#' @return
#' @export
#'
#' @examples
sql_create_master_table_gene <- function(database.path,
                                         #protein_biotype,
                                         hg38,
                                         gene.ids) {
  
  
  logger::log_info(" --> Creating 'gene' master table...")
  con <- dbConnect(RSQLite::SQLite(), database.path)
  
  
  hg38_transcripts_gene <- hg38 %>%
    as_tibble() %>%
    mutate(gene_id = str_sub(gene_id, start = 1, end = 15)) %>%
    dplyr::filter(gene_id %in% gene.ids$gene_id) %>%
    dplyr::count(gene_id, type) %>%
    dplyr::filter(type == "transcript") %>%
    unnest(gene_id) %>%
    dplyr::select(-type) %>%
    dplyr::rename(n_transcripts = n)
  
  
  hg38_genes <- hg38 %>%
    as_tibble() %>%
    mutate(gene_id = str_sub(gene_id, start = 1,end = 15)) %>%
    dplyr::filter(gene_id %in% gene.ids$gene_id) %>%
    dplyr::filter(type == "gene") %>%
    dplyr::select(gene_id, gene_name, gene_width = width, gene_biotype)
  
  
  hg38_tidy <- hg38_transcripts_gene %>%
    inner_join(hg38_genes, by = "gene_id" ) %>% 
    as_tibble()
  
  
  
  ## CREATE GENE_NAME TABLE
  sql_statement <- paste0("CREATE TABLE IF NOT EXISTS 'gene'", 
                          "(id INTEGER PRIMARY KEY NOT NULL,
                          gene_id TEXT NOT NULL,
                          gene_name TEXT,
                          gene_width INTEGER NOT NULL, 
                          n_transcripts INTEGER NOT NULL,
                          gene_biotype TEXT NOT NULL)")
  res <- DBI::dbSendQuery(con, sql_statement)
  DBI::dbClearResult(res)
  
  ## POPULATE GENE TABLE
  hg38_tidy_final <- hg38_tidy %>% 
    as_tibble() %>%
    tibble::rowid_to_column("id")
  
  DBI::dbAppendTable(conn = con,
                     name = "gene", 
                     value = hg38_tidy_final)
  
  
  logger::log_info(" --> gene table created!")
  logger::log_info(" --> ", hg38_tidy_final %>% nrow(), " genes stored!")
  
  
  
  ## CREATE INDEXES ---------------------------------------------------
  
  query <- paste0("CREATE UNIQUE INDEX 'index_gene_id' ON 'gene'(id)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  query <- paste0("CREATE UNIQUE INDEX 'index_gene_ensembl_id' ON 'gene'(gene_id,gene_name)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  DBI::dbDisconnect(con)
}