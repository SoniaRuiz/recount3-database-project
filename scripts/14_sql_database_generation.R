#' Title
#' SQL database helper function to create the different tables sequentially
#' In case of removing different tables, useful to control which tables are removed 
#' (i.e. only child tables or also master tables) 
#' @param database.path Path to the .sql database file
#' @param recount3.project.IDs List of recount3 projects to analyse
#' @param project.name Name given to the project 
#' @param gtf.version Version of the reference transcriptome to use. In this case it has been used '105' corresponding
#' to Ensembl v105
#' @param remove.all Boolean to reset the database (i.e. whether all tables of the database should be removed)
#'
#' @return
#' @export
#'
#' @examples
sql_database_generation <- function(database.path,
                                    recount3.project.IDs, 
                                    remove.all = NULL,
                                    database.folder,
                                    results.folder,
                                    gtf.version,
                                    discard.minor.introns = F) {
  
  
  print(paste0(Sys.time(), " --> ", database.path, "..."))
  
  con <- dbConnect(RSQLite::SQLite(), database.path)
  tables <- DBI::dbListTables(conn = con)
  
  #logger::log_info("Database tables:",paste0(tables %>% print()))

  
  if ( !is.null(remove.all) ) {
    sql_remove_tables(database.path, 
                      all = remove.all)
    tables <- DBI::dbListTables(conn = con)
    tables %>% print()
  }
  
  if ( !any(tables == 'metadata') ) {
    sql_create_master_table_metadata(database.path = database.path,
                                     recount3.project.IDs = recount3.project.IDs,
                                     results.folder = results.folder)
    tables <- DBI::dbListTables(conn = con)
    tables %>% print()
  } else {
    logger::log_info("table 'metadata' exists!")
  }
  
  if ( !any(tables %in% c('intron', 'novel', 'gene', 'transcript')) ) {
    logger::log_info("creating master tables ...")
    sql_create_master_tables(database.path = database.path,
                             gtf.version = gtf.version,
                             database.folder = database.folder,
                             results.folder = results.folder,
                             discard.minor.introns = discard.minor.introns)
    tables <- DBI::dbListTables(conn = con)
    tables %>% print()
  } else {
    logger::log_info("master tables exist!")
  }
  
  logger::log_info("creating child tables ...")
  sql_create_child_tables(database.path,
                          recount3.project.IDs,
                          database.folder,
                          results.folder)
  
  
  logger::log_info("creating combo tables ...")
  # sql_create_combo_tables(database.path,
  #                         recount3.project.IDs,
  #                         database.folder,
  #                         results.folder)
}
