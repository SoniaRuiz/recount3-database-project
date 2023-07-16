#' Get the main tables from the database
#'
#' @param con DBIConnection object used to direct commands to the database
#'   engine.
#'
#' @return List of main tables in the database.
#' @export
getMainTables <- function(con) {
  tables <- DBI::dbListTables(con)
  main_tables <- tables[!grepl("misspliced", tables)] %>% sort()
  return(main_tables)
}
