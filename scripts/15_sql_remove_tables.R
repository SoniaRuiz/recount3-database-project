#' Title
#' Removes either all or only the child tables from the database.
#' @param database.path Local path to the .sqlite database
#' @param all Whether to remove all tables or only the child tables.
#'
#' @return
#' @export
#'
#' @examples
SqlRemoveTables <- function(database.path, all, con = NULL) {
  
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=0")
  
  tables <- dbListTables(con)
  tables %>% print()
  
  for (table in tables) {
    
    if (!all) {
      
      if ( !(table %in% c("gene", "transcript", "intron", "metadata", "novel")) ) {
        dbRemoveTable(conn = con, table)
        print(paste0("Table: '", table, "' has been removed!"))
      }
      
    } else {
      
      dbRemoveTable(conn = con, table)
      print(paste0("Table: '", table, "' has been removed!"))
      
    }
  }
}
