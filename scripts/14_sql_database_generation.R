#' Title
#' SQL database helper function to create the different tables sequentially
#' In case of removing different tables, useful to control which tables are removed 
#' (i.e. only child tables or also master tables) 
#' @param database.sqlite Path to the .sql database file
#' @param recount3.project.IDs List of recount3 projects to analyse
#' @param project.name Name given to the project 
#' @param gtf.path Version of the reference transcriptome to use. In this case it has been used '105' corresponding
#' to Ensembl v105
#' @param remove.all Boolean to reset the database (i.e. whether all tables of the database should be removed)
#'
#' @return
#' @export
#'
#' @examples
SqlDatabaseGeneration <- function(database.sqlite,
                                  recount3.project.IDs,
                                  database.folder,
                                  results.folder,
                                  dependencies.folder,
                                  gtf.path,
                                  max.ent.tool.path,
                                  bedtools.path,
                                  hs.fasta.path,
                                  phastcons.bw.path,
                                  cdts.bw.path, 
                                  mane.gtf.path,
                                  utr.introns.path,
                                  circRNA.path,
                                  miRNA.path,
                                  replace = F,
                                  remove.all = T,
                                  tmp.dir,
                                  discard.minor.introns = F) {
  
  
  print(paste0(Sys.time(), " --> ", database.sqlite, "..."))
  
  con <- dbConnect(RSQLite::SQLite(), database.sqlite)
  tables <- DBI::dbListTables(conn = con)
  
  #logger::log_info("Database tables:",paste0(tables %>% print()))

  if (replace) {
    SqlRemoveTables(database.sqlite, all = remove.all, con)
    tables <- DBI::dbListTables(conn = con)
    tables %>% print()
  }
  
  
  if (!all(c('metadata', 'intron', 'novel', 'gene', 'transcript', 'other') %in% tables)) {
    
    logger::log_info("Creating master tables ...")
    SqlCreateMasterTables(database.sqlite = database.sqlite,
                          gtf.path = gtf.path,
                          database.folder = database.folder,
                          results.folder = results.folder,
                          dependencies.folder = dependencies.folder,
                          recount3.project.IDs = recount3.project.IDs,
                          max.ent.tool.path,
                          bedtools.path,
                          hs.fasta.path,
                          phastcons.bw.path,
                          cdts.bw.path, 
                          mane.gtf.path,
                          utr.introns.path,
                          circRNA.path,
                          miRNA.path,
                          tmp.dir = tmp.dir,
                          discard.minor.introns = discard.minor.introns)
    
    
    tables <- DBI::dbListTables(conn = con)
    tables %>% print()
    
  } else {
    logger::log_info("All master tables exist!")
  }
  
  logger::log_info("Creating child tables ...")
  SqlCreateChildTables(database.sqlite = database.sqlite,
                       database.folder= database.folder,
                       results.folder,
                       recount3.project.IDs)

}