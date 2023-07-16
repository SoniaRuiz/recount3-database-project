main_path <- normalizePath(path = "./")
source(paste0(main_path, "/testthat/helper_Files/helper_Global.R"))
skip_if(!test_DatabaseTables, "Database tables tests not executed. Variable test_DatabaseTables set to FALSE in global options.")
source(paste0(main_path, "/testthat/helper_Files/helper_DatabaseTables.R"))

context("\tTest existence of Main tables")
test_that("Test existence of Main tables", {
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  
  expect_equal(getMainTables(con), expected_main_tables)
  DBI::dbDisconnect(con)
})


context("\tTest pairing of tissue tables")
test_that("Test pairing of tissue tables", {
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  tissue_tables <- getTissueTables(con)
  
  misspliced_tables <- tissue_tables[endsWith(tissue_tables, "_misspliced")]
  nevermisspliced_tables <- tissue_tables[endsWith(tissue_tables, "_nevermisspliced")]
  
  expect_equal(length(misspliced_tables), length(nevermisspliced_tables))
  expect_equal(length(misspliced_tables) + length(nevermisspliced_tables), length(tissue_tables))
  
  misspliced_tissues <- stringr::str_remove_all(misspliced_tables, "_misspliced")
  nevermisspliced_tissues <- stringr::str_remove_all(nevermisspliced_tables, "_nevermisspliced")
  
  expect_equal(misspliced_tissues, nevermisspliced_tissues)
  
  DBI::dbDisconnect(con)
})

clearAllVariables()
