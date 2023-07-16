main_path <- normalizePath(path = "./")
source(paste0(main_path, "/testthat/helper_Files/helper_Global.R"))
skip_if(!test_DataTypes, "Data types and foreign keys tests not executed. Variable test_DataTypes set to FALSE in global options.")
source(paste0(main_path, "/testthat/helper_Files/helper_DatabaseTables.R"))
source(paste0(main_path, "/testthat/helper_Files/helper_DataTypes.R"))

context("\tTest that gene table data types and foreign keys are correctly set up")
test_that("Test that gene table data types and foreign keys are correctly set up", {
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  
  gene_table <- getTableInformation(con, "gene")
  gene_table_types <- gene_table %>% select(name, type) %>% tibble::deframe()
  gene_table_notnull <- gene_table %>% select(name, notnull) %>% tibble::deframe()
  gene_table_pk <- gene_table %>% select(name, pk) %>% tibble::deframe()
  
  expect_equal(gene_table_types, gene_expected_types)
  expect_equal(gene_table_notnull, gene_expected_notnull)
  expect_equal(gene_table_pk, gene_expected_pk)
  DBI::dbDisconnect(con)
})


context("\tTest that 'transcript' table data types and foreign keys are correctly set up")
test_that("Test that 'transcript' table data types and foreign keys are correctly set up", {
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  
  transcript_table <- getTableInformation(con, "transcript")
  transcript_table_types <- transcript_table %>% select(name, type) %>% tibble::deframe()
  transcript_table_notnull <- transcript_table %>% select(name, notnull) %>% tibble::deframe()
  transcript_table_pk <- transcript_table %>% select(name, pk) %>% tibble::deframe()
  
  expect_equal(transcript_expected_types, transcript_table_types)
  expect_equal(transcript_expected_notnull, transcript_table_notnull)
  expect_equal(transcript_expected_pk, transcript_table_pk)
  DBI::dbDisconnect(con)
  
})

context("\tTest that intron table data types and foreign keys are correctly set up")
test_that("Test that intron table data types and foreign keys are correctly set up", {
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  
  intron_table <- getTableInformation(con, "intron")
  intron_table_types <- intron_table %>% select(name, type) %>% tibble::deframe()
  intron_table_notnull <- intron_table %>% select(name, notnull) %>% tibble::deframe()
  intron_table_pk <- intron_table %>% select(name, pk) %>% tibble::deframe()
  
  expect_equal(intron_expected_types, intron_table_types)
  expect_equal(intron_expected_notnull, intron_table_notnull)
  expect_equal(intron_expected_pk, intron_table_pk)
  
  intron_foreign_keys <- getForeignKeyInformation(con, "intron") %>% select(table, from, to) %>% .[1, ] %>% unlist()
  expect_equal(intron_expected_foreign, intron_foreign_keys)
  
  DBI::dbDisconnect(con)
})

context("\tTest that novel table data types and foreign keys are correctly set up")
test_that("Test that novel table data types and foreign keys are correctly set up", {
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  
  novel_table <- getTableInformation(con, "novel")
  novel_table_types <- novel_table %>% select(name, type) %>% tibble::deframe()
  novel_table_notnull <- novel_table %>% select(name, notnull) %>% tibble::deframe()
  novel_table_pk <- novel_table %>% select(name, pk) %>% tibble::deframe()
  
  expect_equal(novel_expected_types, novel_table_types)
  expect_equal(novel_expected_notnull, novel_table_notnull)
  expect_equal(novel_expected_pk, novel_table_pk)
  
  novel_foreign_keys <- getForeignKeyInformation(con, "novel") %>% select(table, from, to) %>% .[1, ] %>% unlist()
  expect_equal(novel_expected_foreign, novel_foreign_keys)
  DBI::dbDisconnect(con)
})



context("\tTest that master table data types and foreign keys are correctly set up")
test_that("Test that master table data types and foreign keys are correctly set up", {
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  
  master_table <- getTableInformation(con, "master")
  master_table_types <- master_table %>% select(name, type) %>% tibble::deframe()
  master_table_notnull <- master_table %>% select(name, notnull) %>% tibble::deframe()
  master_table_pk <- master_table %>% select(name, pk) %>% tibble::deframe()
  
  expect_equal(master_expected_types, master_table_types)
  expect_equal(master_expected_notnull, master_table_notnull)
  expect_equal(master_expected_pk, master_table_pk)
  DBI::dbDisconnect(con)
})

context("\tTest that cluster mis-spliced tables data types and foreign keys are correctly set up")
test_that("Test that cluster mis-spliced tables data types and foreign keys are correctly set up", {
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)

  tissue_tables <- getTissueTables(con)
  misspliced_tables <- tissue_tables[endsWith(tissue_tables, "_misspliced")]
  
  for(table in misspliced_tables){
    # table <- misspliced_tables[1]
    misspliced_table_foreign_keys <- getForeignKeyInformation(con, table) %>% select(table, from, to) %>% arrange(table, from, to)
    misspliced_expected_foreign <- rbind(misspliced_expected_foreign_transcript_id, 
                                         misspliced_expected_foreign_novel_junID, 
                                         misspliced_expected_foreign_ref_junID) %>% tibble::as_tibble() %>% arrange(table, from, to)
    
    misspliced_table_info <- getTableInformation(con, table)
    misspliced_table_types <- misspliced_table_info %>% select(name, type) %>% tibble::deframe()
    misspliced_table_notnull <- misspliced_table_info %>% select(name, notnull) %>% tibble::deframe()
    misspliced_table_pk <- misspliced_table_info %>% select(name, pk) %>% tibble::deframe()
    
    expect_equal(misspliced_expected_types, misspliced_table_types)
    expect_equal(misspliced_expected_notnull, misspliced_table_notnull)
    expect_equal(misspliced_expected_pk, misspliced_table_pk)
    expect_equal(misspliced_expected_foreign, misspliced_table_foreign_keys)
  }
  
  DBI::dbDisconnect(con)
})

context("\tTest that cluster never mis-spliced tables data types and foreign keys are correctly set up")
test_that("Test that cluster never mis-spliced tables data types and foreign keys are correctly set up", {
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)

  tissue_tables <- getTissueTables(con)
  nevermisspliced_tables <- tissue_tables[endsWith(tissue_tables, "_nevermisspliced")]
  
  for(table in nevermisspliced_tables){
    never_table_foreign_keys <- getForeignKeyInformation(con, table) %>% select(table, from, to) %>% arrange(table, from, to)
    never_expected_foreign <- rbind(never_expected_foreign_gene_id, never_expected_foreign_ref_junID) %>% tibble::as_tibble() %>% arrange(table, from, to)
    
    never_table_info <- getTableInformation(con, table)
    never_table_types <- never_table_info %>% select(name, type) %>% tibble::deframe()
    never_table_notnull <- never_table_info %>% select(name, notnull) %>% tibble::deframe()
    never_table_pk <- never_table_info %>% select(name, pk) %>% tibble::deframe()
    
    expect_equal(never_expected_types, never_table_types)
    expect_equal(never_expected_notnull, never_table_notnull)
    expect_equal(never_expected_pk, never_table_pk)
    expect_equal(never_expected_foreign, never_table_foreign_keys)
  }
  
  DBI::dbDisconnect(con)
})

clearAllVariables()
