#################### Valid data types for each table
### Gene table
gene_expected_types <- c(id = "INTEGER", 
                         gene_id = "TEXT", 
                         gene_name = "TEXT", 
                         n_transcripts = "INTEGER", 
                         gene_width = "INTEGER")
gene_expected_notnull <- c(id = 1, gene_id = 1, gene_name = 0, n_transcripts = 1, gene_width = 1)
gene_expected_pk <- c(id = 1, gene_id = 0, gene_name = 0, n_transcripts = 0, gene_width = 0)


### Transcript table
transcript_expected_types <- c(id = "NUMERIC", 
                               transcript_id = "TEXT", 
                               TSL = "NUMERIC", 
                               MANE = "BOOL", 
                               gene_id = "INTEGER")
transcript_expected_notnull <- c(id = 1, transcript_id = 1, TSL = 1, MANE = 1, gene_id = 1)
transcript_expected_pk <- c(id = 1, transcript_id = 0, TSL = 0, MANE = 0, gene_id = 0)
transcript_expected_foreign <- c(table = "gene", from = "gene_id", to = "id")


### Intron table
intron_expected_types <- c(ref_junID = "NUMERIC", 
                           ref_coordinates = "TEXT", 
                           ref_length = "INTEGER", 
                           ref_ss5score = "DOUBLE", 
                           ref_ss3score = "DOUBLE", 
                           ref_cons5score = "DOUBLE", 
                           ref_cons3score = "DOUBLE", 
                           ref_CDTS5score = "DOUBLE", 
                           ref_CDTS3score = "DOUBLE", 
                           u2_intron = "BOOL",
                           u12_intron = "BOOL", 
                           clinvar = "TEXT", 
                           lncRNA = "INTEGER", 
                           protein_coding = "INTEGER", 
                           misspliced = "BOOL", 
                           transcript_id = "INTEGER")
intron_expected_notnull <- c(ref_junID = 1, ref_coordinates = 1, ref_length = 1, ref_ss5score = 1, ref_ss3score = 1, ref_cons5score = 1, ref_cons3score = 1, ref_CDTS5score = 1, 
                             ref_CDTS3score = 1, u2_intron = 0, u12_intron = 0, clinvar = 1, lncRNA = 1, protein_coding = 1, misspliced = 1, transcript_id = 1)

intron_expected_pk <- c(ref_junID = 1, ref_coordinates = 0, ref_length = 0, ref_ss5score = 0, ref_ss3score = 0, ref_cons5score = 0, ref_cons3score = 0, ref_CDTS5score = 0, 
                        ref_CDTS3score = 0, u2_intron = 0, u12_intron = 0, clinvar = 0, lncRNA = 0, protein_coding = 0, misspliced = 0, transcript_id = 0)

intron_expected_foreign <- c(table = "transcript", from = "transcript_id", to = "id")

### Novel table
novel_expected_types <- c(novel_junID = "NUMERIC", 
                          ref_junID = "NUMERIC", 
                          novel_coordinates = "TEXT", 
                          novel_ss5score = "DOUBLE", 
                          novel_ss3score = "DOUBLE", 
                          novel_length = "INTEGER",
                          novel_type = "TEXT",
                          distance = "INTEGER")
novel_expected_notnull <- c(novel_junID = 1, ref_junID = 1, novel_coordinates = 1, novel_ss5score = 1, novel_ss3score = 1, novel_length = 1, novel_type = 1, distance = 1)
novel_expected_pk <- c(novel_junID = 2, ref_junID = 1, novel_coordinates = 0, novel_ss5score = 0, novel_ss3score = 0, novel_length = 0, novel_type = 0, distance = 0)
novel_expected_foreign <- c(table = "intron", from = "ref_junID", to = "ref_junID")

# ### MANE table
# mane_expected_types <- c(seqnames = "TEXT", 
#                          start = "INTEGER", 
#                          end = "INTEGER", 
#                          strand = "TEXT", 
#                          type = "TEXT", 
#                          gene_name = "TEXT", 
#                          transcript_id = "TEXT", 
#                          transcript_name = "TEXT")
# mane_expected_notnull <- c(seqnames = 0, start = 0, end = 0, strand = 0, type = 0, gene_name = 0, transcript_id = 0, transcript_name = 0)
# mane_expected_pk <- c(seqnames = 0, start = 0, end = 0, strand = 0, type = 0, gene_name = 0, transcript_id = 0, transcript_name = 0)

### Master table
master_expected_types <- c(age = "TEXT", 
                           rin = "REAL", 
                           gender = "REAL", 
                           cluster = "TEXT", 
                           smnabtcht = "TEXT", 
                           sample_id = "TEXT", 
                           smafrze = "TEXT", 
                           avg_read_length = "INTEGER", 
                           mapped_read_count = "INTEGER", 
                           SRA_project = "TEXT")
master_expected_notnull <- c(age = 0, rin = 0, gender = 0, cluster = 0, smnabtcht = 0, sample_id = 0, smafrze = 0, avg_read_length = 0, mapped_read_count = 0, SRA_project = 0)
master_expected_pk <- c(age = 0, rin = 0, gender = 0, cluster = 0, smnabtcht = 0, sample_id = 0, smafrze = 0, avg_read_length = 0, mapped_read_count = 0, SRA_project = 0)

### Cluster_misspliced tables
misspliced_expected_types <- c(ref_junID = "INTEGER", 
                               novel_junID = "INTEGER", 
                               ref_n_individuals = "INTEGER",
                               ref_sum_counts = "INTEGER",
                               ref_type = "TEXT",
                               novel_n_individuals = "INTEGER",
                               novel_sum_counts = "INTEGER",
                               MSR_D = "DOUBLE",
                               MSR_A = "DOUBLE",
                               gene_tpm = "DOUBLE",
                               transcript_id = "INTEGER")
misspliced_expected_notnull <- c(ref_junID = 1, novel_junID = 1, ref_n_individuals = 1, ref_sum_counts = 1, ref_type = 1, novel_n_individuals = 1, novel_sum_counts = 1, MSR_D = 1, MSR_A = 1, 
                                 gene_tpm = 0, 
                                 transcript_id = 1)
misspliced_expected_pk <- c(ref_junID = 0, novel_junID = 0, ref_n_individuals = 0, ref_sum_counts = 0, ref_type = 0, novel_n_individuals = 0, novel_sum_counts = 0, MSR_D = 0, MSR_A = 0, 
                            gene_tpm = 0, 
                            transcript_id = 0)
misspliced_expected_foreign_transcript_id <- c(table = "transcript", from = "transcript_id", to = "id")
misspliced_expected_foreign_ref_junID <- c(table = "novel", from = "ref_junID", to = "ref_junID")
misspliced_expected_foreign_novel_junID <- c(table = "novel", from = "novel_junID", to = "novel_junID")





### Cluster_nevermisspliced tables
never_expected_types <- c(ref_junID = "INTEGER", 
                          ref_n_individuals = "INTEGER",
                          ref_sum_counts = "INTEGER",
                          MSR_D = "DOUBLE",
                          MSR_A = "DOUBLE",
                          ref_type = "TEXT",
                          gene_tpm = "DOUBLE",
                          transcript_id = "INTEGER")
never_expected_notnull <- c(ref_junID = 1, ref_n_individuals = 1, ref_sum_counts = 1, MSR_D = 1, MSR_A = 1, ref_type = 1, 
                            gene_tpm = 0, 
                            transcript_id = 1)
never_expected_pk <- c(ref_junID = 0, ref_n_individuals = 0, ref_sum_counts = 0, MSR_D = 0, MSR_A = 0, ref_type = 0, 
                       gene_tpm = 0, 
                       transcript_id = 0)
never_expected_foreign_gene_id <- c(table = "transcript", from = "transcript_id", to = "id")
never_expected_foreign_ref_junID <- c(table = "intron", from = "ref_junID", to = "ref_junID")

#################### Functions

#' Get the table information
#'
#' @param con DBIConnection object used to direct commands to the database
#'   engine.
#' @param table Name of the table to extract the information from.
#'
#' @return Dataframe containing the information of the input table.
#' @export
getTableInformation <- function(con, table) {
  table_info <- DBI::dbGetQuery(con, paste0("PRAGMA table_info('", table, "')")) %>% tibble::as_tibble()
  return(table_info)
}

#' Get the table foreign keys
#'
#' @param con DBIConnection object used to direct commands to the database
#'   engine.
#' @param table Name of the table to extract the foreign keys from.
#'
#' @return Dataframe containing the foreign keys of the input table.
#' @export
getForeignKeyInformation <- function(con, table){
  keys_info <- DBI::dbGetQuery(con, paste0("PRAGMA foreign_key_list('", table, "')")) %>% tibble::as_tibble()
  return(keys_info)
}