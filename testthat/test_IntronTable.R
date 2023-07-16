main_path <- normalizePath(path = "./")
source(paste0(main_path, "/testthat/helper_Files/helper_Global.R"))
skip_if(!test_IntronTable, "Intron table tests not executed. Variable test_IntronTable set to FALSE in global options.")
source(paste0(main_path, "/testthat/helper_Files/helper_IntronTable.R"))

context("\tTest consistent junction IDs with other tables")
test_that("Test consistent junction IDs with other tables", {
  intron_ref_junID <- df_intron %>% pull(ref_junID)
  intron_ref_junID_misspliced <- df_intron %>% filter(misspliced == T) %>% pull(ref_junID) %>% sort()
  intron_ref_junID_nevermisspliced <- df_intron %>% filter(misspliced == F) %>% pull(ref_junID) %>% sort()
  novel_ref_junID <- df_novel %>% pull(ref_junID) %>% unique() %>% sort()
  
  ## All introns categorized as `misspliced` are exactly the same introns that appear as `ref_junID` in `novel` table.
  expect_equal(intron_ref_junID_misspliced, novel_ref_junID)

  ## No intron categorized as `never mis-spliced` is found as `ref_junID` in `novel` table. 
  expect_equal(setdiff(intron_ref_junID_nevermisspliced, novel_ref_junID), intron_ref_junID_nevermisspliced)
  
  ## No duplicates found in `ref_junID`
  expect_equal(intron_ref_junID %>% unique() %>% length, intron_ref_junID %>% length())
})

context("\tTest reference between intron and transcript table")
test_that("Test reference between intron and transcript table", {
  transcript_id <- df_intron %>% pull(transcript_id)
  
  ## All `transcript_id`` are properly referenced to a row in the gene table
  expect_true(all(transcript_id %in% df_transcript$id))
})

context("\tTest that different fields have consistent values")
test_that("Test that different fields have consistent values", {
  u2_intron <- df_intron %>% pull(u2_intron)
  u12_intron <- df_intron %>% pull(u12_intron)
  
  clinvar <- df_intron %>% pull(clinvar) %>% unique()
  clinvar_split <- str_split_fixed(clinvar, "-", n = 2)
  clinvar_types <- clinvar_split[, 1] %>% unique()
  clinvar_values <- clinvar_split[, 2]
  
  #mane <- df_intron %>% pull(MANE)
  #tsl <- df_intron %>% pull(TSL) %>% unique
  lncRNA <- df_intron %>% pull(lncRNA)
  protein_coding <- df_intron %>% pull(protein_coding)
  misspliced <- df_intron %>% pull(misspliced)
  
  ## No intron is classified as u2 and u12 at the same time
  expect_false(any(u2_intron & u12_intron))
  novel_type <- df_novel %>% pull(novel_type) %>% unique() %>% sort()
  
  ## All clinvar structure is either "-" or "[acceptor/donor]-[number]bp"
  expect_true(all(clinvar_types %in% valid_intron_clinvar))
  expect_match(clinvar_values[0:3], regexp = "^$|^[0-9]+bp$")
  
  ## All MANE values are either 0 or 1 (boolean)
  #expect_true(all((mane %>% unique) %in% c(0, 1)))
  
  ## All TSL values are between the valid values (1, 2, 3, 4, 5 and 10)
  #expect_true(all(tsl %in% valid_intron_TSL))
  
  ## All lncRNA values, protein_coding values and their sum are between 0 and 100%
  expect_true(all(lncRNA %>% dplyr::between(0, 100)))
  expect_true(all(protein_coding %>% dplyr::between(0, 100)))
  expect_true(all((protein_coding + lncRNA) %>% dplyr::between(0, 100)))
  
  ## All misspliced values are either 0 or 1 (boolean)
  expect_true(all((misspliced %>% unique) %in% c(0, 1)))
})

context("\tTest ref_coordinates field")
test_that("Test ref_coordinates field", {
  ## We use a custom function to extract the information from the coordinates
  GRdata <- getGRdata(df_intron$ref_junID, df_intron$ref_coordinates, fix_data_types = F, remove_chr = F)
  ref_coordinates <- df_intron$ref_coordinates
  ref_length <- df_intron %>% pull(ref_length)
  novel_coordinates <- df_novel$novel_coordinates
  
  ref_seqnames <- GRdata$seqnames
  ref_strand <- GRdata$strand
  ref_start = GRdata$start
  ref_end = GRdata$end
  
  ## No duplicates found in `ref_coordinates`
  expect_equal(ref_coordinates %>% length(), ref_coordinates %>% unique() %>% length())
  
  ## No `ref_coordinates` value is found in `novel_coordinates`
  expect_length(intersect(ref_coordinates, novel_coordinates), n = 0)
  
  ## No `*` strand found in `ref_coordinates`
  expect_equal(valid_intron_strands, ref_strand %>% unique() %>% sort()) 
  
  ## All `seqnames` in `ref_coordinates` are within the 1-22, X and Y chromosomes
  expect_true(all((ref_seqnames %>% unique()) %in% valid_intron_seqnames))
  
  ## All `start` and `end` positions are integers (or numeric)
  expect_error(as.numeric(ref_start), NA)
  expect_error(as.numeric(ref_end), NA)
  
  ## Intron length is consistent with the coordinates. We use GRanges to
  ## calculate the width (or length) of the intron.
  ref_start <- ref_start %>% as.numeric()
  ref_end <- ref_end %>% as.numeric()
  
  estimated_width <- GRdata %>%
    GenomicRanges::GRanges() %>%
    tibble::as_tibble() %>%
    dplyr::pull(width)
  
  expect_equal(estimated_width, ref_length)
})

context("\tTest that no NAs are found in where they are not allowed")
test_that("Test that no NAs are found in where they are not allowed", {
  ## No single value should be NA or equal to "NA"
  expect_false(any(df_intron %>% is.na()))
  expect_false(any(df_intron == "NA"))
})

context("\tTest that MaxEntScan scores are correctly calculated")
test_that("Test that MaxEntScan scores are correctly calculated", {
  skip_if(!test_MaxEntScan, "No MaxEntScan test. Variable test_MaxEntScan set to FALSE in global options.")
  
  GRdata <- getGRdata(df_intron$ref_junID, df_intron$ref_coordinates)
  GRdata_MaxEntScan <- generateMaxEntScore(GRdata, 
                                           bedtools_path = bedtools_path, 
                                           fasta_path = fasta_path, 
                                           fordownload_path = fordownload_path)
  
  ## Both tables have similar number of rows
  expect_equal(nrow(df_intron), nrow(GRdata_MaxEntScan))
  
  ## They share the same junIDs
  expect_equal(df_intron$ref_junID, GRdata_MaxEntScan$junID)
  
  ## All `MaxEntScan` scores calculated are consistent with the ones found in the `intron` table
  expect_equal(df_intron$ref_ss5score, GRdata_MaxEntScan$ss5score)
  expect_equal(df_intron$ref_ss3score, GRdata_MaxEntScan$ss3score)
})

context("\tTest that the Clinvar data is correctly calculated")
test_that("Test that the Clinvar data is correctly calculated", {
  skip_if(!test_clinvar, "No Clinvar test. Variable test_clinvar set to FALSE in global options.")
  ## Load the intron data
  GRdata <- getGRdata(df_intron$ref_junID, df_intron$ref_coordinates)
  GRdata_clinvar <- GRdata %>%
    mutate(clinvar_type = "-") %>%
    GRanges() 
  
  ## Load the clinvar data
  clinvar_tidy <- readRDS(clinvar_path) %>%
    GenomicRanges::GRanges() %>%
    diffloop::rmchr()
  
  ## Look for overlaps
  overlaps <- GenomicRanges::findOverlaps(query = clinvar_tidy,
                                          subject = GRdata_clinvar,
                                          ignore.strand = F,
                                          type = "any")
  subject_hits <- S4Vectors::subjectHits(overlaps); query_hits <- S4Vectors::queryHits(overlaps)
  
  ## Fill the metadata
  elementMetadata(GRdata_clinvar)[subject_hits, "clinvar_type"] <- ifelse(elementMetadata(clinvar_tidy)[query_hits, "splice_site"] == "none",
                                                                  "-", 
                                                                  paste0(elementMetadata(clinvar_tidy)[query_hits, "splice_site"], 
                                                                         "-",
                                                                         elementMetadata(clinvar_tidy)[query_hits, "distance"], 
                                                                         "bp"))
  
  GRdata_clinvar <- GRdata_clinvar %>% tibble::as_tibble()
  
  ## Both tables have similar number of rows
  expect_equal(nrow(df_intron), nrow(GRdata_clinvar))
  
  ## They share the same junIDs
  expect_equal(df_intron$ref_junID, GRdata_clinvar$junID)
  
  ## All `clinvar` values calculated are consistent with the ones found in the `intron` table
  expect_equal(df_intron$clinvar, GRdata_clinvar$clinvar_type)
})

context("\tTest that the CDTS and conservation scores are correctly obtained")
test_that("Test that the CDTS and conservation scores are correctly obtained", {
  skip_if(!test_conservation_CDTS, "No conservation/CDTS test. Variable test_conservation_CDTS set to FALSE in global options.")
  
  ## Load the reference data
  if(!exists("CNC_CDTS_CONS_gr")) load(CNC_CDTS_CONS_gr_path)
  CNC_CDTS_CONS_gr <- CNC_CDTS_CONS_gr %>% diffloop::rmchr()
  
  ## Load the intron data
  GRdata <- getGRdata(df_intron$ref_junID, df_intron$ref_coordinates)
  GRdata_CDTS <- GRdata %>%
    mutate(CDTS_5ss_mean = 0.0,
           CDTS_3ss_mean = 0.0,
           phastCons20way_5ss_mean = 0.0,
           phastCons20way_3ss_mean = 0.0) %>%
    GRanges() 
  
  ## Scores from the 5'ss
  overlaps <- GenomicRanges::findOverlaps(query = CNC_CDTS_CONS_gr,
                                          subject = GRdata_CDTS %>%
                                            `start<-`(value = start(GRdata_CDTS) - 5) %>%
                                            `end<-`(value = start(GRdata_CDTS) + 35) %>%
                                            `elementMetadata<-`(value = NULL),
                                          ignore.strand = F,
                                          type = "any")
  subject_hits <- S4Vectors::subjectHits(overlaps); query_hits <- S4Vectors::queryHits(overlaps)
  
  ## Get the values for the subject hits
  overlaps_tidy <- overlaps %>%
    tibble::as_tibble() %>%
    mutate(CDTS = CNC_CDTS_CONS_gr[query_hits, ]$CDTS,
           mean_phastCons20way = CNC_CDTS_CONS_gr[query_hits, ]$mean_phastCons20way) %>%
    group_by(subjectHits) %>%
    mutate(CDTS_mean = mean(CDTS, na.rm = T),
              mean_phastCons20way_mean = mean(mean_phastCons20way, na.rm = T))
  
  ## Add the values
  GRdata_CDTS[subject_hits, ]$CDTS_5ss_mean <- overlaps_tidy$CDTS_mean
  GRdata_CDTS[subject_hits, ]$phastCons20way_5ss_mean <- overlaps_tidy$mean_phastCons20way_mean
  
  ## Scores from the 3'ss
  overlaps <- GenomicRanges::findOverlaps(query = CNC_CDTS_CONS_gr %>% diffloop::rmchr(),
                                          subject = GRdata_CDTS %>%
                                            `start<-`(value = end(GRdata_CDTS) - 35) %>%
                                            `end<-`(value = end(GRdata_CDTS) + 5),
                                          ignore.strand = F,
                                          type = "any")
  subject_hits <- S4Vectors::subjectHits(overlaps); query_hits <- S4Vectors::queryHits(overlaps)
  
  ## Get the values for the subject hits
  overlaps_tidy <- overlaps %>%
    tibble::as_tibble() %>%
    mutate(CDTS = CNC_CDTS_CONS_gr[query_hits, ]$CDTS,
           mean_phastCons20way = CNC_CDTS_CONS_gr[query_hits, ]$mean_phastCons20way) %>%
    group_by(subjectHits) %>%
    mutate(CDTS_mean = mean(CDTS),
           mean_phastCons20way_mean = mean(mean_phastCons20way))
  
  ## Add the values
  GRdata_CDTS[subject_hits, ]$CDTS_3ss_mean <- overlaps_tidy$CDTS_mean
  GRdata_CDTS[subject_hits, ]$phastCons20way_3ss_mean <- overlaps_tidy$mean_phastCons20way_mean
  
  ## Convert the data back to tibble
  GRdata_CDTS <- GRdata_CDTS %>% tibble::as_tibble()
  GRdata_CDTS <- GRdata_CDTS %>% 
    mutate(phastCons20way_5ss_mean = phastCons20way_5ss_mean %>% replace(is.na(.), 0),
           phastCons20way_3ss_mean = phastCons20way_3ss_mean %>% replace(is.na(.), 0))
  
  ## Both tables have similar number of rows
  expect_equal(nrow(df_intron), nrow(GRdata_CDTS))
  
  ## They share the same junIDs
  expect_equal(df_intron$ref_junID, GRdata_CDTS$junID)
  
  ## All `CDTS` values calculated are consistent with the ones found in the `intron` table
  expect_equal(df_intron$ref_CDTS5score, GRdata_CDTS$CDTS_5ss_mean)
  expect_equal(df_intron$ref_CDTS3score, GRdata_CDTS$CDTS_3ss_mean)
  
  ## All `conservation` values calculated are consistent with the ones found in the `intron` table
  expect_equal(df_intron$ref_cons5score, GRdata_CDTS$phastCons20way_5ss_mean)
  expect_equal(df_intron$ref_cons3score, GRdata_CDTS$phastCons20way_3ss_mean)
})

# context("\tTest that the TSL values are correctly obtained")
# test_that("Test that the TSL values are correctly obtained", {
#   skip_if(!test_TSL, "No TSL value test. Variable test_TSL set to FALSE in global options.")
#   
#   ## Load the reference transcriptome
#   hg38_transcripts <- rtracklayer::import(reference_transcriptome_path) %>% 
#     as_tibble() %>%
#     filter(type == "transcript") %>%
#     select(transcript_id, transcript_support_level) %>%
#     mutate(transcript_support_level = str_sub(transcript_support_level, 1, 1))
#   
#   ## Load the introns with transcripts IDs
#   df_all_introns_introverse_tidy = readRDS(df_all_introns_introverse_tidy_path)
#   
#   ## Obtain the TSL and select the minimum value for each intron
#   all_intron_tsl <- df_all_introns_introverse_tidy %>%
#     select(ref_junID, tx_id_junction) %>%
#     distinct() %>%
#     tidyr::unnest(tx_id_junction) %>% 
#     left_join(hg38_transcripts,
#               by = c("tx_id_junction" = "transcript_id")) %>%
#     mutate(TSL = transcript_support_level %>% replace(is.na(.) | . == "N", 10) %>% as.integer()) %>%
#     group_by(ref_junID) %>%
#     summarise(TSL = min(TSL))
#   
#   ## Add the calculated values to the reference values
#   df_combined <- df_intron %>% left_join(all_intron_tsl %>% select(ref_junID, TSL), 
#                                          by = c("ref_coordinates" = "ref_junID"))
#   
#   ## Both tables have similar number of rows
#   expect_equal(nrow(df_intron), nrow(all_intron_tsl))
#   
#   ## They share the same junIDs
#   expect_true(length(intersect(df_intron$ref_coordinates, all_intron_tsl$ref_junID)) == length(df_intron$ref_coordinates))
#   
#   ## All `TSL` values calculated are consistent with the ones found in the `intron` table.
#   expect_equal(df_combined$TSL.x, df_combined$TSL.y)
# })

context("\tTest that the biotypes percentage values are corrently obtained")
test_that("Test that the biotypes percentage values are corrently obtained", {
  skip_if(!test_biotypes, "No biotypes percentage test. Variable test_biotypes set to FALSE in global options.")
  
  ## Load the biotype percentages
  df_all_percentage_tidy_merged = readRDS(all_annotated_SR_details_length_105_raw_biotype_path)
  
  
  ## Add the calculated percentages to the reference table
  df_combined <- df_intron %>% 
    left_join(df_all_percentage_tidy_merged,
              by = c("ref_coordinates" = "junID"))
  
  ## All `lncRNA` percentage values calculated are consistent with the ones found in the `intron` table.
  expect_equal(df_combined$lncRNA, df_combined$percent_lncRNA)
  
  ## All `protein_coding` percentage values calculated are consistent with the ones found in the `intron` table.
  expect_equal(df_combined$protein_coding, df_combined$percent_PC)
})

context("\tTest that the MANE information is properly assigned")
test_that("Test that the MANE information is properly assigned", {
  
  ## Extract MANE exons and convert them to introns
  df_mane_exons <- hg_mane %>% as_tibble() %>% dplyr::filter(type == "exon")
  df_mane_introns <- ggtranscript::to_intron(df_mane_exons, "transcript_name") %>%
    mutate(start = start + 1, end = end - 1) %>%
    mutate(intronID = paste0(seqnames, ":", start, "-", end, ":", strand)) %>%
    mutate(MANE = 1) %>%
    distinct(intronID, .keep_all = T)
  
  df_mane_introns <- df_mane_introns %>%
    mutate(transcript_id = str_sub(string = transcript_id,
                                   start = 1,
                                   end = 15)) 
  
  df_mane_introns <- df_mane_introns %>%
    inner_join(y = df_transcript,
              by = c("transcript_id" = "transcript_id"))
  
  df_mane_introns %>%
    dplyr::count(MANE.x)
  
  
  ## Connect to the database
  
  df_intron_tidy <- df_intron %>% 
    left_join(y = df_transcript %>% dplyr::select(id, MANE),
              by = c("transcript_id" = "id")) %>%
    left_join(y = df_mane_introns %>% 
                dplyr::select(id, MANE.x) %>% distinct(id, .keep_all = T),
              by = c("transcript_id" = "id")) %>%
    as_tibble()
  

  
  
  df_intron_tidy[c("MANE.x")][is.na(df_intron_tidy[c("MANE.x")])] <- 0
  
  df_intron_tidy %>%
    dplyr::count(MANE)
  df_intron_tidy %>%
    dplyr::count(MANE.x)
  
  expect_equal(df_intron_tidy %>%
                dplyr::count(MANE) %>% pull(n),
                df_intron_tidy %>%
                  dplyr::count(MANE.x) %>% pull(n))
  
  # df_intron_tidy %>%
  #   dplyr::count(MANE)
  # 
  # query = paste0("ALTER TABLE 'intron' ADD MANE NUMERIC NOT NULL DEFAULT 0");
  # DBI::dbSendQuery(con, query)
  # 
  # query = paste0("UPDATE 'intron' SET MANE = 1 WHERE ref_junID IN (", paste(df_intron_tidy %>%
  #                                                                             filter(MANE == 1) %>%
  #                                                                             pull(ref_junID), collapse = ","), ")")
  # DBI::dbSendQuery(con, query)
  # 
  # 
  # query = paste0("UPDATE 'intron' SET MANE = 0 WHERE ref_junID IN (", paste(df_intron_tidy %>%
  #                                                                             filter(MANE == 0) %>%
  #                                                                             pull(ref_junID), collapse = ","), ")")
  # DBI::dbSendQuery(con, query)
  # 
  
  ###################################################
})
# context("\tTest that the MANE information is properly calculated")
# test_that("Test that the MANE information is properly calculated", {
#   skip_if(!test_MANE, "No MANE test. Variable test_MANE set to FALSE in global options.")
#   
#   ## Load the MANE transcripts
#   hg_MANE <- rtracklayer::import(hg_mane_transcripts_path)
#   hg_MANE_transcripts <- hg_MANE %>% 
#     as_tibble() %>%
#     select(-source, -score, -phase, -gene_id, -gene_type, -tag, -protein_id, -db_xref, -transcript_type, -exon_id, -exon_number, -width) %>%
#     mutate(transcript_id = transcript_id %>% str_sub(1, 15)) %>%
#     tidyr::drop_na() %>%
#     filter(type == "transcript") %>%
#     distinct(transcript_id) %>%
#     pull(transcript_id)
#   
#   ## Load the introns with transcripts IDs
#   df_all_introns_introverse_tidy = readRDS(df_all_introns_introverse_tidy_path)
#   
#   ## Find the MANE introns
#   df_mane_introns <- df_all_introns_introverse_tidy %>%
#     tidyr::unnest(tx_id_junction) %>%
#     mutate(MANE = ifelse(tx_id_junction %in% hg_MANE_transcripts, T, F)) %>%
#     group_by(ref_junID) %>%
#     summarise(MANE = any(MANE))
#   
#   ## Combine with the intron table
#   df_combined <- df_intron %>% 
#     left_join(df_mane_introns,
#               by = c("ref_coordinates" = "ref_junID"))
#   
#   ## All `MANE` classification calculated is consistent with the ones found in the `intron` table
#   expect_equal(df_combined$MANE.x %>% as.logical(), df_combined$MANE.y)
# })

context("\tTest that the misspliced information is properly calculated")
test_that("Test that the misspliced information is properly calculated", {
  intron_ref_junID_misspliced <- df_intron %>% filter(misspliced == T) %>% pull(ref_junID) %>% sort()
  intron_ref_junID_nevermisspliced <- df_intron %>% filter(misspliced == F) %>% pull(ref_junID) %>% sort()
  
  ## Obtain the clusters
  clusters <- getClusters(con)
  
  ## Obtain the misspliced and the never misspliced junctions
  misspliced_junID <- c()
  nevermisspliced_junID <- c()
  for(cluster in clusters){
    misspliced_junID <- c(misspliced_junID, tbl(con, paste0(cluster, "_misspliced")) %>% 
                            pull(ref_junID) %>% unique) %>% unique
    nevermisspliced_junID <- c(nevermisspliced_junID, tbl(con, paste0(cluster, "_nevermisspliced")) %>% 
                                 pull(ref_junID) %>% unique) %>% unique
  }
  nevermisspliced_junID <- setdiff(nevermisspliced_junID, misspliced_junID)
  
  ## All misspliced junctions found in the child tables should be in the intron table as misspliced
  expect_equal(intron_ref_junID_misspliced, misspliced_junID %>% sort())
  
  ## All never misspliced junctions found in the child tables should be in the intron table as never misspliced
  expect_equal(intron_ref_junID_nevermisspliced, nevermisspliced_junID %>% sort())
})

clearAllVariables()

