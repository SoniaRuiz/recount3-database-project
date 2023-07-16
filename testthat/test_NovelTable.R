main_path <- normalizePath(path = "./")
source(paste0(main_path, "/testthat/helper_Files/helper_Global.R"))
skip_if(!test_NovelTable, "Novel table tests not executed. Variable test_NovelTable set to FALSE in global options.")
source(paste0(main_path, "/testthat/helper_Files/helper_NovelTable.R"))

context("\tTest consistent junction IDs with other tables")
test_that("Test consistent junction IDs with other tables", {
  novel_junID <- df_novel %>% pull(novel_junID)
  novel_ref_junID <- df_novel %>% pull(ref_junID)
  intron_ref_junID <- df_intron %>% pull(ref_junID)
  intron_ref_junID_misspliced <- df_intron %>% filter(misspliced == T) %>% pull(ref_junID)
  
  ## All `ref_junID` from the novel table are found in the `ref_junID` from intron table.
  expect_true(all(novel_ref_junID %in% intron_ref_junID))
  
  ## All `ref_junID` from the novel table are found in the mis-spliced `ref_junID` from intron table.
  expect_true(all(novel_ref_junID %in% intron_ref_junID_misspliced))
  
  ## No duplicates found in `novel_junID`
  expect_equal(novel_junID %>% unique() %>% length, novel_junID %>% length())
})

context("\tTest that different fields have consistent values")
test_that("Test that different fields have consistent values", {
  novel_type <- df_novel %>% pull(novel_type) %>% unique() %>% sort()
  
  ## All `novel_types` are found in the valid types
  expect_true(all(novel_type %in% valid_novel_types))
})

context("\tTest novel_coordinates field")
test_that("Test novel_coordinates field", {
  ## We use a custom function to extract the information from the coordinates
  GRdata <- getGRdata(df_novel$novel_junID, df_novel$novel_coordinates, fix_data_types = F, remove_chr = F)
  novel_coordinates <- df_novel$novel_coordinates
  ref_coordinates <- df_intron$ref_coordinates
  
  novel_seqnames <- GRdata$seqnames
  novel_strand <- GRdata$strand
  novel_start = GRdata$start
  novel_end = GRdata$end
  
  ## No duplicates found in `novel_coordinates`
  expect_equal(novel_coordinates %>% length(), novel_coordinates %>% unique() %>% length())
  
  ## No `novel_coordinates` value is found in `ref_coordinates`
  expect_length(intersect(novel_coordinates, ref_coordinates), n = 0)
  
  ## No `*` strand found in `novel_coordinates`
  expect_equal(valid_novel_strands, novel_strand %>% unique() %>% sort()) 
  
  ## All `seqnames` in `novel_coordinates` are within the 1-22, X and Y chromosomes
  expect_true(all((novel_seqnames %>% unique()) %in% valid_novel_seqnames))
  
  ## All `start` and `end` positions are integers (or numeric)
  expect_error(as.numeric(novel_start), NA)
  expect_error(as.numeric(novel_end), NA)
})

context("\tTest that no NAs are found in where they are not allowed")
test_that("Test that no NAs are found in where they are not allowed", {
  ## No single value should be NA or equal to "NA"
  expect_false(any(df_novel %>% is.na()))
  expect_false(any(df_novel == "NA"))
})

context("\tTest that MaxEntScan scores are correctly calculated")
test_that("Test that MaxEntScan scores are correctly calculated", {
  skip_if(!test_MaxEntScan, "No MaxEntScan test. Variable test_MaxEntScan set to FALSE in global options.")
  
  GRdata <- getGRdata(df_novel$novel_junID, df_novel$novel_coordinates)
  GRdata_MaxEntScan <- generateMaxEntScore(GRdata, 
                                           bedtools_path = bedtools_path, 
                                           fasta_path = fasta_path, 
                                           fordownload_path = fordownload_path)
  
  ## Both tables have similar number of rows
  expect_equal(nrow(df_novel), nrow(GRdata_MaxEntScan))
  
  ## They share the same junIDs
  expect_equal(df_novel$novel_junID, GRdata_MaxEntScan$junID)
  
  ## The calculated MaxEntScore are consistent
  expect_equal(df_novel$novel_ss5score, GRdata_MaxEntScan$ss5score)
  expect_equal(df_novel$novel_ss3score, GRdata_MaxEntScan$ss3score)
})

context("\tDistances and coordinates betweeen novel and reference")
test_that("Distances and coordinates betweeen novel and reference", {
  df_novel_expanded <- getGRdata(df_novel$novel_junID, df_novel$novel_coordinates)
  df_novel_expanded <- df_novel_expanded %>% 
    dplyr::rename("novel_junID" = "junID",
                  "novel_seqnames" = "seqnames",
                  "novel_start" = "start",
                  "novel_end" = "end",
                  "novel_strand" = "strand") %>%
    left_join(df_novel %>% 
                select(novel_junID, ref_junID, novel_type, distance),
              by = "novel_junID")
  
  df_intron_expanded <- getGRdata(df_intron$ref_junID, df_intron$ref_coordinates) %>%
    dplyr::rename("ref_junID" = "junID",
                  "ref_seqnames" = "seqnames",
                  "ref_start" = "start",
                  "ref_end" = "end",
                  "ref_strand" = "strand")
  
  df_combined <- df_novel_expanded %>% 
    left_join(df_intron_expanded, by = "ref_junID")
  
  ## The strands of the reference junction in `novel` table and `intron` table are the same
  expect_equal(df_combined$novel_strand, df_combined$ref_strand)
  
  ## The `novel_type` is consistent with the fixed starting and end position of the novel junction and reference intron
  expect_true(df_combined %>% filter(novel_type == "novel_donor" & novel_strand == "+") %>% mutate(test = (novel_end == ref_end)) %>% pull(test) %>% all())
  expect_true(df_combined %>% filter(novel_type == "novel_donor" & novel_strand == "-") %>% mutate(test = (novel_start == ref_start)) %>% pull(test) %>% all())
  expect_true(df_combined %>% filter(novel_type == "novel_acceptor" & novel_strand == "+") %>% mutate(test = (novel_start == ref_start)) %>% pull(test) %>% all())
  expect_true(df_combined %>% filter(novel_type == "novel_acceptor" & novel_strand == "-") %>% mutate(test = (novel_end == ref_end)) %>% pull(test) %>% all())
  
  df_combined <- df_combined %>%
    dplyr::rowwise() %>%
    dplyr::mutate(calculated_distance = measureDistance(novel_start, novel_end, ref_start, ref_end, novel_type, novel_strand))
  
  ## The distances are properly calculated
  expect_equal(df_combined$distance, df_combined$calculated_distance)
})

clearAllVariables()