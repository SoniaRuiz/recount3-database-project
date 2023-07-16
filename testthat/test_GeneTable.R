main_path <- normalizePath(path = "./")
source(paste0(main_path, "/testthat/helper_Files/helper_Global.R"))
skip_if(!test_GeneTable, "Gene table tests not executed. Variable test_GeneTable set to FALSE in global options.")
source(paste0(main_path, "/testthat/helper_Files/helper_GeneTable.R"))

context("\tTest that genes exists in the reference transcriptome")
test_that("Test that genes exists in the reference transcriptome", {
  reference_genes <- hg38 %>% as_tibble() %>% select(gene_id, gene_name) %>% distinct()
  df_combined <- df_gene %>% left_join(reference_genes, by = "gene_id")
  
  ## All gene IDs are found in the reference transcriptome
  expect_true(all(df_gene$gene_id %in% reference_genes$gene_id))
  
  ## All gene names are found in the reference transcriptome
  expect_true(all(df_gene$gene_name %in% reference_genes$gene_name))
  
  ## All gene names and IDs are consistent with the reference transcriptome
  expect_equal(df_combined$gene_name.x, df_combined$gene_name.y)
})

context("\tTest that all transcripts are correct")
test_that("Test that all transcripts are correct", {
  hg38_transcripts <- hg38 %>%
    as_tibble() %>%
    filter(gene_id %in% df_gene$gene_id & type == "transcript") %>%
    dplyr::count(gene_id, name = "n_transcripts")

  df_combined = df_gene %>% left_join(hg38_transcripts, by = "gene_id")

  ## All gene transcripts are consistent with the reference transcriptome
  expect_equal(df_combined$n_transcripts.x, df_combined$n_transcripts.y)
})

context("\tTest that gene widths are correct")
test_that("Test that gene widths are correct", {
  hg38_widths <- hg38 %>% 
    as_tibble() %>%
    filter(type == "gene") %>%
    select(width, gene_id)
  
  df_combined <- df_gene %>% left_join(hg38_widths, by = "gene_id")
  
  ## All gene widths are consistent with the reference transcriptome
  expect_equal(df_combined$gene_width, df_combined$width)
})

context("\tTest that no NAs are found in where they are not allowed")
test_that("Test that no NAs are found in where they are not allowed", {
  ## No single value should be NA or equal to "NA"
  expect_false(any(df_gene %>% select(-gene_name) %>% is.na()))
  expect_false(any((df_gene %>% select(-gene_name)) == "NA"))
})

clearAllVariables()

