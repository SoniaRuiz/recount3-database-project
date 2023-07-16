main_path <- normalizePath(path = "./")
source(paste0(main_path, "/testthat/helper_Files/helper_Global.R"))
skip_if(!test_ClusterTables, "Cluster tables tests not executed. Variable test_ClusterTables set to FALSE in global options.")
source(paste0(main_path, "/testthat/helper_Files/helper_ClusterTables.R"))

context("\tTest junction IDs and their reference with intron and novel tables")
test_that("Test junction IDs and their reference with intron and novel tables", {
  ## Loop for every valid cluster
  for(i in seq(length(test_clusters))){
    cluster <- test_clusters[i]
    df_misspliced <- tbl(con, paste0(cluster, "_misspliced")) %>% collect()
    df_never <- tbl(con, paste0(cluster, "_nevermisspliced")) %>% collect()
    
    ## Extract the reference introns IDs
    misspliced_ref_junID <- df_misspliced %>% pull(ref_junID)
    never_ref_junID <- df_never %>% pull(ref_junID)
    intron_ref_junID <- df_intron %>% pull(ref_junID)
    
    ## All reference introns are found in the intron table
    expect_true(all(misspliced_ref_junID %in% intron_ref_junID), info = paste0("Error in cluster ", cluster))
    expect_true(all(never_ref_junID %in% intron_ref_junID), info = paste0("Error in cluster: ", cluster))
    
    ## Extract the novel junction IDs
    misspliced_novel_junID <- df_misspliced %>% pull(novel_junID)
    novel_junID <- df_novel %>% pull(novel_junID)
    
    ## All novel junctions are found in the novel table
    expect_true(all(misspliced_novel_junID %in% novel_junID), info = paste0("Error in cluster ", cluster))
    
    ## All novel_junID and ref_junID pairs are also found in the novel table
    df_combined <- df_misspliced %>% 
      select(novel_junID, ref_junID) %>%
      left_join(df_novel %>% select(novel_junID, ref_junID),
                by = "novel_junID")
    
    expect_equal(df_combined$ref_junID.x, df_combined$ref_junID.y, info = paste0("Error in cluster ", cluster))
  }
})

context("\tTest reference between cluster and transcript table")
test_that("Test reference between cluster and transcript table", {
  ## Loop for every valid cluster
  for(i in seq(length(test_clusters))){
    cluster <- test_clusters[i]
    df_misspliced <- tbl(con, paste0(cluster, "_misspliced")) %>% collect()
    df_never <- tbl(con, paste0(cluster, "_nevermisspliced")) %>% collect()
    
    misspliced_transcript_id <- df_misspliced %>% pull(transcript_id)
    never_transcript_id <- df_never %>% pull(transcript_id)
    
    ## All gene_id in misspliced table are properly referenced to a row in the transcript table
    expect_true(all(misspliced_transcript_id %in% df_transcript$id))
    
    ## All gene_id in nevermisspliced table are properly referenced to a row in the transcript table
    expect_true(all(never_transcript_id %in% df_transcript$id))
  }
})

context("\tTest that all reads, introns and never mis-spliced information are properly extracted from files")
test_that("Test that all reads, introns and never mis-spliced information are properly extracted from files", {
  skip_if(!test_cluster_data, "No data verification for cluster tables. Variable test_cluster_data set to FALSE in global options.")
  ## Loop for every valid cluster
  for(i in seq(length(test_clusters))){
    cluster <- test_clusters[i]
    
    df_misspliced <- tbl(con, paste0(cluster, "_misspliced")) %>% collect()
    df_never <- tbl(con, paste0(cluster, "_nevermisspliced")) %>% collect()
    
    ## Load the Split Read counts and samples for the cluster
    split_read_counts_path <- generateClusterSplitReadsPath(cluster, results_folder, gtf_version, main_project)
    split_read_counts <<- getClusterSplitReads(split_read_counts_path)
    samples <- names(split_read_counts %>% select(-junID))
    
    ## Load the Annotated SR details
    annotated_SR_details_path <- generateAnnotatedSRdetailsPath(cluster, results_folder, 
                                                                gtf_version, main_project)
    annotated_SR_details <<- getAnnotatedSR(annotated_SR_details_path)
    
    ## Fix the * strand in the previous files
    fixUnknownStrand(split_read_counts, annotated_SR_details)
    
    #################### Test that all reads are properly extracted
    context(paste0("\t\tCluster ", cluster, ". Test that all reads are properly extracted"))
    ## Get the coordinates of every novel junction and reference intron
    df_merged_misspliced <- df_misspliced %>%
      left_join(df_novel %>% select(novel_junID, novel_coordinates),
                by = "novel_junID") %>%
      left_join(df_intron %>% select(ref_junID, ref_coordinates),
                by = "ref_junID")
    
    df_merged_never <- df_never %>%
      left_join(df_intron %>% select(ref_junID, ref_coordinates),
                by = "ref_junID")
    
    ## All novel and ref coordinates are found in the split reads
    novel_idx <- match(df_merged_misspliced$novel_coordinates, split_read_counts$junID)
    ref_idx_misspliced <- match(df_merged_misspliced$ref_coordinates %>% unique, split_read_counts$junID)
    ref_idx_never <- match(df_merged_never$ref_coordinates, split_read_counts$junID)
    
    expect_false(any(is.na(novel_idx)))
    expect_false(any(is.na(ref_idx_misspliced)))
    expect_false(any(is.na(ref_idx_never)))
    
    ## All novel_sum_counts and novel_n_individuals are properly calculated
    split_counts_novel <- split_read_counts[novel_idx, ] %>%
      mutate(sum = rowSums(across(where(is.numeric))),
             individuals = rowSums(across(where(is.numeric) & !sum) > 0)) %>%
      select(junID, sum, individuals)
    
    ## All novel junctions have at least 2 reads
    expect_false(any(split_counts_novel$sum == 1))
    expect_false(any(df_merged_misspliced$novel_sum_counts == 1))
    
    expect_equal(df_merged_misspliced$novel_sum_counts, split_counts_novel$sum)
    expect_equal(df_merged_misspliced$novel_n_individuals, split_counts_novel$individuals)
    
    ## All misspliced ref_sum_counts and ref_n_individuals are properly calculated
    split_counts_ref <- split_read_counts[ref_idx_misspliced %>% unique, ] %>%
      mutate(sum = rowSums(across(where(is.numeric))),
             individuals = rowSums(across(where(is.numeric) & !sum) > 0)) %>%
      select(junID, sum, individuals)
    
    
    ## All mis-spliced annotated introns have at least 2 reads
    expect_false(any(split_counts_ref$sum == 1))
    expect_false(any(df_merged_misspliced$ref_sum_counts == 1))
    
    expect_equal(df_merged_misspliced %>% distinct(ref_junID, .keep_all = T) %>% pull(ref_sum_counts), split_counts_ref$sum)
    expect_equal(df_merged_misspliced %>% distinct(ref_junID, .keep_all = T) %>% pull(ref_n_individuals), split_counts_ref$individuals)
    
    ## All never misspliced ref_sum_counts and ref_n_individuals are properly calculated
    split_counts_ref_never <- split_read_counts[ref_idx_never, ] %>%
      mutate(sum = rowSums(across(where(is.numeric))),
             individuals = rowSums(across(where(is.numeric) & !sum) > 0)) %>%
      select(junID, sum, individuals)
    
    ## All never mis-spliced annotated introns have at least 2 reads
    expect_false(any(split_counts_ref_never$sum == 1))
    expect_false(any(df_merged_never$ref_sum_counts == 1))
    
    expect_equal(df_merged_never$ref_sum_counts, split_counts_ref_never$sum)
    expect_equal(df_merged_never$ref_n_individuals, split_counts_ref_never$individuals)
    
    #################### Test that all junctions are properly obtained from annotated_SR_details
    context(paste0("\t\tCluster ", cluster, ". Test that all junctions are properly obtained from annotated_SR_details"))
    
    novel_coordinates <- df_novel %>% 
      filter(novel_junID %in% df_misspliced$novel_junID) %>% 
      pull(novel_coordinates)
    ref_coordinates_misspliced <- df_intron %>%
      filter(ref_junID %in% unique(df_misspliced$ref_junID)) %>%
      pull(ref_coordinates)
    ref_coordinates_never <- df_intron %>%
      filter(ref_junID %in% df_never$ref_junID) %>%
      pull(ref_coordinates)
    
    ## All novel junctions must be novel acceptor and novel donor in annotated SR details
    expect_true(all(novel_coordinates %in% (annotated_SR_details %>% 
                                              filter(type %in% c("novel_acceptor", "novel_donor")) %>%
                                              pull(junID))))
    
    ## All misspliced annotated introns must be "annotated" in annotated SR details
    expect_true(all(ref_coordinates_misspliced %in% (annotated_SR_details %>% 
                                                       filter(type == "annotated") %>%
                                                       pull(junID))))
    
    
    ## All nevermisspliced annotated introns must be "annotated" in annotated SR details
    expect_true(all(ref_coordinates_never %in% (annotated_SR_details %>% 
                                                  filter(type == "annotated") %>%
                                                  pull(junID))))
    
    
    #################### Test that never misspliced annotated introns are trully never misspliced
    context(paste0("\t\tCluster ", cluster, ". Test that never misspliced annotated introns are trully never misspliced"))
    
    ## Change annotated SR to accelerate the process
    annotated_SR_details <- annotated_SR_details %>%
      as_tibble() %>%
      select(junID, seqnames, start, end, width, strand, type)
    
    ## Calculate all distances for the cluster
    tmp_dir = tempdir()
    doParallel::registerDoParallel(num_cores)
    cluster_distances <- foreach(j = 1:length(samples), .combine = "rbind") %dopar% {
    #for(j in seq(length(samples))){
      #if(i %% 10 == 0) print(paste0(i, " of ", length(samples)))
      sample <- samples[j]
      
      ## Extract the counts for the particular sample
      split_read_counts_sample <- split_read_counts %>%
        dplyr::select(junID, !!sym(sample)) %>%
        dplyr::filter(!!sym(sample) != 0) %>%
        mutate(total=rowSums(select_if(., is.numeric))) %>% 
        filter(total >= 2) %>%
        dplyr::select(-total)
      
      ## Merge the annotated introns with the reads in the current sample
      annotated_SR_details_sample <- annotated_SR_details %>%
        inner_join(
          y = split_read_counts_sample,
          by = "junID"
        ) %>%
        dplyr::rename(counts = all_of(sample)) %>%
        GenomicRanges::GRanges()
      
      ## Annotated introns
      all_annotated <- annotated_SR_details_sample[annotated_SR_details_sample$type == "annotated"]
      all_annotated_forward <- all_annotated[all_annotated@strand == "+"]
      all_annotated_reverse <- all_annotated[all_annotated@strand == "-"]
      
      ## Novel donors
      all_donor <- annotated_SR_details_sample[annotated_SR_details_sample$type == "novel_donor"]
      all_donor_forward <- all_donor[all_donor@strand == "+"]
      all_donor_reverse <- all_donor[all_donor@strand == "-"]
      
      ## Novel acceptors
      all_acceptor <- annotated_SR_details_sample[annotated_SR_details_sample$type == "novel_acceptor"]
      all_acceptor_forward <- all_acceptor[all_acceptor@strand == "+"]
      all_acceptor_reverse <- all_acceptor[all_acceptor@strand == "-"]
      
      ## Generation of the distance dataframes
      df_nd_f <- getDistancesDataFrame(all_donor_forward, all_annotated_forward, "end", sample, "novel_donor")
      df_nd_r <- getDistancesDataFrame(all_donor_reverse, all_annotated_reverse, "start", sample, "novel_donor")
      df_na_f <- getDistancesDataFrame(all_acceptor_forward, all_annotated_forward, "start", sample, "novel_acceptor")
      df_na_r <- getDistancesDataFrame(all_acceptor_reverse, all_annotated_reverse, "end", sample, "novel_acceptor")
      
      sample_distances <- rbind(df_nd_f, df_nd_r, df_na_f, df_na_r)
      sample_distances %>% saveRDS(paste0(tmp_dir, "/distances_sample_", j, ".rds"))
      return()
    }
    
    cluster_distances <- purrr::map_df(1:length(samples), function(j){
      file_name <- paste0(tmp_dir, "/distances_sample_", j, ".rds")
      
      if(file.exists(file_name)){
        df <- readRDS(file_name)
        file.remove(file_name)
        return(df)
      }else{
        print("ERROR")
        break
      }
    })
    
    cluster_distances %>% head
    miss_spliced_junctions <- cluster_distances %>%
      pull(ref_junID) %>%
      unique()
    
    ## List of all annotated introns
    all_junctions <- annotated_SR_details %>%
      filter(type == "annotated" & strand != "*") %>%
      pull(junID) %>%
      unique()
    
    ## Difference between all the annotated introns and the mis-spliced introns
    never_miss_spliced_junctions <- setdiff(all_junctions, miss_spliced_junctions)
    
    ## Obtaining the coordinates of every never misspliced intron in the cluster
    df_never_coordinates <- df_never %>% 
      left_join(df_intron %>% 
                  select(ref_junID, ref_coordinates), 
                by = "ref_junID")
    
    ## All never misspliced introns from the cluster are not found as misspliced
    ## when calculating the distances again.
    expect_true(all(df_never_coordinates$ref_coordinates %in% never_miss_spliced_junctions))
    
    # #################### Test that TPM values are properly calculated
    # context(paste0("\t\tCluster ", cluster, ". Test that TPM values are properly calculated"))
    # 
    # ## Get the TPM dataframe and calculate the median of each row
    # tpm_path <-  generateTPMpath(cluster, results_folder, gtf_version, main_project)
    # tpm <- getClusterTPM(tpm_path) %>%
    #   select(gene_id = gene, all_of(samples)) %>%
    #   mutate(tpm_median = matrixStats::rowMedians(as.matrix(.[-1]))) %>%
    #   select(gene_id, tpm_median) %>% 
    #   group_by(gene_id) %>% 
    #   summarize_all(max)
    # 
    # ## Merge the calculated tpm values with the database information
    # df_genes_tpm <- df_gene %>% 
    #   select(id, gene_id) %>%
    #   left_join(tpm, by = "gene_id")
    #   
    # df_tpm_misspliced <- df_misspliced %>%
    #   select(id = gene_id, gene_tpm) %>%
    #   left_join(df_genes_tpm, by = "id")
    # 
    # df_tpm_never <- df_never %>%
    #   select(id = gene_id, gene_tpm) %>%
    #   left_join(df_genes_tpm, by = "id")
    # 
    # ## All gene_tpm values should match the obtained from the tpm results
    # expect_equal(df_tpm_misspliced$gene_tpm, df_tpm_misspliced$tpm_median)
    # expect_equal(df_tpm_never$gene_tpm, df_tpm_never$tpm_median)
    
    ## Delete the cluster variables after execution
    rm(split_read_counts, annotated_SR_details, envir = .GlobalEnv)
    rm(list = ls())
  }
})

context("\tTest that Mis-Splicing Ratio is properly calculated")
test_that("Test that Mis-Splicing Ratio is properly calculated", {
  ## Loop for every valid cluster
  for(i in seq(length(test_clusters))){
    cluster <- test_clusters[i]
    df_misspliced <- tbl(con, paste0(cluster, "_misspliced")) %>% collect()
    df_never <- tbl(con, paste0(cluster, "_nevermisspliced")) %>% collect()
    
    ## Get all reference introns and novel junctions
    df_MSR <- df_misspliced %>% 
      left_join(df_novel %>% select(novel_junID, novel_type),
                by = "novel_junID") %>% 
      group_by(ref_junID, novel_type) %>%
      mutate(MSR = sum(novel_sum_counts)/(sum(novel_sum_counts) + ref_sum_counts)) %>%
      tidyr::spread(key = novel_type, value = MSR, fill = 0) %>%
      group_by(ref_junID) %>%
      mutate(MSR_Donor = max(novel_donor, na.rm = T),
             MSR_Acceptor = max(novel_acceptor, na.rm = T)) %>%
      select(-novel_donor, -novel_acceptor) %>%
      ungroup()
    
    ## All MSR ratios are properly calculated
    expect_equal(df_MSR$MSR_A, df_MSR$MSR_Acceptor, info = paste0("Error in cluster ", cluster))
    expect_equal(df_MSR$MSR_D, df_MSR$MSR_Donor, info = paste0("Error in cluster ", cluster))
    
    ## Test the consistency of the values
    MSR_D <- df_misspliced %>% pull(MSR_D)
    MSR_A <- df_misspliced %>% pull(MSR_A)
    
    MSR_D_never <- df_never %>% pull(MSR_D)
    MSR_A_never <- df_never %>% pull(MSR_A)
    
    ## All misspliced MSR are between 0 and 1
    expect_true(all(MSR_D %>% dplyr::between(0, 1)), info = paste0("Error in cluster ", cluster))
    expect_true(all(MSR_A %>% dplyr::between(0, 1)), info = paste0("Error in cluster ", cluster))
    
    ## No misspliced MSR at 1
    expect_false(any(MSR_D == 1), info = paste0("Error in cluster ", cluster))
    expect_false(any(MSR_A == 1), info = paste0("Error in cluster ", cluster))
    
    ## Never misspliced MSR are equal to 0
    expect_false(any(MSR_D_never != 0), info = paste0("Error in cluster ", cluster))
    expect_false(any(MSR_A_never != 0), info = paste0("Error in cluster ", cluster))
  }
})

context("\tTest that reference intron type is consistent with MSR")
test_that("Test that reference intron type is consistent with MSR", {
  ## Loop for every valid cluster
  for(i in seq(length(test_clusters))){
    cluster <- test_clusters[i]
    df_misspliced <- tbl(con, paste0(cluster, "_misspliced")) %>% collect()
    df_never <- tbl(con, paste0(cluster, "_nevermisspliced")) %>% collect()
    
    df_ref_type <- df_misspliced %>% 
      rowwise %>%
      mutate(test_ref_type = missplicingClass(MSR_D, MSR_A))
    
    ## All ref_types are as expected
    expect_equal(df_ref_type$ref_type, df_ref_type$test_ref_type)
    expect_true(all(df_never$ref_type == "never"))
  }
})

# context("\tTest that no NAs are found in where they are not allowed")
# test_that("Test that no NAs are found in where they are not allowed", {
#   ## Loop for every valid cluster
#   for(i in seq(length(test_clusters))){
#     cluster <- test_clusters[i]
#     df_misspliced <- tbl(con, paste0(cluster, "_misspliced")) %>% collect()
#     df_never <- tbl(con, paste0(cluster, "_nevermisspliced")) %>% collect()
#     
#     ## No single value should be NA or equal to "NA"
#     expect_false(any(df_misspliced %>% select(-gene_tpm) %>% is.na()))
#     expect_false(any((df_misspliced %>% select(-gene_tpm)) == "NA"))
#     
#     expect_false(any(df_never %>% select(-gene_tpm) %>% is.na()))
#     expect_false(any((df_never %>% select(-gene_tpm)) == "NA"))
#   }
# })

context("\tTest that the set of introns in intron table and the set in the cluster tables are the same")
test_that("Test that the set of introns in intron table and the set in the cluster tables are the same", {
  ## Obtain the misspliced and the never misspliced introns
  misspliced_junID <- c()
  nevermisspliced_junID <- c()
  for(cluster in getClusters(con)){
    misspliced_junID <- c(misspliced_junID, tbl(con, paste0(cluster, "_misspliced")) %>% pull(ref_junID) %>% unique) %>% unique
    nevermisspliced_junID <- c(nevermisspliced_junID, tbl(con, paste0(cluster, "_nevermisspliced")) %>% pull(ref_junID) %>% unique) %>% unique
  }
  nevermisspliced_junID <- setdiff(nevermisspliced_junID, misspliced_junID)
  
  ## All misspliced introns found in the cluster tables should be in the intron table as misspliced
  expect_equal(df_intron %>% filter(misspliced == T) %>% pull(ref_junID) %>% sort(), misspliced_junID %>% sort())
  
  ## All never misspliced introns found in the cluster tables should be in the intron table as never misspliced
  expect_equal(df_intron %>% filter(misspliced == F) %>% pull(ref_junID) %>% sort(), nevermisspliced_junID %>% sort())
})

clearAllVariables()
