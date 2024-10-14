##############################################################
## CODE Adapted from:
## https://github.com/guillermo1996/ENCODE_Metadata_Extraction
##############################################################



#' Title
#' Downloads the metadata from each shRNA RBP knockdown experiment from the ENCODE platform
#' @param experiment_type 
#' @param results_path 
#' @param dependencies_path 
#' @param required_cell_lines 
#' @param valid_genome_annotation 
#' @param valid_file_format 
#' @param valid_output_type 
#' @param valid_nucleic_acid_type 
#' @param download_method 
#'
#' @return
#' @export
#'
#' @examples
ENCODEDownloadMetadata <- function(experiment_type,
                                   results_path,
                                   dependencies_path,
                                   required_cell_lines = c("HepG2", "K562"),
                                   valid_genome_annotation = "V29",
                                   valid_file_format = "bam",
                                   valid_output_type = "alignments",
                                   valid_nucleic_acid_type = "polyadenylated mRNA",
                                   #download_method = "experiments",
                                   download_method ="gene_silencing_series") {
  
  
  
  
  ## Valid or required values ----
  
  
  #valid_file_format <- "alignments"
  #valid_output_type <- "alignments"
  # valid_output_type <- c("minus strand signal of unique reads",
  #                        "plus strand signal of unique reads")
  # valid_nucleic_acid_type <- "polyadenylated mRNA"
  # overwrite_db <- T
  # #download_method = "experiments"
  # download_method ="gene_silencing_series"
  
  
  
  ### Output files
  output_json <- file.path(results_path, "response.json")
  output_search <- file.path(results_path, paste0("all_", experiment_type, "_experiments.tsv"))
  output_metadata <- file.path(results_path, paste0("metadata_", experiment_type, "_samples.tsv"))
  
  
  
  
  ### LOAD RBP NAMES
  
  # input_target_gene_categories <- readr::read_delim(file = paste0(here::here(), "/Additional_Files/Target_gene_categories.tsv"), show_col_types = F) %>% 
  #   dplyr::select(name, `Splicing regulation`, Spliceosome, `Novel RBP`, `Exon Junction Complex`) %>% 
  #   filter_at(vars(c(`Splicing regulation`, Spliceosome, `Novel RBP`, `Exon Junction Complex`)),  any_vars(. > 0))
  # 
  # input_target_gene_NMD <-  readr::read_delim(file = paste0(here::here(), "/Additional_Files/NMD.txt"), show_col_types = F, delim = "\n")%>% 
  #   dplyr::rename("Name" = "Name\t") %>% 
  #   dplyr::mutate(across(Name, str_replace, "\t", "")) %>% 
  #   dplyr::pull(Name)
  # 
  # 
  # valid_target_genes <- c(input_target_gene_categories$name, input_target_gene_NMD) %>% unique() %>% sort()
  
  # 2. Pipeline ----
  
  ## Download the search data ----
  URL = "https://www.encodeproject.org/search/?type=Experiment&assay_title=shRNA+RNA-seq&control_type!=*&status=released&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&limit=all&format=json"
  
  ## LONG READ
  # URL = "https://www.encodeproject.org/search/?type=Experiment&control_type!=*&assay_term_name=long+read+RNA-seq&status=released&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&files.platform.term_name=Pacific+Biosciences+Sequel+II&replicates.library.nucleic_acid_term_name=polyadenylated+mRNA&assembly=GRCh38&limit=all&format=json"
  ## CRISPR
  #URL = "https://www.encodeproject.org/search/?type=GeneSilencingSeries&searchTerm=gene+silencing&target.investigated_as=RNA+binding+protein&target.investigated_as=transcription+factor&assay_term_name=CRISPR+genome+editing+followed+by+RNA-seq&related_datasets.replicates.library.biosample.applied_modifications.method=CRISPR&organism.scientific_name=Homo+sapiens&assembly=GRCh38&related_datasets.files.file_type=bam&related_datasets.files.file_type=bigWig&award.rfa=ENCODE4&status=released&limit=all&format=json"
  ## shRNA
  
  
  
  ## Query the ENCODE API to download a list of experiments in json format ---
  response_data <- GetUrlResponse(URL, output_file = output_json)
  
  
  
  
  ## Summarize the list of experiments from ENCODE search --------------------
  summary_df <- GenerateSummary(response_data, 
                                valid_target_genes = NULL,
                                output_file = output_search)
  
  
  
  
  
  ## Extract the metadata ----
  metadata_df <- GenerateMetadata(summary_df, 
                                  download_method = download_method,
                                  required_cell_lines = required_cell_lines,
                                  valid_file_format = valid_file_format,
                                  valid_genome_annotation = valid_genome_annotation,
                                  valid_output_type = valid_output_type,
                                  valid_nucleic_acid_type = valid_nucleic_acid_type,
                                  output_file = output_metadata,
                                  overwrite_db = T)
  
  
  ## Add category information ----
  if ( str_detect(string = metadata_df[1,]$assay, pattern = "shRNA", negate = F) ) {
    
    
    metadata_df <- AddTargetGeneCategory(metadata_df = output_metadata,
                                         input_Category = file.path(dependencies_path, "RBPs_subgroups.xlsx"),
                                         input_NMD = file.path(dependencies_path, "NMD.txt"),
                                         output_file = output_metadata)
    
    metadata_df <- metadata_df %>%
      dplyr::filter(if_any(c("Splicing.regulation", Spliceosome, "Novel.RBP", "Exon.Junction.Complex", NMD), ~ . != 0))
    
  }
  
  return(metadata_df)
  
}



## HELPER FUNCTIONS ----------------------------------------------------------------------


#' Extract API response
#'
#' Given an URL, it requests the information from the ENCODE API.
#'
#' @param url Link to the ENCODE page to extract the information.
#' @param output_file (Optional) Path to where the json downloaded file will be
#'   stored.
#'
#' @return A list object containing the json downloaded from the ENCODE portal.
#' @export
GetUrlResponse <- function(url, output_file = NULL) {
  response <- httr::GET(url)
  r <- httr::content(response, as = "text", encoding = "UTF-8")
  
  if(!is.null(output_file)){
    write(jsonlite::prettify(r, indent = 4), output_file)
  }
  
  return(jsonlite::fromJSON(r, flatten = T))
}

#' Generates a summary of the ENCODE experiment search
#'
#' From the generated list object built in [GetUrlResponse()], it reads every
#' experiment found and extract their experiment ID, Gene Silencing Series (gss)
#' ID and cell line. All this information is summarized in a dataframe, where
#' every row is an experiment.
#'
#' @param response_data List object build from [GetUrlResponse()] and a search
#'   URL.
#' @param output_file (Optional) Path to where the json downloaded file will be
#'   stored.
#'
#' @return Data.frame containing a summary of every experiment found in the
#'   ENCODE search.
#' @export
GenerateSummary <- function(response_data,
                            valid_target_genes = c(),
                            output_file = NULL) {
  summary_df <- tibble()
  
  for (i in seq(nrow(response_data$`@graph`))) {
    # i <- 1
    
    experiment <- response_data$`@graph`[i, ]
    
    if ( experiment$assay_term_name[[1]] == "CRISPR genome editing followed by RNA-seq" ) {
      
      biosample_ontology <- experiment$biosample_ontology[[1]]$term_name
      target_gene <- experiment$target[[1]]$label
      gene_silencing_series <- experiment$accession
      
      ## Extract the Gene Silencing Series
      sample_id <- sapply(experiment$related_datasets[[1]]$`@id`, function(x) str_split(string = x, pattern = "/")[[1]][3])
      
      
      tmp_df <- tibble::tibble(target_gene = target_gene,
                               experiment_id = sample_id,
                               cell_line = biosample_ontology,
                               gene_silencing_series = gene_silencing_series)
      summary_df <- rbind(summary_df, tmp_df)
      
    } 
    else if ( experiment$assay_term_name[[1]] == "long read RNA-seq" ) {
      
      ## Only process if not audit warnings
      if ( is.null(experiment$audit.NOT_COMPLIANT[[1]]) ) {
        
        
        biosample_ontology <- experiment$biosample_ontology.term_name
        type <- experiment$biosample_ontology.classification
        gene_silencing_series <- experiment$accession
        
        ## Extract the Gene Silencing Series
        sample_id <- sapply(experiment$files[[1]]$`@id`, function(x) str_split(string = x, pattern = "/")[[1]][3])
        
        
        tmp_df <- tibble::tibble(target_gene = biosample_ontology,
                                 experiment_id = sample_id,
                                 cell_line = type,
                                 gene_silencing_series = gene_silencing_series)
        summary_df <- rbind(summary_df, tmp_df)  
        
        
      }
      
    } else {
      
      biosample_ontology <- experiment$biosample_ontology.term_name
      target_gene <- experiment$target.label
      sample_id <- experiment$accession
      
      ## Extract the Gene Silencing Series
      gene_silencing_series <- sapply(experiment$related_series, function(x) x$accession)
      
      if(length(gene_silencing_series) > 1){
        logger::log_warn("More than one gene silencing series found for target gene ", target_gene, " and biosample ", biosample_ontology)
        gene_silencing_series <- gene_silencing_series[which(sapply(experiment$related_series[[1]]$`@type`, 
                                                                    function(x) "GeneSilencingSeries" %in% x))] 
      }
      
      tmp_df <- tibble::tibble(target_gene = target_gene,
                               experiment_id = sample_id,
                               cell_line = biosample_ontology,
                               gene_silencing_series = gene_silencing_series)
      summary_df <- rbind(summary_df, tmp_df)
    }
    
  }
  
  summary_df <- summary_df %>% 
    dplyr::arrange(target_gene, cell_line)
  
  summary_df %>%
    dplyr::count(target_gene,cell_line) %>%
    as.data.frame()
  
  if ( experiment$assay_term_name[[1]] != "long read RNA-seq" && 
       is.null(valid_target_genes) ) {
    
    # valid_target_genes <- summary_df %>%
    #   dplyr::count(target_gene) %>%
    #   filter(n == 4) %>%
    #   pull(target_gene)
  }
  
  
  if( !is.null(valid_target_genes) && length(valid_target_genes) > 0 ){
    summary_df <- summary_df %>% dplyr::filter(target_gene %in% valid_target_genes)
  }
  
  
  if( !is.null(output_file) ){
    write.table(summary_df %>%
                  tidyr::unnest(gene_silencing_series, keep_empty = T), 
                output_file, sep = "\t", row.names = F)
  }
  
  return(summary_df)
}


#' Generates the metadata data.frame
#'
#' Given a summary data.frame with the experiments information and the filters
#' to apply in the extraction process, it generates a data.frame with each
#' sample as a row and the sample's metadata as the columns.
#'
#' @param summary_df Data.frame containing a summary of every experiment found
#'   in the ENCODE search.
#' @param download_method (Optional) Whether to use the Gene Silencing Series
#'   (gss) ID or experiment ID to access the experiment's information. If
#'   available, the gss method is recomended since it requires less calls to the
#'   API. Valid inputs: "gene_silencing_series", "experiments". Defaults to
#'   "gene_silencing_series".
#' @param required_cell_lines (Optional) Required cell lines to extract
#'   information about the experiments. It requires that a least one experiment
#'   for each cell line provided is present in the summary data.frame. Defaults
#'   to "HepG2" and "K562".
#' @param valid_nucleic_acid_type (Optional) Required nucleic acid type of the
#'   experiment. Defaults to "polyadenylated mRNA".
#' @param valid_genome_annotation (Optional) Required gene annotation version of
#'   the sample to extract its metadata. Defaults to "V29".
#' @param valid_file_format (Optional) Required output file format of the
#'   sample. Defaults to "BAM".
#' @param valid_output_type (Optional) Required output type of the sample.
#'   Defaults to "alignments".
#' @param output_file (Optional) Path to where the json downloaded file will be
#'   stored.
#' @param overwrite_db (Optional) If set to TRUE, it will read first the
#'   metadata found in "output_file" and will try to not redownload the metadata
#'   of previously successfull runs.
#'
#' @return Data.frame where every row is an ENCODE sample and the columns
#'   contain their metadata.
#' @export
GenerateMetadata <- function(summary_df,
                             download_method,
                             required_cell_lines = NULL, 
                             valid_nucleic_acid_type = "polyadenylated mRNA",
                             valid_genome_annotation = "V29",
                             valid_file_format = "bam",
                             valid_output_type = "alignments",
                             output_file = "",
                             overwrite_db = T){
  
  ## Check the overwrite_db configuration. If set to FALSE, we will try to read the
  ## output file and extract the previous information from it.
  if(output_file == "") overwrite_db = T
  if(!overwrite_db & file.exists(output_file)){
    df_previous <- read.csv(output_file, sep = "\t") #%>% 
    #  tibble::as_tibble() %>% 
  }else{
    df_previous <- tibble() 
    overwrite_db = T
  }
  
  ## Get the valid target genes by cell lines
  valid_target_genes <- summary_df %>% 
    dplyr::group_by(target_gene) %>%
    dplyr::summarise(cell_lines = list(cell_line))
  
  if ( !is.null(required_cell_lines) ) {
    valid_target_genes <- valid_target_genes %>%
      dplyr::rowwise() %>%
      dplyr::filter(all(required_cell_lines %in% cell_lines)) %>%
      dplyr::pull(target_gene)
  } else {
    valid_target_genes <- valid_target_genes %>%
      dplyr::pull(target_gene)
  }
  
  metadata_df <- tibble()
  
  ## Loop through every target gene
  for(iter_tg in seq(length(valid_target_genes))){
    
    # iter_tg <- 1
    iter_target_gene <- valid_target_genes[iter_tg]
    logger::log_info("Starting target ", iter_target_gene, ":")
    
    ## Local target gene information
    summary_tg_df <- summary_df %>% dplyr::filter(target_gene == iter_target_gene)
    
    ## If we find a total of 8 entries in the previous search, we skip this
    ## target gene. 
    if(!overwrite_db){
      if (df_previous %>% dplyr::filter(target_gene == iter_target_gene) %>% nrow() == 8){
        logger::log_info("\t Already found in output file.")
        metadata_df <- rbind(metadata_df, df_previous %>% dplyr::filter(target_gene == iter_target_gene))
        next
      }
    }
    
    ## Download the target gene metadata
    if (download_method == "gene_silencing_series") {
      tg_df <- LoopGeneSilencingSeries(summary_tg_df,
                                       valid_nucleic_acid_type,
                                       valid_file_format,
                                       valid_output_type,
                                       valid_genome_annotation)
    }else if(download_method == "experiments"){
      tg_df <- LoopExperiments(summary_tg_df,
                               valid_nucleic_acid_type,
                               valid_file_format,
                               valid_output_type,
                               valid_genome_annotation)
    }else{
      logger::ERROR("No valid download method provided. Only gene_silencing_series or experiments are allowed.")
    }
    
    if (tg_df %>% nrow == 0) {
      next;
    }
    
    ## Add final information and sort the columns
    tg_df <- tg_df %>%
      dplyr::mutate(target_gene = iter_target_gene, .before = cell_line) %>%
      dplyr::filter(nucleic_acid_type == valid_nucleic_acid_type)
    
    column_order <- c("target_gene", "experiment_type", "cell_line", 
                      "gene_silencing_series", "experiment_id", 
                      "experiment_doi", "sample_id", "rin", "read_depth", 
                      "bio_rep", "tech_rep", "sex", "age", "life_stage", 
                      "gene_quantification_id", "file_format", 
                      "output_type", "genome_annotation", 
                      "mapped_run_type", "lab", "assay", "cellosaurus")
    
    tg_df <- tg_df %>% 
      dplyr::select(c(intersect(column_order, names(.)), setdiff(names(.), column_order)))
    
    ## If the required cell lines are not found, ignore the target gene
    #if(!all.equal(tg_df$cell_line %>% unique %>% sort, required_cell_lines %>% sort)){
    #  next
    #}
    
    ## Appends the data and store to disk
    metadata_df <- rbind(metadata_df, tg_df)
    if(output_file != ""){
      write.table(metadata_df, output_file, sep = "\t", row.names = F, quote = FALSE)
    }
  }
  
  return(metadata_df)
}

#' Extract metadata from the Gene Silencing Series (gss)
#'
#' @param summary_tg_df Local target gene information from the summary
#'   data.frame generated with [GenerateSummary()].
#' @param valid_nucleic_acid_type (Optional) Required nucleic acid type of the
#'   experiment. Defaults to "polyadenylated mRNA".
#' @param valid_file_format (Optional) Required output file format of the
#'   sample. Defaults to "BAM".
#' @param valid_output_type (Optional) Required output type of the sample.
#'   Defaults to "alignments".
#' @param valid_genome_annotation (Optional) Required gene annotation version of
#'   the sample to extract its metadata. Defaults to "V29".
#'
#' @return Data.frame containing the metadata for a particular gss.
#' @export
LoopGeneSilencingSeries <- function(summary_tg_df,
                                    valid_nucleic_acid_type,
                                    valid_file_format,
                                    valid_output_type,
                                    valid_genome_annotation){
  
  ## Loop through every gene silencing series found for the target gene
  tg_df <- foreach(iter_gss = seq(nrow(summary_tg_df)), .combine = dplyr::bind_rows) %do%{
    
    gss <- summary_tg_df[iter_gss, ]
    
    ## Gene silencing series information
    gss_id <- gss$gene_silencing_series
    gss_cell_line <- gss$cell_line
    logger::log_info("\t Starting Gene Silencing Series ", gss_id, " (cell line = ", gss_cell_line, "):")
    
    ## Get the API response for the gene silencing series
    response_gss <- GetUrlResponse(paste0("https://www.encodeproject.org/gene-silencing-series/", gss_id, "?format=json"))
    
    ## Information about the case and control samples
    gss_experiments <- response_gss$related_datasets
    
    ## Loop through every experiment found in the gene silencing series
    gss_df <- foreach(iter_experiment = seq(nrow(gss_experiments)), .combine = dplyr::bind_rows) %do%{
      
      experiment <- gss_experiments[iter_experiment, ]
      
      experiment_id <- experiment$accession
      experiment_type <- GetExperimentType(experiment, gss$target_gene, gss_cell_line) 
      
      logger::log_info("\t\t Starting experiment ", experiment_id, " (type = ", experiment_type, ").")
      
      ## Requirements metadata
      experiment_additional_info <- GetAdditionalInformation(experiment)
      
      ## Main metadata
      experiment_sample_files <- GetSampleFiles(files = experiment$files[[1]],
                                                valid_file_format = valid_file_format, 
                                                valid_output_type = valid_output_type, valid_genome_annotation)
      
      if ( nrow(experiment_sample_files) == 0 ) {
        next;
      }
      
      if(experiment_additional_info$nucleic_acid_type != valid_nucleic_acid_type) 
        return(cbind(experiment_sample_files, experiment_additional_info))
      
      ## Other metadata
      experiment_rin <- GetRin(experiment$replicates[[1]])
      experiment_read_depth <- GetReadDepth(experiment$analyses[[1]], valid_genome_annotation)
      experiment_donor_info <- GetDonorInfo(experiment$replicates[[1]])
      experiment_documents <- GetDocumentFiles(experiment$replicates[[1]], gss$target_gene, experiment_type)
      experiment_doi <- experiment$doi
      experiment_gene_quantifications <- GetGeneQuantificationFiles(experiment$files[[1]])
      
      ## Combine all information
      experiment_combined <- experiment_sample_files %>% 
        dplyr::left_join(experiment_rin, by = "bio_rep") %>%
        dplyr::left_join(experiment_read_depth, by = "bio_rep") %>%
        dplyr::left_join(experiment_donor_info, by = "bio_rep") %>% 
        dplyr::left_join(experiment_gene_quantifications %>% dplyr::select(gene_quantification_id, bio_rep), by = "bio_rep") %>%
        dplyr::left_join(experiment_documents, by = "bio_rep") %>%
        dplyr::mutate(experiment_type = experiment_type, experiment_id = experiment_id, experiment_doi = experiment_doi, 
                      .before = "bio_rep") %>%
        cbind(experiment_additional_info)
      
      return(experiment_combined)
    }
    
    ## Add the relevant information
    
    if (nrow(gss_df) > 0)
      gss_df <- gss_df %>%
      dplyr::mutate(cell_line = gss_cell_line, 
                    gene_silencing_series = gss_id, .before = sample_id)
  }
  
  return(tg_df)
}


#' Gets the experiment type (e.g. control/case)
#'
#' @param related_dataset Data.frame with the information about the experiment
#'   as provided by the ENCODE API.
#' @param target_gene Target gene of study.
#' @param gss_cell_line Cell line of the gss.
#'
#' @return The experiment type (e.g. control/case)
#' @export
GetExperimentType <- function(related_dataset,
                              target_gene,
                              gss_cell_line){
  if(!"control_type" %in% names(related_dataset)){
    logger::log_warn(" No control experiment found for target gene ", target_gene, " and cell line ", gss_cell_line, ".")
    experiment_type = "case"
  }else{
    experiment_type <- ifelse(is.na(related_dataset$control_type), "case", "control")
  }
  
  return(experiment_type)
}

#' Gets the experiment additional information.
#'
#' @param related_dataset Data.frame with the information about the experiment
#'   as provided by the ENCODE API.
#'
#' @return A data.frame with the metadata about the laboratory, the assay, the
#'   cellosaurus, the nucleic acid type, the extraction method, the
#'   fragmentation method, the size selection method and the strand specificity.
#' @export
GetAdditionalInformation <- function(related_dataset){
  
  sample_lab <- related_dataset$lab.title
  
  if(is.null(sample_lab)){
    sample_lab <- related_dataset$lab$title
  }
  
  sample_assay <- related_dataset$assay_term_name
  
  if (sample_assay == "long read RNA-seq") {
    
    sample_cellosaurus <- related_dataset$biosample_ontology$dbxrefs %>% 
      unlist() %>% 
      unique()
    
    sample_nucleic_acid_type <- related_dataset$replicate$library$nucleic_acid_term_name %>% 
      unique()
    
    sample_extraction_method <- NULL
    
    sample_fragmentation_method <- NULL
    
    sample_size_selection_method <- related_dataset$file_size
    
    sample_strand_specificity <-  related_dataset$replicate$library$strand_specificity
    
  } else {
    
    sample_cellosaurus <- related_dataset$biosample_ontology.dbxrefs %>% 
      unlist() %>% 
      unique()
    if(is.null(sample_cellosaurus)){
      sample_cellosaurus <- related_dataset$biosample_ontology$dbxrefs %>% 
        unlist() %>% 
        unique()
    }
    
    sample_nucleic_acid_type <- related_dataset$replicates %>% 
      data.frame() %>% 
      tibble::as_tibble() %>% 
      dplyr::pull(library.nucleic_acid_term_name) %>% 
      unique()
    sample_extraction_method <- related_dataset$replicates %>% 
      data.frame() %>% 
      tibble::as_tibble() %>% 
      dplyr::pull(library.extraction_method) %>% 
      unique()
    sample_fragmentation_method <- related_dataset$replicates %>% 
      data.frame() %>% 
      tibble::as_tibble() %>% 
      dplyr::pull(library.fragmentation_methods) %>% 
      unlist() %>% 
      unique()
    sample_size_selection_method <- related_dataset$replicates %>% 
      data.frame() %>% 
      tibble::as_tibble() %>% 
      dplyr::pull(library.library_size_selection_method) %>% 
      unique()
    sample_strand_specificity <- related_dataset$replicates %>% 
      data.frame() %>% 
      tibble::as_tibble() %>% 
      dplyr::pull(library.strand_specificity) %>% 
      unique()
  }
  
  
  tibble::tibble(lab = sample_lab,
                 assay = sample_assay,
                 cellosaurus = sample_cellosaurus,
                 nucleic_acid_type = sample_nucleic_acid_type,
                 extraction_method = sample_extraction_method,
                 fragmentation_method = sample_fragmentation_method,
                 size_selection_method = sample_size_selection_method,
                 strand_specificity = sample_strand_specificity) %>%
    return()
}



#' Gets the experiment's RIN information
#'
#' @param replicates Data.frame containing information about the isogenic
#'   replicates within the experiment.
#'
#' @return Data.frame with the RIN information for each experiment's isogenic
#'   replicates.
#' @export
GetRin <- function(replicates) {
  if(!"library.rna_integrity_number" %in% names(replicates)){
    logger::log_warn("\t\t\t No RIN found. Defaulted to NA.")
    rin_info <- replicates %>% 
      dplyr::select(biological_replicate_number) %>% 
      dplyr::mutate(rin = NA)
  }else{
    rin_info <- replicates %>% 
      dplyr::select(biological_replicate_number, library.rna_integrity_number)
  }
  names(rin_info) <- c("bio_rep", "rin")
  
  return(rin_info)
}

#' Gets the experiment's sample read depth
#'
#' @param analyses Data.frame containing information about the different
#'   analyses executed within the experiment.
#' @param valid_genome_annotation (Optional) Required gene annotation version of
#'   the sample to extract its metadata. Defaults to "V29".
#'
#' @return Data.frame containing the read depth for each isogenic replicate.
#' @export
GetReadDepth <- function(analyses,
                         valid_genome_annotation = "V29") {
  read_depth_info <- analyses %>% 
    dplyr::filter(genome_annotation == valid_genome_annotation) %>%
    dplyr::pull(`quality_metrics_report.Read depth`) %>%
    magrittr::extract2(1) %>% 
    dplyr::select(biological_replicates, metric) %>%
    tidyr::unnest(biological_replicates)
  names(read_depth_info) <- c("bio_rep", "read_depth")
  
  return(read_depth_info)
}

#' Gets the sample donor information
#'
#' @param replicates Data.frame containing information about the isogenic
#'   replicates within the experiment.
#'
#' @return Data.frame containing information about the sample donor.
#' @export
GetDonorInfo <- function(replicates) {
  donor_info <- replicates %>% 
    dplyr::select(biological_replicate_number, library.biosample.sex, library.biosample.age, library.biosample.life_stage)
  names(donor_info) <- c("bio_rep", "sex", "age", "life_stage")
  
  return(donor_info)
}

#' Get documents related to the experiment
#'
#' @param replicates Data.frame containing information about the isogenic
#'   replicates within the experiment.
#' @param target_gene Target gene of study.
#' @param experiment_type Experiment type (e.g. control/case).
#'
#' @return Data.frame containing the information about the biosample and characterization documents.
#' @export
GetDocumentFiles <- function(replicates, 
                             target_gene, 
                             experiment_type){
  if(experiment_type == "case"){
    documents_info <- foreach(iter_replicate = seq(nrow(replicates)), .combine = "rbind") %do%{
      replicate <- replicates[iter_replicate, ]
      
      bio_rep <- replicate$biological_replicate_number
      accession <- replicate$library.accession
      documents <- replicate$library.biosample.documents %>% unlist()
      if(is.null(documents)) documents <- NA
      
      aliases <- replicate$library.biosample.aliases %>% unlist()
      if(length(aliases) > 1){
        aliases <- aliases[which(!grepl(":BG[a-zA-Z]LV", aliases))]
      }
      aliases <- stringr::str_split(aliases, ":|,", simplify = T)[, 2]
      
      tibble::tibble(bio_rep = bio_rep, 
                     biosample = accession, 
                     document = documents, 
                     biosample_alias = aliases)
    }
  }else{
    documents_info <- tibble(bio_rep = c(1, 2))
  }
  
  return(documents_info)
}

#' Gets the gene quantification files from the experiment
#'
#' @param files Data.frame containing all the files found for the experiment.
#' @param valid_file_format (Optional) Required output file format of the
#'   sample. Defaults to "TSV".
#' @param valid_output_type (Optional) Required output type of the sample.
#'   Defaults to "gene quantifications".
#' @param valid_genome_annotation (Optional) Required gene annotation version of
#'   the sample to extract its metadata. Defaults to "V29".
#'
#' @return Data.frame containing the information of the gene quantifications
#'   files found within the experiment.
#' @export
GetGeneQuantificationFiles <- function(files,
                                       valid_file_format = "tsv",
                                       valid_output_type = "gene quantifications",
                                       valid_genome_annotation = "V29") {
  tsv_files_info <- files %>% 
    dplyr::filter(file_format == valid_file_format) %>%
    dplyr::filter(output_type == valid_output_type) %>%
    dplyr::filter(genome_annotation == valid_genome_annotation) %>%
    dplyr::select(accession, biological_replicates, file_format, output_type, genome_annotation, technical_replicates) %>%
    tidyr::unnest(c(biological_replicates, technical_replicates)) %>%
    dplyr::rename("bio_rep" = "biological_replicates",
                  "tech_rep" = "technical_replicates",
                  "gene_quantification_id" = "accession")
  
  return(tsv_files_info)
}


#' Extract metadata from the case/control experiments
#'
#' @param summary_tg_df Local target gene information from the summary
#'   data.frame generated with [GenerateSummary()].
#' @param valid_nucleic_acid_type (Optional) Required nucleic acid type of the
#'   experiment. Defaults to "polyadenylated mRNA".
#' @param valid_file_format (Optional) Required output file format of the
#'   sample. Defaults to "BAM".
#' @param valid_output_type (Optional) Required output type of the sample.
#'   Defaults to "alignments".
#' @param valid_genome_annotation (Optional) Required gene annotation version of
#'   the sample to extract its metadata. Defaults to "V29".
#'
#' @return Data.frame containing the metadata for case and control experiments.
#' @export
LoopExperiments <- function(summary_tg_df,
                            valid_nucleic_acid_type,
                            valid_file_format,
                            valid_output_type,
                            valid_genome_annotation){
  ## Loop through every gene silencing series found for the target gene
  tg_df <- foreach(iter_exp = seq(nrow(summary_tg_df)), .combine = dplyr::bind_rows) %do%{
    
    # iter_exp <-1
    exp <- summary_tg_df[iter_exp, ]
    
    ## Experiment information
    exp_id <- exp$experiment_id
    gss_id <- exp$gene_silencing_series
    exp_cell_line <- exp$cell_line
    
    logger::log_info("\t Starting Experiment ", exp_id, " (cell line = ", exp_cell_line, "):")
    
    ## Get the API response for the experiment
    response_case <- GetUrlResponse(paste0("https://www.encodeproject.org/experiment/", exp_id, "?format=json"))
    control_id = response_case$possible_controls$accession
    
    if (is.null(control_id)) {
      control_id = response_case$replicate$experiment$possible_controls %>% unlist
    }
    
    if (!is.null(response_case)) {
      case_df <- ExtractMetadataExperiment(experiment = response_case, 
                                           experiment_type = "case", 
                                           valid_nucleic_acid_type, 
                                           valid_file_format,
                                           valid_output_type,
                                           valid_genome_annotation)
      control_df <- NULL
    }
    if (!is.null(control_id)) {
      
      response_control = GetUrlResponse(paste0("https://www.encodeproject.org/experiment/", control_id, "?format=json"))
      control_df <- ExtractMetadataExperiment(response_control, 
                                              "control", 
                                              valid_nucleic_acid_type, 
                                              valid_file_format,
                                              valid_output_type,
                                              valid_genome_annotation)
    }
    
    ## Combine cases and controls
    exp_df <- dplyr::bind_rows(case_df, control_df)
    
    ## Add the relevant information
    exp_df <- exp_df %>% 
      dplyr::mutate(cell_line = exp_cell_line, gene_silencing_series = gss_id, .before = sample_id)
  }
  
  return(tg_df)
}



#' Extracts the metadata for a particular experiment
#'
#' @param experiment Data.frame containing the information of the experiment as
#'   obtained from the ENCODE portal.
#' @param experiment_type Experiment type (e.g. control/case).
#' @param valid_nucleic_acid_type (Optional) Required nucleic acid type of the
#'   experiment. Defaults to "polyadenylated mRNA".
#' @param valid_file_format (Optional) Required output file format of the
#'   sample. Defaults to "BAM".
#' @param valid_output_type (Optional) Required output type of the sample.
#'   Defaults to "alignments".
#' @param valid_genome_annotation (Optional) Required gene annotation version of
#'   the sample to extract its metadata. Defaults to "V29".
#'
#' @return Data.frame containing the metadata for a particular experiment.
#' @export 
ExtractMetadataExperiment <- function(experiment, 
                                      experiment_type,
                                      valid_nucleic_acid_type,
                                      valid_file_format,
                                      valid_output_type,
                                      valid_genome_annotation){
  
  experiment_id = experiment$accession
  
  logger::log_info("\t\t Starting experiment ", experiment_id, " (type = ", experiment_type, ").")
  
  ## Requirements metadata
  experiment_additional_info <- GetAdditionalInformation(related_dataset = experiment)
  
  ## Main metadata
  experiment_sample_files <- GetSampleFiles(files = experiment$files[[1]], 
                                            valid_file_format, valid_output_type, valid_genome_annotation)
  
  if(experiment_additional_info$nucleic_acid_type != valid_nucleic_acid_type) 
    return(cbind(experiment_sample_files, experiment_additional_info))
  
  ## Other metadata
  experiment_rin <- GetRin(experiment$replicates)
  experiment_read_depth <- GetReadDepth(experiment$analyses, valid_genome_annotation)
  experiment_donor_info <- GetDonorInfo(experiment$replicates)
  experiment_documents <- GetDocumentFiles(experiment$replicates, gss$target_gene, experiment_type)
  experiment_doi <- experiment$doi
  experiment_gene_quantifications <- GetGeneQuantificationFiles(experiment$files)
  
  ## Combine all information
  experiment_combined <- experiment_sample_files %>% 
    dplyr::left_join(experiment_rin, by = "bio_rep") %>%
    dplyr::left_join(experiment_read_depth, by = "bio_rep") %>%
    dplyr::left_join(experiment_donor_info, by = "bio_rep") %>% 
    dplyr::left_join(experiment_gene_quantifications %>% 
                       dplyr::select(gene_quantification_id, bio_rep), 
                     by = "bio_rep") %>%
    dplyr::left_join(experiment_documents, by = "bio_rep") %>%
    dplyr::mutate(experiment_type = experiment_type, experiment_id = experiment_id, experiment_doi = experiment_doi, .before = "bio_rep") %>%
    dplyr::cbind(experiment_additional_info)
  
}


#' Gets the experiment's sample files
#'
#' @param files Data.frame containing all the files found for the experiment.
#' @param valid_file_format (Optional) Required output file format of the
#'   sample. Defaults to "BAM".
#' @param valid_output_type (Optional) Required output type of the sample.
#'   Defaults to "alignments".
#' @param valid_genome_annotation (Optional) Required gene annotation version of
#'   the sample to extract its metadata. Defaults to "V29".
#'
#' @return A data.frame with the metadata about the sample files, biological
#'   replicate, output type, etc.
#' @export
GetSampleFiles <- function(files,
                           valid_file_format = "bam",
                           valid_output_type = "alignments",
                           valid_genome_annotation = "V29") {
  
  sample_files_info <- files %>% 
    dplyr::filter(file_format == valid_file_format) %>%
    dplyr::filter(output_type %in% valid_output_type) %>%
    dplyr::filter(genome_annotation == valid_genome_annotation) %>%
    dplyr::select(accession, biological_replicates, file_format, output_type, genome_annotation, technical_replicates, mapped_run_type) %>%
    tidyr::unnest(c(biological_replicates, technical_replicates)) %>%
    dplyr::rename("bio_rep" = "biological_replicates",
                  "tech_rep" = "technical_replicates",
                  "sample_id" = "accession")
  
  return(sample_files_info)
}

#' Adds the functional category to the target genes
#'
#' Given the metadata dataframe, it adds a column to register the functional
#' category of the different target genes. A column is added for each category
#' being considered: Splicing regulation, Spliceosome, Exon junction complex or
#' Nononsense-mediated decay.
#'
#' @param metadata_df Data.frame where every row is an ENCODE sample and the
#'   columns contain their metadata.
#' @param input_Category Path to the file containing the functional category
#'   information obtained from Nostrand et. al. publication
#'   (https://www.nature.com/articles/s41586-020-2077-3).
#' @param input_NMD Path to the file containing the NMD functional category
#'   information as a list of genes from this category.
#' @param output_file Path to the output file where the metadata data.frame will
#'   be stored.
#'
#' @return Data.frame where every row is an ENCODE sample and the columns
#'   contain their metadata, including their functional category.
#' @export
AddTargetGeneCategory <- function(metadata_df,
                                  input_Category = "",
                                  input_NMD = "",
                                  output_file = ""){
  
  
  metadata_df <- readr::read_delim(file = metadata_df, delim = "\t", show_col_types = F)
  target_RBPs_metadata <- readxl::read_xlsx(path = input_Category, sheet = 1) %>% drop_na()
  
  
  metadata_df <- metadata_df %>% 
    dplyr::left_join(target_RBPs_metadata,# %>% 
                     #dplyr::select(name, `Splicing regulation`, Spliceosome, `Novel RBP`, `Exon Junction Complex`) %>%
                     #dplyr::rename("Splicing_regulation" = "Splicing regulation",
                     #               "Novel_RBP" = "Novel RBP",
                     #               "Exon_junction_complex" = "Exon Junction Complex"),
                     by = c("target_gene" = "name"))
  
  if(input_NMD != ""){
    
    NMD_list <- readr::read_delim(input_NMD, show_col_types = F, delim = "\n") %>% 
      dplyr::rename("Name" = "Name\t") %>% 
      dplyr::mutate(across(Name, str_replace, "\t", "")) %>% 
      dplyr::pull(Name)
    
    metadata_df <- metadata_df %>%
      dplyr::mutate(NMD = ifelse(target_gene %in% NMD_list, 1, 0))
  }
  
  if(output_file != ""){
    write.table(metadata_df, output_file, sep = "\t", row.names = F, quote = FALSE)
  }
  
  return(metadata_df)
}

## 28