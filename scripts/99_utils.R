#' Remove the junctions shorter than 25bp.
#'
#' @param input.SR.details Dataframe object with the relevant junctions.
#'
#' @return Junctions bigger than 25bp.
#' @export
RemoveShortJunctions <- function(input.SR.details) {
  logger::log_info("\t\t Removing junctions shorter than 25bp.")
  output_SR_details <- input.SR.details %>%
    dplyr::filter(width >= 25)
  
  return(output_SR_details)
  
}

#' Title
#' Removes all split reads with coordinates overlapping any of the regions published within the ENCODE blacklist 
#' @param GRdata Genome Ranges object with the split reads to analyse
#' @param blacklist.path Local path to the .bed file containing the ENCODE blacklist 
#'
#' @return The 'GRdata' object without split reads overlapping the ENCODE blacklist 
#' @export
#'
#' @examples
RemoveEncodeBlacklistRegions <- function(GRdata, blacklist.path) {
  
  
  if (!exists("encode_blacklist_hg38")) {
    encode_blacklist_hg38 <- rtracklayer::import(con = blacklist.path) %>% diffloop::rmchr()
  } else {
    print("'encode_blacklist_hg38' file already loaded!")
  }
  
  overlaped_junctions <- GenomicRanges::findOverlaps(query = encode_blacklist_hg38, 
                                                     subject = GRdata %>% diffloop::rmchr(),
                                                     type = "any",
                                                     ignore.strand = F)
  
  ## JuncID indexes to be removed: they overlap with a black region
  indexes <- S4Vectors::subjectHits(overlaped_junctions)
  
  if (length(indexes) > 0) {
    print(paste0(length(unique(indexes)), " junctions overlap with a ENCODE blacklist region"))
    GRdata <- GRdata[-indexes, ]
    print(paste0(length(GRdata), " junctions after removing overlaps with ENCODE BlackList regions!"))
  }else{
    print("No junctions overlapping with an ENCODE blacklist region")
  }
  
  ## Return tidied data
  return(GRdata)
}

#' Loads the reference genome into memory
#'
#' The different versions can be downloaded
#' \href{http://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/}{here}.
#'
#' @param gtf.path Path to the reference genome GTF file.
#'
#' @return Connection to the reference genome DB.
#' @export
LoadEdb <- function(gtf.path) {
  if (!exists("edb")) {
    logger::log_info("\t\t Loading the reference genome.")
    edb <<- ensembldb::ensDbFromGtf(gtf.path, outfile = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
    edb <<- ensembldb::EnsDb(x = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
  } else {
    logger::log_info("\t\t Variable 'edb' already loaded!")
  }
  return(edb)
}


#' Annotate and categorize every junction
#'
#' @param GRdata GRanges class object with the relevant junctions.
#' @param edb The connection to the reference genome DB.
#'
#' @return Annotated input by dasper.
#' @export 
AnnotateDasper <- function(GRdata, edb) {
  logger::log_info("\t\t Annotating using dasper::junction_annot().")
  GRdata <- dasper::junction_annot(GRdata,
                                   ref = edb,
                                   ref_cols = c("gene_id", "gene_name", "symbol", "tx_id"),
                                   ref_cols_to_merge = c("gene_id", "gene_name", "tx_id")
  )
  
  return(GRdata)
}


#' Remove the junctions with uncategorized classification.
#'
#' @param input.SR.details Dataframe object with the relevant junctions.
#'
#' @return Junctions categorized as 'annotated', 'novel_donor', 'novel_acceptor', 'novel combo' or 'novel exon skip' 
#' @export
RemoveUncategorizedJunctions <- function(input.SR.details) {
  
  logger::log_info("\t\t Removing split reads not classified as 'annotated', 'novel_donor', 'novel_acceptor', 'novel combo' or 'novel exon skip' ...")

  ## Subset columns
  input.SR.details <- input.SR.details %>% 
    dplyr::select(any_of(c("junID", "seqnames", "start", "end", "width", "strand",
                           "gene_id_junction", "in_ref", "type", "tx_id_junction",
                           "annotated", "left_motif", "right_motif", "n_projects")))
  
  ## Only use annotated introns, novel donor and novel acceptor junctions
  output_SR_details <- input.SR.details %>%
    filter(type %in% c("annotated", "novel_donor", "novel_acceptor", "novel_combo", "novel_exon_skip"))
  
  return(output_SR_details)
  
}


#' Remove the junctions from ambiguous genes.
#'
#' @param input_SR_details Dataframe object with the relevant junctions.
#'
#' @return Junctions assigned to only one gene.
#' @export
RemoveAmbiguousJunctions <- function(input_SR_details, database.folder) {
  
  logger::log_info("\t\t Removing junctions associated to more than one gene.")

  
  input_SR_details <- input_SR_details %>%
    as_tibble() %>%
    distinct(junID, .keep_all = T) %>% 
    rowwise() %>%
    mutate(ambiguous = ifelse(gene_id_junction %>% unlist() %>% length() > 1, T, F))
  
  ambiguous_introns <- input_SR_details %>% dplyr::filter(ambiguous == T)
  
  logger::log_info("\t\t Removing ", nrow(ambiguous_introns)," ambiguous junctions!")
  
  saveRDS(object = ambiguous_introns, file = file.path(database.folder, "all_ambiguous_jxn.rds"))
  
  return(input_SR_details %>%
           filter(ambiguous == F))
  
}


#' Title
#' Calculates the cummulative number of reads and number of samples for a given junction
#' @param split.read.counts dataframe of reads counts per junction. The first column correspond to the junction ID and the rest of the columns
#' correspond to the samples. Each value represents the number of reads for a given junction in a given sample.
#' @param samples Vector of samples to consider
#' @param junIDs List of junction IDs to consider
#'
#' @return
#' @export
#'
#' @examples
GenerateCoverage <- function(split.read.counts, samples, junIDs) {
  
  # stopifnot(
  #   "Still there are split reads with less than 2 supportive reads" =
  #     split.read.counts %>% 
  #     mutate(sumCounts = rowSums(select(., !contains("junID")))) %>%
  #     filter(sumCounts <= 2) %>% 
  #     nrow() == 0
  # )
  
  split_read_counts_intron <- split.read.counts %>%
    dplyr::filter(junID %in% junIDs) %>%
    dplyr::select(junID, all_of(samples %>% as.character())) 
  
  split_read_counts_intron[,"n_individuals"] <- (matrixStats::rowCounts(split_read_counts_intron[, -c(1), drop=FALSE] > 0, na.rm = T)) 
  split_read_counts_intron <- split_read_counts_intron %>% as.data.frame()
  
  split_read_counts_intron[,"sum_counts"] <- Matrix::rowSums(split_read_counts_intron[,-c(split_read_counts_intron %>% ncol(),1), drop=FALSE], na.rm = T)
  split_read_counts_intron <- split_read_counts_intron %>% as.data.frame()
  
  split_read_counts_intron <- split_read_counts_intron[, c(1,(split_read_counts_intron %>% ncol() - 1),(split_read_counts_intron %>% ncol())), drop=FALSE]
  
  
  if (any(split_read_counts_intron[, "n_individuals"] < 1)) {
    print("Error: some ref junctions do not present any read across any of the samples.")
    break;
  }
  
  split_read_counts_intron %>% return()
}


#' Title
#' Splits a string with the formatting: "chrX:start-end:strand" into separated elements 'seqnames', 'start', 'end', 'strand'
#' @param coordinates 
#'
#' @return
#' @export
#'
#' @examples
GetGenomicCoordinates <- function(coordinates) {
  
  map_df(coordinates, function(coordinate) {
    # coordinate <- df_gene_splicing$novel_coordinates[1]
    chr_junc <- coordinate %>%
      str_sub(start = 1,
              end = str_locate_all(string = coordinate, pattern = ":")[[1]][1,2]-1)
    start_junc <- coordinate %>%
      str_sub(start = str_locate_all(string = coordinate, pattern = ":")[[1]][1,2]+1,
              end = str_locate_all(string = coordinate, pattern = "-")[[1]][1,2]-1)
    end_junc <- coordinate %>%
      str_sub(start = str_locate_all(string = coordinate, pattern = "-")[[1]][1,2]+1,
              end = str_locate_all(string = coordinate, pattern = ":")[[1]][2,2]-1)
    strand_junc <- coordinate %>%
      str_sub(start = str_locate_all(string = coordinate, pattern = ":")[[1]][2,2]+1,
              end = coordinate %>% stringr::str_count())
    
    data.frame(ID = coordinate,
               seqnames = chr_junc,
               start = start_junc %>% as.integer(),
               end = end_junc %>% as.integer(),
               strand = strand_junc) %>%
      return()
    
  })
}


#' Title
#' Function to calculate the TPM value per gene
#' @param rse RangedSummarizedExperiment-class object
#' @param ref.tidy Reference transcriptome
#'
#' @return
#' @export
#'
#' @examples
GenerateTPM <- function(rse, ref.tidy) {
  
  
  # Remove anything after . in ensembl id
  rownames(rse) <- rownames(rse) %>% 
    str_remove("\\..*")
  
  
  # Convert to tpm, which is calculated by:
  # 1. Divide the read counts by the length of each gene in kilobases (i.e. RPK)
  # 2. Count up all the RPK values in a sample and divide this number by 1,000,000.
  # 3. Divide the RPK values by the “per million” scaling factor.
  
  message("Calculating RPK values ...")
  srp_rpk <- 
    rse %>% 
    SummarizedExperiment::assay() %>%
    as_tibble(rownames = "gene") %>% 
    tidyr::pivot_longer( ## equivalent to gather
      cols = -c("gene"),
      names_to = "recount_id",
      values_to = "counts"
    ) %>% 
    dplyr::inner_join(
      ref.tidy %>% 
        as_tibble() %>% 
        dplyr::select(gene_id, width),
      by = c("gene" = "gene_id")
    ) %>% # 1. Divide the read counts by the length of each gene in kilobases (i.e. RPK)
    dplyr::mutate(
      rpk = counts/width
    ) 
  srp_rpk %>% head()
  
  message("Transforming RPK into TPM values ...")
  tpm <- 
    srp_rpk %>% 
    # 2. Count up all the RPK values in a sample and divide this number by 1,000,000.
    dplyr::inner_join(
      srp_rpk %>% 
        dplyr::group_by(recount_id) %>% 
        dplyr::summarise(
          scaling_factor = sum(rpk)/1e6
        ),
      by = "recount_id"
    ) %>% # 3. Divide the RPK values by the “per million” scaling factor.
    dplyr::mutate(
      tpm = rpk/scaling_factor
    ) %>%
    dplyr::select(gene, recount_id, tpm) %>% 
    spread(key = recount_id, value = tpm)
  
  
  return(tpm)
}



#' Title
#' Converts a string dataframe into a numeric dataframe
#' @param sample.metadata 
#' @param samples 
#'
#' @return
#' @export
#'
#' @examples
TidySampleMetadata <- function(sample.metadata, samples) {
  
  
  
  age_numeric <- as.numeric(factor(as.matrix(sample.metadata$gtex.age))) 
  sample.metadata$gtex.age <- age_numeric
  
  sample.metadata <- sample.metadata %>%
    as.data.table() %>%
    tibble::column_to_rownames(var = "external_id")
  
  covariates <- c("gtex.age", "gtex.smcenter", "gtex.smtsd",
                  "gtex.smgebtch", "gtex.smgebtchd",
                  "gtex.smnabtch", "gtex.smnabtchd", "gtex.smnabtcht",
                  "gtex.dthhrdy", "gtex.sex", "gtex.smrin")
  
  
  smcenter <- as.numeric(factor(as.matrix(sample.metadata$gtex.smcenter))) 
  sample.metadata$gtex.smcenter <- smcenter
  
  gtex.smgebtch <- as.numeric(factor(as.matrix(sample.metadata$gtex.smgebtch))) 
  sample.metadata$gtex.smgebtch <- gtex.smgebtch
  
  gtex.smgebtchd <- as.numeric(factor(as.matrix(sample.metadata$gtex.smgebtchd))) 
  sample.metadata$gtex.smgebtchd <- gtex.smgebtchd
  
  gtex.smnabtch <- as.numeric(factor(as.matrix(sample.metadata$gtex.smnabtch))) 
  sample.metadata$gtex.smnabtch <- gtex.smnabtch
  
  gtex.smnabtchd <- as.numeric(factor(as.matrix(sample.metadata$gtex.smnabtchd))) 
  sample.metadata$gtex.smnabtchd <- gtex.smnabtchd
  
  gtex.smnabtcht <- as.numeric(factor(as.matrix(sample.metadata$gtex.smnabtcht))) 
  sample.metadata$gtex.smnabtcht <- gtex.smnabtcht
  
  gtex.smtsd <- as.numeric(factor(as.matrix(sample.metadata$gtex.smtsd))) 
  sample.metadata$gtex.smtsd <- gtex.smtsd
  
  
  # 8. Return covariates ---------------------------
  
  return(t(sample.metadata %>%
             dplyr::select(all_of(covariates))))
}


#' Title
#' Calculates the mode value from a vector of numbers
#' @param vector 
#'
#' @return
#' @export
#'
#' @examples
GetMode <- function(vector) {
  uniqv <- unique(vector)
  uniqv[which.max(tabulate(match(vector, uniqv)))]
}