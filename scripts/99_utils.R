#' Remove the junctions shorter than 25bp.
#'
#' @param input.SR.details Dataframe object with the relevant junctions.
#'
#' @return Junctions bigger than 25bp.
#' @export
RemoveShortJunctions <- function(input.SR.details) {
  logger::log_info("Removing junctions shorter than 25bp...")
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
  
  
  if (file.exists(blacklist.path)) {
    encode_blacklist_hg38 <- rtracklayer::import(con = blacklist.path) 
  } else {
    stop(blacklist.path, " does not exist!")
  }
  
  overlaped_junctions <- GenomicRanges::findOverlaps(query = encode_blacklist_hg38, subject = GRdata,
                                                     type = "any", ignore.strand = F)
  
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

#' Collapse gene IDs that belong to the same isomatch gene group
#'
#' The URA GTF carries, per transcript, the gene_id of the source record that was selected for it
#' (Ensembl > GENCODE > RefSeq > CHESS). Inside one isomatch gene group (ISOM_GENE) transcripts can
#' therefore carry gene_ids from different sources, e.g. an Ensembl gene and a CHESS-only gene at the
#' same locus. A junction at that locus would then be annotated to two genes and treated as ambiguous.
#'
#' Rule: a gene is mapped to the gene of the highest-priority source in its group, but only if that
#' gene is unique. If the group has several top-priority genes (real overlapping genes), or the gene is
#' itself top-priority in any of its groups, the gene is kept as it is.
#'
#' @param tx Tibble with one row per URA transcript: gene_id, ISOM_GENE, rank (1 = highest priority).
#'
#' @return Named character vector: original gene_id -> collapsed gene_id.
#' @export
CollapseGeneIds <- function(tx) {
  
  gene_rank <- tx %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(gene_rank = min(rank), .groups = "drop")
  
  group_genes <- tx %>%
    dplyr::distinct(gene_id, ISOM_GENE) %>%
    dplyr::left_join(gene_rank, by = "gene_id") %>%
    dplyr::group_by(ISOM_GENE) %>%
    dplyr::mutate(best_rank = min(gene_rank),
                  n_best = sum(gene_rank == best_rank),
                  best_gene = gene_id[gene_rank == best_rank][1],
                  is_best = gene_rank == best_rank) %>%
    dplyr::ungroup()
  
  per_gene <- group_genes %>%
    dplyr::mutate(blocked = is_best | n_best > 1,
                  target = ifelse(!is_best & n_best == 1, best_gene, NA_character_)) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(blocked = any(blocked),
                     targets = list(unique(stats::na.omit(target))),
                     .groups = "drop")
  
  gene_map <- stats::setNames(per_gene$gene_id, per_gene$gene_id)
  to_collapse <- !per_gene$blocked & lengths(per_gene$targets) == 1
  gene_map[per_gene$gene_id[to_collapse]] <- unlist(per_gene$targets[to_collapse])
  
  ## a target can itself be collapsed through another group: follow the chain (rank strictly decreases)
  for (i in seq_len(4)) {
    gene_map <- stats::setNames(unname(gene_map[gene_map]), names(gene_map))
  }
  
  return(gene_map)
}


#' Components of an undirected graph given as an edge list
#' @keywords internal
.ComponentsFromEdges <- function(nodes, from, to) {
  parent <- seq_along(nodes)
  a <- match(from, nodes)
  b <- match(to, nodes)
  find <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  for (i in seq_along(a)) {
    ra <- find(a[i]); rb <- find(b[i])
    if (ra != rb) parent[rb] <- ra
  }
  vapply(seq_along(nodes), find, integer(1))
}


#' Merge gene IDs that are the same gene annotated by different sources
#'
#' Second stage after CollapseGeneIds(). Two (collapsed) genes are merged when
#' (1) they are annotated for the same multi-exon merged transcript (ISOM_SRC lists the records of
#'     every source that has the identical intron chain, so their gene IDs are the same gene), or
#' (2) they have the same gene_name, are on the same chromosome and strand and their extents overlap.
#'     Overlap is required because many small RNA genes share a name (Y_RNA, U6, ...) across unrelated loci.
#' Every merged group is mapped to the gene of the highest-priority source (ties: gene_id order).
#'
#' @param tx Tibble with one row per URA transcript: gene_id, gene_name, seqnames, start, end, strand,
#' rank (1 = highest priority), ISOM_SRC, ISOM_EXONS.
#' @param gene_map Named vector from CollapseGeneIds(): gene_id -> collapsed gene_id.
#'
#' @return Named character vector: original gene_id -> merged gene_id.
#' @export
MergeEquivalentGenes <- function(tx, gene_map) {
  
  ## Node table: one row per collapsed gene
  nodes <- tx %>%
    dplyr::mutate(node = unname(gene_map[gene_id])) %>%
    dplyr::group_by(node) %>%
    dplyr::summarise(rank = min(rank),
                     seqnames = dplyr::first(seqnames), strand = dplyr::first(strand),
                     start = min(start), end = max(end),
                     gene_name = dplyr::first(gene_name[gene_id == node]),
                     .groups = "drop")
  
  ## 1. Cross-references: genes of the different records of one multi-exon merged transcript
  records <- tx %>%
    dplyr::filter(as.integer(ISOM_EXONS) > 1, grepl("|", ISOM_SRC, fixed = TRUE)) %>%
    dplyr::select(transcript_id, ISOM_SRC) %>%
    tidyr::separate_rows(ISOM_SRC, sep = "\\|") %>%
    tidyr::separate(ISOM_SRC, into = c("sample", "tx", "gene"), sep = ":", extra = "drop", fill = "right") %>%
    dplyr::mutate(gene = ifelse(sample == "S3", sub("\\.[0-9]+", "", gene), gene),
                  node = unname(gene_map[gene])) %>%
    dplyr::filter(!is.na(node))
  
  edges_xref <- records %>%
    dplyr::distinct(transcript_id, node) %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::mutate(anchor = dplyr::first(node)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(node != anchor) %>%
    dplyr::distinct(from = anchor, to = node)
  
  ## 2. Same symbol, same chromosome and strand, overlapping extents
  edges_symbol <- nodes %>%
    dplyr::filter(!is.na(gene_name), strand %in% c("+", "-")) %>%
    dplyr::group_by(gene_name, seqnames, strand) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::arrange(start, .by_group = TRUE) %>%
    dplyr::mutate(running_end = dplyr::lag(cummax(end)),
                  new_cluster = is.na(running_end) | start > running_end,
                  cluster = cumsum(new_cluster)) %>%
    dplyr::group_by(gene_name, seqnames, strand, cluster) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::mutate(anchor = dplyr::first(node)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(node != anchor) %>%
    dplyr::distinct(from = anchor, to = node)
  
  edges <- dplyr::bind_rows(edges_xref, edges_symbol)
  nodes$component <- .ComponentsFromEdges(nodes$node, edges$from, edges$to)
  
  ## Representative of each component: highest-priority source, then gene_id
  reps <- nodes %>%
    dplyr::group_by(component) %>%
    dplyr::arrange(rank, node, .by_group = TRUE) %>%
    dplyr::mutate(representative = dplyr::first(node), size = dplyr::n()) %>%
    dplyr::ungroup()
  
  node_map <- stats::setNames(reps$representative, reps$node)
  merged_map <- stats::setNames(unname(node_map[gene_map]), names(gene_map))
  
  largest <- reps %>% dplyr::distinct(representative, size) %>% dplyr::arrange(dplyr::desc(size)) %>% utils::head(5)
  logger::log_info(sum(reps$node != reps$representative), " further gene IDs merged (",
                   nrow(edges_xref), " cross-reference links, ", nrow(edges_symbol), " same-symbol links). ",
                   "Largest merged groups: ", paste0(largest$representative, " (", largest$size, ")", collapse = ", "))
  
  return(merged_map)
}


#' Loads the URA reference annotation into memory
#'
#' Builds a TxDb from the enriched URA GTF (script.sh + build_biotype_lookup.sh), with the gene_id of
#' every transcript collapsed by CollapseGeneIds(), and a table with one row per gene.
#'
#' @param gtf.path Path to the enriched URA GTF.
#' @param merge.equivalent.genes If TRUE (default), genes that are the same gene annotated by different sources
#' are merged (see MergeEquivalentGenes()) after the isomatch-group collapse.
#'
#' @return List with `txdb` (TxDb used by dasper, collapsed gene IDs, Ensembl-style seqlevels),
#' `genes` (tibble: gene_id, gene_name, gene_biotype, rank; uncollapsed IDs, i.e. the gene table) and
#' `gene_map` (named vector: gene_id -> collapsed gene_id).
#' The result is cached in the global variable `ura_reference` and reloaded automatically if the cached
#' TxDb connection is no longer valid or `merge.equivalent.genes` changed; rm(ura_reference) forces a reload.
#' @export
LoadUraReference <- function(gtf.path, merge.equivalent.genes = TRUE) {
  ## A TxDb keeps an SQLite connection that is dead after a workspace restore: reuse the cached
  ## reference only if its connection is valid and it was built with the same merge setting.
  cached_ok <- exists("ura_reference", envir = .GlobalEnv) &&
    isTRUE(tryCatch(DBI::dbIsValid(AnnotationDbi::dbconn(get("ura_reference", envir = .GlobalEnv)$txdb)),
                    error = function(e) FALSE)) &&
    identical(get("ura_reference", envir = .GlobalEnv)$merge.equivalent.genes, merge.equivalent.genes)
  
  if (!cached_ok) {
    logger::log_info("Loading the reference annotation.")
    
    gr <- rtracklayer::import(gtf.path, format = "gtf",
                              feature.type = c("transcript", "exon"),
                              colnames = c("type", "gene_id", "transcript_id", "gene_name",
                                           "gene_biotype", "ISOM_GENE", "ISOM_ID_SRC", "ISOM_SRC", "ISOM_EXONS"))
    
    source_rank <- c(S2 = 1L, S3 = 2L, S4 = 3L, S1 = 4L)
    tx <- gr[gr$type == "transcript"] %>%
      tibble::as_tibble() %>%
      dplyr::mutate(seqnames = as.character(seqnames), strand = as.character(strand))
    stopifnot("The GTF has no ISOM_GENE / ISOM_ID_SRC attributes: was it produced by the current script.sh?" =
                all(c("ISOM_GENE", "ISOM_ID_SRC") %in% colnames(tx)),
              "Unknown ISOM_ID_SRC value" = all(tx$ISOM_ID_SRC %in% names(source_rank)))
    tx$rank <- unname(source_rank[tx$ISOM_ID_SRC])
    
    genes <- tx %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(gene_name = dplyr::first(gene_name),
                       gene_biotype = dplyr::first(gene_biotype),
                       rank = min(rank), .groups = "drop")
    
    ## collapse gene IDs inside the reference seen by dasper
    gene_map <- CollapseGeneIds(tx)
    logger::log_info(sum(gene_map != names(gene_map)), " of ", length(gene_map),
                     " gene IDs collapsed into a higher-priority gene of the same isomatch group.")
    if (merge.equivalent.genes) {
      gene_map <- MergeEquivalentGenes(tx, gene_map)
    }
    gr$gene_id <- unname(gene_map[gr$gene_id])
    stopifnot("gene_id missing from the collapse map" = !anyNA(gr$gene_id))
    
    ## the rest of the pipeline uses Ensembl-style chromosome names (1, 2, ..., X, Y, MT)
    GenomeInfoDb::seqlevelsStyle(gr) <- "Ensembl"
    txdb <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr))
    
    ura_reference <<- list(txdb = txdb, genes = genes, gene_map = gene_map,
                           merge.equivalent.genes = merge.equivalent.genes)
  } else {
    logger::log_info("Variable 'ura_reference' already loaded!")
  }
  return(ura_reference)
}


#' Annotate and categorize every junction
#'
#' @param GRdata GRanges class object with the relevant junctions.
#' @param ref The reference returned by LoadUraReference().
#'
#' @return Annotated input by dasper, with the same columns as before:
#' gene_id_junction, gene_name_junction and tx_id_junction (transcript_id, version included).
#' @export
AnnotateDasper <- function(GRdata, ref) {
  logger::log_info("Annotating using dasper::junction_annot().")
  GRdata <- dasper::junction_annot(GRdata,
                                   ref = ref$txdb,
                                   ref_cols = c("gene_id", "tx_name"),
                                   ref_cols_to_merge = c("gene_id", "tx_name")
  )
  
  ## downstream code expects tx_id_junction (a TxDb has no EnsDb-style tx_id)
  S4Vectors::mcols(GRdata)$tx_id_junction <- S4Vectors::mcols(GRdata)$tx_name_junction
  S4Vectors::mcols(GRdata)$tx_name_junction <- NULL
  
  ## gene names are not part of a TxDb: add them from the URA gene table
  gene_names <- stats::setNames(ref$genes$gene_name, ref$genes$gene_id)
  gene_name_junction <- lapply(GRdata$gene_id_junction, function(g) unname(gene_names[g]))
  S4Vectors::mcols(GRdata)$gene_name_junction <- IRanges::CharacterList(gene_name_junction)
  
  return(GRdata)
}


#' Remove the junctions with uncategorized classification.
#'
#' @param input.SR.details Dataframe object with the relevant junctions.
#'
#' @return Junctions categorized as 'annotated', 'novel_donor', 'novel_acceptor', 'novel combo' or 'novel exon skip' 
#' @export
RemoveUncategorizedJunctions <- function(input.SR.details) {
  
  logger::log_info("Removing split reads not classified as 'annotated', 'novel_donor', 'novel_acceptor', 'novel combo' or 'novel exon skip' ...")
  
  ## Subset columns
  input.SR.details <- input.SR.details %>% 
    dplyr::select(any_of(c("junID", "seqnames", "start", "end", "width", "strand",
                           "gene_id_junction", "in_ref", "type", "tx_id_junction",
                           "annotated", "left_motif", "right_motif", "n_projects")))
  
  ## Only use annotated introns, novel donor and novel acceptor junctions
  output_SR_details <- input.SR.details %>%
    filter(!(type %in% c("unannotated"))) #"annotated", "novel_donor", "novel_acceptor", "novel_combo", "novel_exon_skip"
  
  return(output_SR_details)
  
}


#' Remove the junctions from ambiguous genes.
#'
#' @param input_SR_details Dataframe object with the relevant junctions.
#'
#' @return Datagrame object of junctions with an additional column indicating whether the jxn presents ambiguity or not 
#' ie ambiguity = the jxn has been annotated to more than one single gene
#' @export
RemoveAmbiguousJunctions <- function(input.SR.details, database.folder) {
  
  logger::log_info("Removing junctions associated to more than one gene...")
  
  input.SR.details <- input.SR.details %>%
    as_tibble() %>%
    distinct(junID, .keep_all = T) %>% 
    rowwise() %>%
    mutate(ambiguous = ifelse(gene_id %>% unlist() %>% length() > 1, T, F)) %>%
    ungroup()
  
  ambiguous_introns <- input.SR.details %>% dplyr::filter(ambiguous == TRUE)
  
  logger::log_info("There are ", nrow(ambiguous_introns)," ambiguous junctions!")
  
  saveRDS(object = ambiguous_introns, file = file.path(database.folder, "all_ambiguous_jxn.rds"))
  
  return(input.SR.details)
  
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



# gtf.path = "/rds/project/rds-2Hstq49W5EY/sruiz/reference/ensembl/Homo_sapiens.GRCh38.116.chr.gtf"
# utr.introns.path = "/rds/user/sg2173/hpc-work/moved_from_home_dir/PROJECTS/srRNAseq/recount3-database-project/dependencies/UTR_introns_hg38.116.rds"
ExtractUTRintronsFromReference <- function(gtf.path, utr.introns.path) {
  
  gtf_hg38 <- rtracklayer::import(file.path(gtf.path))
  gtf_hg38_tb <- gtf_hg38 %>% as_tibble()
  
  gtf_hg38_tb_utrs <- gtf_hg38_tb %>% filter(type %in% c("five_prime_utr", "three_prime_utr"))
  
  gtf_hg38_utr_introns_gr <- gtf_hg38_tb_utrs %>% 
    mutate(utr_type = type) %>% 
    ggtranscript::to_intron(., group_var = "transcript_id") %>% 
    dplyr::select(seqnames, start, end, strand, width, type, utr_type, gene_id, gene_name, gene_biotype, transcript_id, transcript_name, transcript_biotype) %>% 
    GRanges()
  
  gtf_hg38_utr_introns_gr$junID <- as.character(gtf_hg38_utr_introns_gr)
  
  # Evaluate using PSEN1 as example
  # gtf_hg38_tb %>% filter(type %in% c("five_prime_utr", "three_prime_utr"), gene_name == "PSEN1", transcript_id == "ENST00000324501") %>% as.data.frame()
  
  gtf_hg38_utr_introns_gr %>%
    as_tibble %>% 
    mutate(start = start + 1 ,end = end - 1) %>% ## Convert to intron coordinates
    GRanges %>%
    saveRDS(., file.path(utr.introns.path))
  
  gtf_hg38_utr_introns_gr %>% return()
}


#' Introns lying between two segments of the same UTR of a transcript
#'
#' @param segments Tibble with seqnames, strand, transcript_id, utr_type, seg_start, seg_end:
#' the UTR segments of each transcript (one per exon that has a UTR part).
#' @return Tibble with seqnames, start, end, strand, utr_type, transcript_id (intron coordinates).
#' @keywords internal
.UtrIntronsFromSegments <- function(segments) {
  segments %>%
    dplyr::group_by(transcript_id, utr_type) %>%
    dplyr::arrange(seg_start, .by_group = TRUE) %>%
    dplyr::mutate(next_start = dplyr::lead(seg_start)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(next_start)) %>%
    dplyr::transmute(seqnames, start = seg_end + 1L, end = next_start - 1L, strand, utr_type, transcript_id)
}


#' UTR segments of every coding transcript, derived from its exons and CDS
#'
#' 5' UTR = the part of the exons upstream of the first coding base (in transcription direction),
#' 3' UTR = the part downstream of the last coding base. The coding region is the CDS plus the stop codon
#' (Ensembl, GENCODE and RefSeq annotate it separately from the CDS; CHESS has no stop_codon features).
#'
#' @param exons,cds Tibbles with seqnames, start, end, strand, transcript_id (cds: CDS and stop_codon features).
#' @keywords internal
.UtrSegmentsFromExonsCds <- function(exons, cds) {
  cds_range <- cds %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::summarise(cds_min = min(start), cds_max = max(end), .groups = "drop")
  
  seg <- exons %>%
    dplyr::filter(strand %in% c("+", "-")) %>%
    dplyr::inner_join(cds_range, by = "transcript_id")
  
  lower <- seg %>%
    dplyr::filter(start < cds_min) %>%
    dplyr::mutate(seg_start = start,
                  seg_end = pmin(end, cds_min - 1L),
                  utr_type = ifelse(strand == "+", "five_prime_utr", "three_prime_utr"))
  upper <- seg %>%
    dplyr::filter(end > cds_max) %>%
    dplyr::mutate(seg_start = pmax(start, cds_max + 1L),
                  seg_end = end,
                  utr_type = ifelse(strand == "+", "three_prime_utr", "five_prime_utr"))
  
  dplyr::bind_rows(lower, upper) %>%
    dplyr::select(seqnames, strand, transcript_id, utr_type, seg_start, seg_end)
}


#' Rename RefSeq chromosome accessions (NC_0000NN.x, NC_012920.x) to Ensembl-style names and drop the rest
#' @keywords internal
.RefseqToEnsemblSeqlevels <- function(gr) {
  lv <- GenomeInfoDb::seqlevels(gr)
  num <- sub("^NC_0000(\\d{2})\\.\\d+$", "\\1", lv)
  is_nc <- grepl("^NC_0000\\d{2}\\.\\d+$", lv)
  new <- rep(NA_character_, length(lv))
  new[is_nc] <- ifelse(num[is_nc] == "23", "X", ifelse(num[is_nc] == "24", "Y", as.character(as.integer(num[is_nc]))))
  new[grepl("^NC_012920\\.\\d+$", lv)] <- "MT"
  keep <- !is.na(new)
  gr <- GenomeInfoDb::keepSeqlevels(gr, lv[keep], pruning.mode = "coarse")
  GenomeInfoDb::renameSeqlevels(gr, stats::setNames(new[keep], lv[keep]))
}


#' Build the list of UTR introns from the four source annotations
#'
#' The URA GTF has no UTR/CDS information, so the UTR introns are built from the source annotations
#' (CHESS, Ensembl, GENCODE, RefSeq). For every source, UTRs are derived in the same way from the exons,
#' CDS and stop codon of each coding transcript, and the introns lying between two segments of the same UTR are kept.
#' The introns of the four sources are pooled (identical coordinates, strand and UTR type are listed once).
#' For Ensembl, the derivation is checked against its explicit five_prime_utr / three_prime_utr features.
#'
#' @param gtf.paths Named character vector with the paths of the source GTFs: chess, ensembl, gencode, refseq.
#' @param utr.introns.path Where to save the resulting GRanges (RDS), e.g. UTR_introns_hg38.URA.rds.
#'
#' @return GRanges of UTR introns (Ensembl-style seqlevels) with metadata columns utr_type
#' (five_prime_utr, three_prime_utr or both) and sources.
#' @export
# BuildUtrIntronsFromSources(
#   gtf.paths = c(chess = "/home/sg2173/rds/rds-rytenlab-2Hstq49W5EY/sruiz/reference/chess3/chess3.1.3.GRCh38.gtf",
#               ensembl = "/home/sg2173/rds/rds-rytenlab-2Hstq49W5EY/sruiz/reference/ensembl/Homo_sapiens.GRCh38.116.chr.gtf",
#               gencode = "/home/sg2173/rds/rds-rytenlab-2Hstq49W5EY/sruiz/reference/gencode/gencode.v49.annotation.gtf.gz",
#               refseq = "/home/sg2173/rds/rds-rytenlab-2Hstq49W5EY/sruiz/reference/refseq/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"),
#               utr.introns.path = file.path(args$dependencies_folder, "UTR_introns_hg38.URA.rds"))
BuildUtrIntronsFromSources <- function(gtf.paths, utr.introns.path) {
  
  stopifnot("gtf.paths must be a named vector: chess, ensembl, gencode, refseq" =
              setequal(names(gtf.paths), c("chess", "ensembl", "gencode", "refseq")))
  
  per_source <- lapply(c("chess", "ensembl", "gencode", "refseq"), function(src) {
    
    logger::log_info("Reading ", src, " exons and CDS ...")
    feature_types <- c("exon", "CDS", "stop_codon")
    if (src == "ensembl") feature_types <- c(feature_types, "five_prime_utr", "three_prime_utr")
    
    gr <- rtracklayer::import(gtf.paths[[src]], format = "gtf", feature.type = feature_types,
                              colnames = c("type", "transcript_id"))
    
    if (src == "refseq") {
      gr <- .RefseqToEnsemblSeqlevels(gr)
    } else {
      gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse")
      GenomeInfoDb::seqlevelsStyle(gr) <- "Ensembl"
    }
    
    tb <- gr %>% tibble::as_tibble() %>%
      dplyr::mutate(seqnames = as.character(seqnames), strand = as.character(strand),
                    start = as.integer(start), end = as.integer(end))
    rm(gr)
    
    segments <- .UtrSegmentsFromExonsCds(exons = tb %>% dplyr::filter(type == "exon"),
                                         cds = tb %>% dplyr::filter(type %in% c("CDS", "stop_codon")))
    introns <- .UtrIntronsFromSegments(segments) %>% dplyr::mutate(source = src)
    logger::log_info(src, ": ", nrow(introns), " UTR introns derived from exons and CDS (",
                     dplyr::n_distinct(introns$transcript_id), " coding transcripts with UTR introns).")
    
    ## Ensembl: compare with its explicit UTR features
    if (src == "ensembl") {
      explicit_segments <- tb %>%
        dplyr::filter(type %in% c("five_prime_utr", "three_prime_utr")) %>%
        dplyr::transmute(seqnames, strand, transcript_id, utr_type = type, seg_start = start, seg_end = end)
      explicit <- .UtrIntronsFromSegments(explicit_segments)
      key <- function(x) paste(x$seqnames, x$start, x$end, x$strand, x$utr_type)
      k_derived <- unique(key(introns)); k_explicit <- unique(key(explicit))
      logger::log_info("Ensembl check: ", length(k_explicit), " UTR introns from its UTR features, ",
                       length(k_derived), " derived from exons/CDS, ", length(intersect(k_explicit, k_derived)),
                       " in both (", length(setdiff(k_explicit, k_derived)), " only in the UTR features, ",
                       length(setdiff(k_derived, k_explicit)), " only derived).")
    }
    
    gc()
    introns %>% dplyr::select(seqnames, start, end, strand, utr_type, source)
  })
  
  utr_introns <- dplyr::bind_rows(per_source) %>%
    dplyr::group_by(seqnames, start, end, strand) %>%
    dplyr::summarise(n_types = dplyr::n_distinct(utr_type),
                     first_type = dplyr::first(utr_type),
                     sources = paste(sort(unique(source)), collapse = ","),
                     .groups = "drop") %>%
    dplyr::mutate(utr_type = ifelse(n_types > 1, "both", first_type)) %>%
    dplyr::select(seqnames, start, end, strand, utr_type, sources)
  
  logger::log_info(nrow(utr_introns), " distinct UTR introns; ", sum(utr_introns$utr_type == "both"),
                   " are labelled 'both' (5' UTR in some transcripts and 3' UTR in others).")
  
  utr_introns_gr <- GenomicRanges::makeGRangesFromDataFrame(utr_introns, keep.extra.columns = TRUE)
  saveRDS(utr_introns_gr, utr.introns.path)
  logger::log_info("UTR introns saved in ", utr.introns.path)
  
  return(utr_introns_gr)
}
