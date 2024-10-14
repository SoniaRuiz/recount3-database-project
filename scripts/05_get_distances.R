#' Title
#' 
#' @param project.id 
#' @param cluster 
#' @param samples 
#' @param split.read.counts 
#' @param all.split.reads.details 
#' @param folder.name 
#' @param replace 
#'
#' @return
#' @export
#'
#' @examples
GetDistances <- function(project.id,
                         cluster, 
                         samples,
                         split.read.counts,
                         all.split.reads.details,
                         folder.name,
                         replace) {
  
  num_sample <- 1
  
  ProcessSample <- function(sample) {
    
    file_path <- file.path(folder.name, paste0(cluster, "_", sample, "_distances.rds"))
    
    if (file.exists(file_path) && !replace) {
      stop(message(cluster, " sample '", sample, "' exists!"))
    }
    
    if (sum(colnames(split.read.counts) == sample) != 1) {
      stop(message("Error: sample '", sample, "' not found within the split-reads file!"))
    }
    
    print(paste0(Sys.time(), " - ", project.id, " - ", cluster, " sample '", num_sample, "'"))
    
    split_read_counts_sample <- split.read.counts %>%
      as_tibble() %>%
      dplyr::select(junID, all_of(as.character(sample))) %>%
      drop_na() %>%
      filter(.data[[sample]] > 0)
    
    split_reads_details_sample <- all.split.reads.details %>%
      inner_join(split_read_counts_sample, by = "junID") %>%
      dplyr::rename(counts = all_of(as.character(sample))) %>%
      as_tibble()
    
    RunQC(split.reads.details.sample = split_reads_details_sample, 
          split.read.counts.sample = split_read_counts_sample, 
          sample)
    
    distances <- bind_rows(
      ComputeDistances(junction.class = "novel_donor", split.reads.details.sample = split_reads_details_sample, sample),
      ComputeDistances(junction.class = "novel_acceptor", split.reads.details.sample = split_reads_details_sample, sample)
    ) %>% dplyr::select(-sample)
    
    saveRDS(distances, file_path)
    
    num_sample <<- num_sample + 1
    if (num_sample %% 50 == 0) gc()
  }
  
  for (sample in samples) {
    # sample <- samples[1]
    ProcessSample(sample)
  }
}

RunQC <- function(split.reads.details.sample, split.read.counts.sample, sample) {
  
  index <- sample(nrow(split.reads.details.sample), 1)
  data <- split.reads.details.sample[index, c("junID", "counts")]
  data_col <- split.read.counts.sample %>% filter(junID == data$junID)
  
  if (data_col[[sample]] != data$counts) {
    stop("Error: QC failed!")
  }
  
}

ComputeDistances <- function(junction.class, split.reads.details.sample, sample) {
  
  strand_types <- c("+", "-")
  
  results <- lapply(strand_types, function(strand_type) {
    
    # strand_type <- strand_types[1]
    
    junctions <- split.reads.details.sample %>% filter(type == junction.class, strand == strand_type) %>% GRanges()
    annotated <- split.reads.details.sample %>% filter(type == "annotated", strand == strand_type) %>% GRanges()
    
    if ((junction.class == "novel_donor" && strand_type == "+") || (junction.class == "novel_acceptor" && strand_type == "-")){
      overlaps <- GenomicRanges::findOverlaps(query = junctions, subject = annotated, ignore.strand = FALSE, type = "end")
    } else {
      overlaps <- GenomicRanges::findOverlaps(query = junctions, subject = annotated, ignore.strand = FALSE, type = "start")
    }

    novel_junctions <- junctions[queryHits(overlaps),] %>% as_tibble()
    ref_junctions <- annotated[subjectHits(overlaps),] %>% as_tibble()
    
    if (nrow(novel_junctions) == 0) return(NULL)
    
    distance <- if ((junction.class == "novel_donor" && strand_type == "+") || (junction.class == "novel_acceptor" && strand_type == "-")) {
      ref_junctions$start - novel_junctions$start
    } else {
      novel_junctions$end - ref_junctions$end
    }
    
    CreateDistanceDataFrame(novel.junctions = novel_junctions, ref.junctions = ref_junctions, distance, sample, strand_type)
  })
  
  bind_rows(results)
}

CreateDistanceDataFrame <- function(novel.junctions, ref.junctions, distance, sample, type) {
  
  df <- data.frame(
    sample = sample,
    type = type,
    distance = distance,
    novel_junID = novel.junctions$junID,
    novel_counts = novel.junctions$counts,
    novel_seq = novel.junctions$seqnames,
    novel_start = novel.junctions$start,
    novel_end = novel.junctions$end,
    novel_strand = novel.junctions$strand,
    ref_junID = ref.junctions$junID,
    ref_counts = ref.junctions$counts,
    ref_seq = ref.junctions$seqnames,
    ref_start = ref.junctions$start,
    ref_end = ref.junctions$end,
    ref_strand = ref.junctions$strand,
    stringsAsFactors = FALSE
  )
  
  
  df %>%
    group_by(novel_junID) %>%
    filter(abs(distance) == min(abs(distance))) %>%
    filter(ref_counts == max(ref_counts)) %>%
    ungroup()
  
}
