# library(tidyverse)
# dependencies.folder = "~/rds/hpc-work/SR/recount3-database-project/dependencies"

GenerateClinvarData <- function(dependencies.folder){
  
  ## Download "clinvar.vcf" from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
  clinvar_data <- read.table(file = file.path(dependencies.folder, "clinvar.vcf"), header = F, sep = "\t") %>% as_tibble
  
  clinvar_likely_pathogenic <- clinvar_data %>%
    filter((str_detect(string = V8, pattern = "CLNSIG=Likely_pathogenic") | str_detect(string = V8, pattern = "CLNSIG=Pathogenic")) & 
             (str_detect(string = V8, pattern = "splice_acceptor_variant") | str_detect(string = V8, pattern = "splice_donor_variant") | str_detect(string = V8, pattern = "missense_variant")))
  
  clinvar_likely_pathogenic_tidy <- clinvar_likely_pathogenic %>%
    separate_rows(V8, sep = "\\;\\s*") %>%
    separate(V8, into = c('col1', 'col2'), sep = "=") %>% 
    pivot_wider(names_from = col1, values_from = col2) 
  
  clivar_tidy_to_save <- clinvar_likely_pathogenic_tidy %>%
    dplyr::rename(seqnames = V1,
                  start = V2, 
                  ID = V3,
                  REF =V4,
                  ALT = V5,
                  QUAL = V6,
                  FILTER = V7)%>%
    mutate(end = start)
  
  saveRDS(object = clivar_tidy_to_save, file = file.path(dependencies.folder, "clinvar_splicing_pathogenic.rds"))
  
  return(clivar_tidy_to_save)
  
}

AddClinvarData <- function(db.introns,
                           dependencies.folder) {
  
  clinvar_gr <- if (!file.exists(file.path(dependencies.folder, "clinvar_splicing_pathogenic.rds"))) {
    GenerateClinvarData(dependencies.folder)
  } else {
    readRDS(file = paste0(dependencies.folder, "/clinvar_splicing_pathogenic.rds")) 
  }
  clinvar_gr <- clinvar_gr %>% GenomicRanges::GRanges()
  
  seqlevelsStyle(clinvar_gr) <- "UCSC"
  elementMetadata(clinvar_gr)[, "ID"] <- (clinvar_gr) %>% as.character()
  
  df_all_introns_gr <- db.introns %>% mutate(clinvar = F) %>% GRanges() 
  
  ## Find overlaps between clinvar mutations and mis-splicing ratios
  overlaps <- GenomicRanges::findOverlaps(query = clinvar_gr, subject = df_all_introns_gr, type = "any")
  
  elementMetadata(clinvar_gr)[queryHits(overlaps), "intron_ID"] <- subjectHits(overlaps)
  elementMetadata(df_all_introns_gr)[subjectHits(overlaps), "clinvar"] <- T
  
  clinvar_hits <- (clinvar_gr)[queryHits(overlaps),] %>%
    as_tibble() %>%
    dplyr::select(seqnames,start,end, width, strand, ID, CLNSIG, CLNVC, MC, intron_ID) %>%
    dplyr::group_by(intron_ID) %>%
    mutate(ID_list = paste(ID, collapse = ",")) %>%
    mutate(CLNSIG_list = paste(CLNSIG, collapse = ",")) %>%
    mutate(CLNVC_list = paste(CLNVC, collapse = ",")) %>%
    mutate(MC_list = paste(MC, collapse = ",")) %>%
    ungroup() 
  
  elementMetadata(df_all_introns_gr)[subjectHits(overlaps), "clinvar"] <- T
  
  logger::log_info(df_all_introns_gr %>% as_tibble() %>% filter(clinvar == T) %>%  nrow(), " introns containing ClinVar variants")
  
  return(df_all_introns_gr)
  
}
