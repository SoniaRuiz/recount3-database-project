library(tidyverse)
library(GenomicRanges)
library(DBI)

####################################################
## CONNECT TO THE SPLICING DATABASE ################
####################################################


## CONNECT TO THE DATABASE ------------------------------

setwd(normalizePath("."))

gtf_version <- 105
project_name <- "SRP100948"

database_folder <- paste0(getwd(), "/database/", project_name, "/", gtf_version, "/")
results_folder <- paste0(getwd(), "/results/", project_name, "/", gtf_version, "/")

results_path <- paste0(results_folder, "/_paper/results/")
figures_path <- paste0(results_folder, "/_paper/figures/")





database_path <- paste0(database_folder, "/", project_name, ".sqlite")


con <- dbConnect(RSQLite::SQLite(), database_path)
dbListTables(con)


## GET FROM MASTER TABLE
query = paste0("SELECT * FROM 'master'")
df_metadata <- dbGetQuery(con, query) %>% as_tibble()

all_projects <- df_metadata$SRA_project %>% unique
all_clusters <- df_metadata$cluster %>% unique()


get_mode <- function(data) {
  uniqv <- unique(data)
  uniqv[which.max(tabulate(match(data, uniqv)))]
}


custom_ggtheme <-  theme(text = element_text(size = 7, family="Arial", colour = "black"),
                         
                         axis.ticks = element_line(colour = "black", linewidth = 2),
                         axis.text = element_text(size = 7, family="Arial", colour = "black"),
                         axis.line = element_line(colour = "black"),
                         axis.title = element_text(size = 7, family="Arial", colour = "black"),
                         axis.text.y = element_text(size = 7, family="Arial", colour = "black"),
                         axis.text.x = element_text(size = 7, family="Arial", colour = "black", 
                                                    hjust = 0.5, vjust = 0.5),
                         strip.text = element_text(size = 7, family="Arial", colour = "black"),
                         
                         legend.text = element_text(size = "7", family="Arial", colour = "black"),
                         legend.title = element_blank(),
                         legend.position = "top",
                         legend.box = "vertical")

########################################
## FUNCTIONS 
########################################


## SECTION 1 ---------------------------------------------

#' Title
#' Obtain general stats about the database
#' @return
#' @export
#'
#' @examples
get_database_stats <- function() {
  
 
  ## Methods: number of samples by RIN number
  df_metadata %>% nrow() 
  
  df_metadata %>%
    filter(rin >= 8) %>% nrow()
  df_metadata %>%
    filter(rin >= 7) %>% nrow()
  df_metadata %>%
    filter(rin >= 6) %>% nrow()
  
  #if ( df_metadata %>% filter(rin < 6) %>% nrow() > 1 ) {
  #  print("ERROR! Some of the samples considered present a RIN number lower than 6!")
  #  break;
  #}
  

  query <- paste0("SELECT * from 'intron'")
  db_introns <- dbGetQuery(con, query) 
  db_introns %>% head()
  db_introns %>% distinct(ref_junID) %>% nrow()
  db_introns %>%
    dplyr::count(misspliced)
  
  
  ## Collectively, we detected 3,865,268 unique novel junctions, equating to 14 novel junctions per annotated intron. 
  query <- paste0("SELECT * from 'novel'")
  db_novel <- dbGetQuery(con, query) 
  db_novel %>% head
  db_novel %>% nrow() 
  
  db_novel %>% 
    dplyr::count(novel_type)
  
  (db_novel %>% nrow()) / (db_introns %>% filter(misspliced==1) %>% nrow())
  
  
}

#' Title
#' Get an overview of the metadata of the project
#' @return
#' @export
#'
#' @examples
plot_metadata <- function() {
  

  
  ## Num samples
  plot_num_samples <- ggplot(df_metadata %>% 
                               dplyr::count(SRA_project, cluster) %>%
                               arrange(cluster , n) %>%
                               mutate(SRA_project = fct_inorder(SRA_project)) ) +
    geom_bar(aes(y = SRA_project, x = n, fill = cluster),
             stat = "identity", position = position_dodge()) + 
    theme_light() +
    scale_fill_hue() +
    labs(y = "", x = "Num. samples" ) +
    guides(fill = guide_legend(title = "Age group: ", ncol = 3, nrow = 1)) +
    custom_ggtheme
  
  plot_num_samples
  
  
  ## Gender
  plot_gender <- ggplot(df_metadata %>% 
                          dplyr::rename(gender = Sex) %>%
                          mutate(gender = gender %>% as.character()) %>%
                          dplyr::count(SRA_project, gender) %>%
                          arrange(gender , n) %>%
                          mutate(SRA_project = fct_inorder(SRA_project)) ) +
    geom_bar(aes(y = SRA_project, x = n, group = gender, fill = gender),
             stat = "identity", position = "dodge") + 
    theme_light() +
    labs(y = "", x = "Num. samples" ) +
    scale_fill_manual(labels = c("Male", "Female"), values = c("1", "2"), palette=scales::hue_pal()) +
    guides(fill = guide_legend(title = "Gender: ", ncol = 2, nrow = 1)) +
    custom_ggtheme
  
  plot_gender
  
  
  ## RIN
  plot_rin <- ggplot(df_metadata ) +
    geom_density(aes(x = rin_score, fill = cluster), alpha = 0.5) + 
    theme_light() +
    labs(x = "RIN" ) +
    scale_fill_hue() +
    #guides(fill = guide_legend(title = "Age group: ", ncol = 3, nrow = 1)) +
    custom_ggtheme
  
  plot_rin
  
  
  ## Mapped read depth
  plot_read_depth <- ggplot(df_metadata ) +
    geom_density(aes(x = all_mapped_reads, fill = cluster), alpha = 0.5) + 
    theme_light() +
    scale_fill_hue() +
    labs(x = "Mapped Read Count" ) +
    guides(fill = guide_legend(title = "Age group: ", ncol = 3, nrow = 1)) +
    custom_ggtheme
  
  plot_read_depth
  
  ggpubr::ggarrange(plot_num_samples,plot_gender,
                    plot_rin,plot_read_depth,
                    labels = c("a","b","c","d"))
  
  dir.create(path = figures_path, recursive = T,showWarnings = T)
  ggplot2::ggsave(paste0(figures_path, "/metadata.png"), width = 180, height = 100, units = "mm", dpi = 300)
  
  
}

#' Title
#' Get common introns between cases and controls (ONLY FRONTAL CORTEX) and subsample them by read coverage
#' @return
#' @export
#'
#' @examples
get_common_subsample_introns <- function() {
  
  ###############################
  ## GET DATA FOR CONTROL
  ###############################
  
  cluster_id <- "control"
  
  print(paste0(Sys.time(), " - ", project_id, " - ", cluster_id))
  
  query <- paste0("SELECT DISTINCT ref_junID, MSR_D, MSR_A, ref_type, ref_sum_counts, ref_n_individuals FROM '", 
                  cluster_id, "_", project_id, "_nevermisspliced'")
  introns <- dbGetQuery(con, query) %>% as_tibble()
  introns %>% nrow()
  query <- paste0("SELECT DISTINCT ref_junID, MSR_D,  MSR_A, ref_type, ref_sum_counts, ref_n_individuals FROM '", 
                  cluster_id, "_", project_id, "_misspliced'")
  introns_control <- dbGetQuery(con, query) %>% as_tibble()
  introns_control %>% nrow()
  introns_control <- rbind(introns, introns_control)
  #introns_control <- dbGetQuery(con, query) %>% as_tibble()
  #query <- paste0("SELECT DISTINCT ref_junID, protein_coding, lncRNA FROM 'intron' WHERE ref_junID IN (",
  #                paste(introns$ref_junID, collapse = ","),")")
  #introns_control <- introns %>%
  #  left_join(y = dbGetQuery(con, query) %>% as_tibble(),
  #            by = "ref_junID") %>% 
  #  as_tibble() 
  
  ###############################
  ## GET DATA FOR PD
  ###############################
  
 
  cluster_id <- "AD"
  
  print(paste0(Sys.time(), " - ", project_id, " - ", cluster_id))
  
  query <- paste0("SELECT DISTINCT ref_junID, MSR_D, MSR_A, ref_type, ref_sum_counts, ref_n_individuals FROM '", 
                  cluster_id, "_", project_id, "_nevermisspliced'")
  introns <- dbGetQuery(con, query) %>% as_tibble()
  introns %>% nrow()
  query <- paste0("SELECT DISTINCT ref_junID, MSR_D,  MSR_A, ref_type, ref_sum_counts, ref_n_individuals FROM '", 
                  cluster_id, "_", project_id, "_misspliced'")
  introns_pd <- dbGetQuery(con, query) %>% as_tibble()
  introns_pd %>% nrow()
  introns_pd <- rbind(introns, introns_pd)
  #introns_pd <- dbGetQuery(con, query) %>% as_tibble()
  
  
  
  ##############################################################################
  ## TIDY DATAFRAME - ONLY USE COMMON INTRONS
  ##############################################################################
  
  common_introns <- intersect(introns_pd$ref_junID,
                              introns_control$ref_junID) %>% unique()
  common_introns %>% length()
  
  ## TIDY DATAFRAME
  df_introns_tidy <- rbind(introns_pd %>%
                             mutate(sample_type = "AD"), 
                           introns_control %>%
                             mutate(sample_type = "control"))  %>%
    filter(ref_junID %in% common_introns) %>%
    #group_by(sample_type, ref_junID) %>%
    mutate(mean_coverage = ref_sum_counts/ref_n_individuals) %>%
    mutate(mean_coverage = mean_coverage %>% log10())# %>%
  #ungroup()
  df_introns_tidy
  
  
  ## PLOT COVERAGE OF SHARED INTRONS BEFORE SUBSAMPLING
  
  plot_BS <- ggplot(data = df_introns_tidy) +
    geom_density(mapping = aes(x = mean_coverage, fill = sample_type), alpha = 0.8) +
    ggtitle("Before subsampling") +
    custom_ggtheme +
    scale_fill_discrete(name = "Transcript biotype: ") +
    xlab("log10 mean read coverage")
  
  plot_BS
  
  
  
  
  #####################################
  #### START SUBSAMPLIING
  #####################################
  
  
  print(paste0(Sys.time(), " - starting subsampling ... "))
  
  ## Subsampling introns to control by similarity in mean read coverage
  m.out <- MatchIt::matchit(sample_type ~ mean_coverage,
                            data = df_introns_tidy %>%
                              mutate(sample_type = sample_type %>% as.factor()),
                            distance = df_introns_tidy$mean_coverage,
                            method = "nearest",
                            caliper = c(mean_coverage = 0.0005),
                            std.caliper = FALSE)
  subsample <- MatchIt::match.data(m.out)
  subsample %>% dplyr::count(sample_type)
  
  ## PLOT COVERAGE AFTER SUBSAMPLE
  
  plot_AS <- ggplot(data = subsample) +
    geom_density(mapping = aes(x = mean_coverage, fill = sample_type), alpha = 0.8) +
    ggtitle("After subsampling") +
    #ggtitle("Mean read coverage per annotated intron across all samples\nfrom 54 GTEx v8 tissues - Subsampling performed.") +
    scale_fill_discrete(name = "Transcript biotype: ") +
    custom_ggtheme +
    xlab("log10 mean read coverage")
  
  ggpubr::ggarrange(plot_BS,
                    plot_AS,
                    labels = c("a", "b"),
                    common.legend = T)
  

  dir.create(path = figures_path, recursive = T,showWarnings = T)
  ggplot2::ggsave(paste0(figures_path, "/subsampling.png"), width = 180, height = 100, units = "mm", dpi = 300)
  
  saveRDS(object = subsample,
          file = paste0(results.folder, "/common_subsampled_introns.rds"))
}

#' Title
#' Test for differences in MSR_D and MSR_A between the introns from AD vs control samples.
#' Using one-tailed paired Wilcoxon test
#' @return
#' @export
#'
#' @examples
test_MSR_differences <- function()  {
  
  
  common_introns_subsample <- readRDS(file = paste0(results_folder, "/common_subsampled_introns.rds"))
 
  
  #########################################
  ## TEST DONOR - AFTER SUBSAMPLING
  #########################################

  wilcox.test(x = common_introns_subsample %>% filter(sample_type == "AD") %>% pull(MSR_D),
              y = common_introns_subsample %>% filter(sample_type == "control") %>% pull(MSR_D),
              alternative = "greater",
              paired = T,
              correct = T)
  
  common_introns_subsample %>% filter(sample_type == "AD") %>% pull(MSR_D) %>% summary
  common_introns_subsample %>% filter(sample_type == "control") %>% pull(MSR_D) %>% summary
  
  
  
  ########################################
  ## TEST ACCEPTOR
  ########################################

  
  wilcox.test(x = common_introns_subsample %>% filter(sample_type == "AD") %>% pull(MSR_A),
              y = common_introns_subsample %>% filter(sample_type == "control") %>% pull(MSR_A),
              alternative = "greater",
              correct = T)
  
  
  common_introns_subsample %>% filter(sample_type == "AD") %>% pull(MSR_A) %>% summary
  common_introns_subsample %>% filter(sample_type == "control") %>% pull(MSR_A) %>% summary
  
  
  #######################################
  ## EXACT PAIRING
  #######################################
  
  common_introns_subsample_column <- common_introns_subsample %>% 
    filter(sample_type == "AD") %>% 
    dplyr::select(ref_junID, MSR_A_AD = MSR_A) %>%
    inner_join(y = common_introns_subsample %>% 
                 filter(sample_type == "control") %>% 
                 dplyr::select(ref_junID, MSR_A_control = MSR_A) )
  
  wilcox.test(x = common_introns_subsample_column$MSR_A_AD,
              y = common_introns_subsample_column$MSR_A_control,
              paired = T,
              alternative = "greater",
              correct = T)

}


#' Title
#' Visualise differences in MSR_D and MSR_A values in introns from AD vs control samples
#' Only using common, subsampled introns
#' @return
#' @export
#'
#' @examples
plot_distances <- function() {
  
  limit_bp <- 30
  common_introns_subsample <- readRDS(file = paste0(results_folder, "/common_subsampled_introns.rds")) %>%
    dplyr::select(-distance,-weights,-subclass)
  
  ###############################
  ## GET novel junction data and join with common introns
  ###############################
  
  df_misspliced <- map_df(all_clusters, function(cluster_id) {
    
    query <- paste0("SELECT tissue.novel_junID, tissue.ref_junID, novel.novel_type, novel.distance 
                  FROM '", cluster_id, "_", project_name, "_misspliced' AS tissue 
                  INNER JOIN 'novel' ON novel.novel_junID = tissue.novel_junID")
    df_misspliced <- dbGetQuery(con, query) %>% 
      as_tibble() %>%
      filter(ref_junID %in% common_introns_subsample$ref_junID) %>%
      mutate(cluster = cluster_id)
  })
  

  ################################
  ## PLOT DISTANCES
  ################################
  
  ggplot(data = df_misspliced %>%
           filter(abs(distance) <= limit_bp)) + 
    geom_histogram(aes(x = distance, fill = novel_type),
                   bins = limit_bp * 2,
                   binwidth = 1,
                   position = "stack") +
    xlab("Distance (in bp)") +
    ylab("Unique novel junctions") +
    facet_grid(novel_type~cluster) +
    theme_light() +
    scale_x_continuous(limits = c(limit_bp,(limit_bp * -1)),
                       breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6),
                       trans = "reverse") +
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      breaks = c("novel_donor", "novel_acceptor")) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1 )) +
    custom_ggtheme +
    theme(legend.position = "none") 
  
  ggplot2::ggsave(paste0(figures_path, "/distances.png"), width = 180, height = 100, units = "mm", dpi = 300)
  
}


#' Title
#' Visualise differences in modulo3 values in introns from AD vs control samples
#' Only using common, subsampled introns
#' @return
#' @export
#'
#' @examples
plot_modulo <- function() {
  
  common_introns_subsample <- readRDS(file = paste0(results_folder, "/common_subsampled_introns.rds")) %>%
    dplyr::select(-distance,-weights,-subclass)
  
  common_introns_subsample$ref_junID %>% unique() %>% length()
  
  ## Calculate modulo3 from MANE transcripts in 100bp distance
  df_modulo <- map_df(all_clusters, function(cluster_id) {
      
    # cluster_id <- all_clusters[1]
    print(paste0(Sys.time(), " - ", cluster_id))
    
    #####################################
    ## Get the novel junctions from the current tissue
    #####################################
    
    query <- paste0("SELECT novel_junID 
                    FROM '", cluster_id, "_", project_name, "_misspliced'")
    introns <- dbGetQuery(con, query) %>% as_tibble()
      
    ## Get the distance location
    query <- paste0("SELECT * FROM 'novel' WHERE novel_junID IN (",
                    paste(introns$novel_junID, collapse = ","),")")
    introns <- introns %>%
      left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                by = "novel_junID") %>% 
      as_tibble() 
      
    ## Add the transcript and MANE info
    query <- paste0("SELECT intron.ref_junID, intron.protein_coding, transcript.MANE
                    FROM 'intron' 
                    INNER JOIN transcript
                    ON intron.transcript_id = transcript.id
                    WHERE ref_junID IN (", paste(introns$ref_junID, collapse = ","),") 
                    AND transcript.MANE = 1")
    introns <- introns %>%
      inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
                by = "ref_junID") %>% 
      as_tibble() 
    
    df_novel_tidy <- introns %>%
      distinct(novel_junID, .keep_all = T) %>%
      filter(abs(distance) <= 100, MANE == 1) %>% 
      mutate(novel_type = str_replace(string = novel_type,
                                      pattern = "_",
                                      replacement = " ")) %>%
      mutate(type_p = ifelse(distance < 0, paste0(novel_type," intron"), paste0(novel_type," exon"))) %>% 
      mutate(modulo = abs(distance) %% 3) %>%
      ## Only common, subsampled introns
      filter(ref_junID %in% common_introns_subsample$ref_junID) 
      
    df_novel_tidy <- df_novel_tidy %>% 
      group_by(modulo) %>%
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>%
      mutate(cluster = cluster_id)
    
    return(df_novel_tidy)
  })
  
  
  df_modulo$modulo = factor(df_modulo$modulo, 
                            levels = c( "0", "1", "2"))
  
  
  ################
  ## DENSITY PLOT
  ################

  df_modulo_tidy <- df_modulo %>%
    mutate(freq = freq * 100)
  
  ggplot(df_modulo_tidy, aes(x = cluster, y = freq)) +
    geom_bar(stat = "identity") +
    xlab("Sample type") +
    ylab("% of novel junctions") +
    facet_grid(~modulo) +
    theme_light() +
    custom_ggtheme +
    theme(legend.position = "none") 
  
  
  ggplot2::ggsave(paste0(figures_path, "/modulo.png"), width = 180, height = 100, units = "mm", dpi = 300)
  
}


plot_unique_junctions <- function() {
  
  ## Load common introns
  common_introns_subsample <- readRDS(file = paste0(results_folder, "/common_subsampled_introns.rds")) %>%
    dplyr::select(-distance,-weights,-subclass)
  common_introns_subsample$ref_junID %>% unique() %>% length()
  

  df_unique_junctions <- map_df(all_clusters, function(cluster_id) {

    print(cluster_id)
  
    ####################
    ## GET THE INTRONS
    ####################
    
    query <- paste0("SELECT DISTINCT ref_junID FROM '", 
                    cluster_id, "_", project_name, "_nevermisspliced'")
    introns <- dbGetQuery(con, query) %>% as_tibble()
    query <- paste0("SELECT DISTINCT ref_junID FROM '", 
                    cluster_id, "_", project_name, "_misspliced'")
    introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble()) %>%
      ## Only common, subsampled introns
      filter(ref_junID %in% common_introns_subsample$ref_junID)
    
    
    
    ###########################
    ## GET THE NOVEL JUNCTIONS
    ###########################
    
    query <- paste0("SELECT ref_junID, novel_junID
                    FROM '", cluster_id, "_", project_name, "_misspliced'")
    novel_junctions <- dbGetQuery(con, query) %>% as_tibble() 
    query <- paste0("SELECT novel_junID, novel_type FROM 'novel' WHERE novel_junID IN (", 
                    paste(novel_junctions$novel_junID, collapse = ","), ")")
    novel_junctions <- novel_junctions %>%
      left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                by = "novel_junID") %>% 
      as_tibble() %>%
      ## Only common, subsampled introns
      filter(ref_junID %in% common_introns_subsample$ref_junID)
    
    
    ###########################
    ## GET THE PROPORTIONS
    ###########################
    
    
    annotated_junc <- introns %>% distinct(ref_junID) %>% nrow()
    donor_junc <- novel_junctions %>% filter(novel_type == "novel_donor") %>% distinct(novel_junID) %>% nrow()
    acceptor_junc <- novel_junctions %>% filter(novel_type == "novel_acceptor") %>% distinct(novel_junID) %>% nrow()
    
    annotated_prop <- annotated_junc/(annotated_junc + donor_junc + acceptor_junc)
    donor_prop <- donor_junc/(annotated_junc + donor_junc + acceptor_junc)
    acceptor_prop <- acceptor_junc/(annotated_junc + donor_junc + acceptor_junc)
    
    
    
    ## Return the data.frame
    return(data.frame(cluster = cluster_id,
                      annotated_junc = annotated_junc ,
                      donor_junc = donor_junc,
                      acceptor_junc = acceptor_junc,
                      annotated_prop = annotated_prop,
                      donor_prop = donor_prop,
                      acceptor_prop = acceptor_prop))
    
    
    
  })
  

  ######################
  ## BAR PLOT
  ######################
  
  
  df_unique_junctions <- df_unique_junctions %>%
    dplyr::select(cluster, 
                  donor = donor_prop, 
                  acceptor = acceptor_prop, 
                  annotated_intron = annotated_prop) %>%
    tidyr::gather(key = "type", value = "prop", -cluster ) 
  
  df_unique_junctions <- df_unique_junctions %>%
    filter(type != "annotated_intron") %>%
    mutate(prop = prop * 100)
  
  df_unique_junctions$type = factor(df_unique_junctions$type, 
                                      levels = c("donor","acceptor"))
  
  
  ggplot(df_unique_junctions, aes(x = cluster, y = prop, fill = type)) + 
    geom_bar(stat = "identity") +
    facet_grid(~type) +
    theme_light() +
    ylab("% unique novel junctions") +
    xlab("") +   
    custom_ggtheme  +
    theme(axis.line = element_line(colour = "black"), 
          text = element_text(colour = "black", size = 12),
          legend.position = "top") +
    ggsci::scale_fill_npg()
  
  
  ggplot2::ggsave(paste0(figures_path, "/unique_novel_junctions.png"), width = 180, height = 100, units = "mm", dpi = 300)
  
  
  
}

#' Title
#' Visualise differences in number of novel reads in common introns from AD vs control samples
#' Only using common, subsampled introns
#' @return
#' @export
#'
#' @examples
plot_unique_reads <- function() {
    
  ## Load common introns
  common_introns_subsample <- readRDS(file = paste0(results_folder, "/common_subsampled_introns.rds")) %>%
    dplyr::select(-distance,-weights,-subclass)
  common_introns_subsample$ref_junID %>% unique() %>% length()
    
    
  df_mean_counts <- map_df(all_clusters, function(cluster_id) {
      
    # project_id <- "BRAIN"
    # cluster_id <- "Brain - Frontal Cortex (BA9)"
    print(cluster_id)
    
    ####################
    ## GET THE INTRONS
    ####################
    
    query <- paste0("SELECT DISTINCT ref_junID, ref_sum_counts, ref_n_individuals FROM '", 
                    cluster_id, "_", project_name, "_nevermisspliced'")
    introns <- dbGetQuery(con, query) %>% as_tibble()
    query <- paste0("SELECT DISTINCT ref_junID, ref_sum_counts, ref_n_individuals FROM '", 
                    cluster_id, "_", project_name, "_misspliced'")
    introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble()) %>%
      ## Only common, subsampled introns
      filter(ref_junID %in% common_introns_subsample$ref_junID)
    
    
    
    ###########################
    ## GET THE NOVEL JUNCTIONS
    ###########################
    
    query <- paste0("SELECT ref_junID, novel_junID, novel_sum_counts, novel_n_individuals 
                    FROM '", cluster_id, "_", project_name, "_misspliced'")
    novel_junctions <- dbGetQuery(con, query) %>% as_tibble() 
    query <- paste0("SELECT novel_junID, novel_type FROM 'novel' WHERE novel_junID IN (", 
                    paste(novel_junctions$novel_junID, collapse = ","), ")")
    novel_junctions <- novel_junctions %>%
      left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                by = "novel_junID") %>% 
      as_tibble() %>%
      ## Only common, subsampled introns
      filter(ref_junID %in% common_introns_subsample$ref_junID)
    
    
    ###########################
    ## GET THE PROPORTIONS
    ###########################
    
    
    annotated <- introns %>%
      dplyr::distinct(ref_junID, .keep_all = T) %>%
      pull(ref_sum_counts) %>% 
      sum()
    
    acceptor <- novel_junctions %>%
      filter(novel_type == "novel_acceptor") %>%
      dplyr::distinct(novel_junID, .keep_all = T) %>%
      pull(novel_sum_counts) %>% 
      sum()
    
    donor <- novel_junctions %>%
      filter(novel_type == "novel_donor") %>%
      dplyr::distinct(novel_junID, .keep_all = T) %>%
      pull(novel_sum_counts) %>% 
      sum()
    
    annotated_p = annotated * 100 / (annotated + acceptor + donor)
    acceptor_p = acceptor * 100 / (annotated + acceptor + donor)
    donor_p = donor * 100 / (annotated + acceptor + donor)
    
    return(data.frame(cluster = cluster_id,
                      type = c("annotated","acceptor", "donor"),
                      prop_counts = c(annotated_p, acceptor_p, donor_p),
                      sum_counts = c(annotated, acceptor, donor)))
    
  
    
  })
  
  ######################
  ## BAR PLOT
  ######################
  
  df_mean_counts_violin <- df_mean_counts %>%
    filter(type != "annotated") 
  
  df_mean_counts_violin$type = factor(df_mean_counts_violin$type, 
                                      levels = c("donor", "acceptor"))
  
  
  ggplot(df_mean_counts_violin, aes(x = cluster, y = prop_counts, fill = type)) + 
    geom_bar(stat = "identity") +
    facet_grid(~type) +
    theme_light() +
    ylab("% cumulative novel read counts") +
    xlab("") +   
    custom_ggtheme  +
    theme(axis.line = element_line(colour = "black"), 
          text = element_text(colour = "black", size = 12),
          legend.position = "top") +
    ggsci::scale_fill_npg()
  
  
  ggplot2::ggsave(paste0(figures_path, "/read_counts.png"), width = 180, height = 100, units = "mm", dpi = 300)
  
  
  
}

##################################
## CALLS
##################################



# compare_tissues_somatic_mutations(project_id1 = "BRAIN",
#                                   cluster_id1 = "Brain - Nucleus accumbens (basal ganglia)",
#                                   project_id2 = "BLOOD",
#                                   cluster_id2 = "Whole Blood")
  