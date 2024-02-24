library(tidyverse)
library(GenomicRanges)
library(DBI)

####################################################
## CONNECT TO THE SPLICING DATABASE ################
####################################################


## CONNECT TO THE DATABASE ------------------------------

setwd(normalizePath("."))
# setwd(normalizePath("/mnt/PROJECTS/recount3-database-project/"))

gtf_version <- 105
supportive_reads <- 1
#project_id <- "SRP058181" 
project_id <- "SRP100948"
project_name <- paste0(project_id, "_", supportive_reads, "read")

database_folder <- paste0(getwd(), "/database/", project_name, "/", gtf_version, "/")
results_folder <- paste0(getwd(), "/results/", project_name, "/", gtf_version, "/")

results_path <- paste0(results_folder, "/_paper/results/")
figures_path <- paste0(results_folder, "/_paper/figures/")


database_path <- paste0(database_folder, "/", project_name, ".sqlite")


con <- dbConnect(RSQLite::SQLite(), database_path)
dbListTables(con)


## QUERY MASTER TABLES 

query = paste0("SELECT * FROM 'metadata'")
df_metadata <- dbGetQuery(con, query) %>% distinct(external_id, .keep_all = T) %>% as_tibble()
all_projects <- df_metadata$SRA_project %>% unique

query <- paste0("SELECT * FROM 'intron'")
master_introns <- dbGetQuery(con, query) %>% as_tibble()

query <- paste0("SELECT * FROM 'novel'")
master_novel_junctions <- dbGetQuery(con, query) %>% as_tibble()

query <- paste0("SELECT * FROM 'transcript'")
master_transcripts <- dbGetQuery(con, query) %>% as_tibble()

query <- paste0("SELECT * FROM 'gene'")
master_gene <- dbGetQuery(con, query) %>% as_tibble()

all_projects <- df_metadata$SRA_project %>% unique
all_clusters <- df_metadata$cluster %>% unique()


## UTILS

get_mode <- function(data) {
  uniqv <- unique(data)
  uniqv[which.max(tabulate(match(data, uniqv)))]
}


custom_ggtheme <-  theme(text = element_text(size = 7, family="Arial", colour = "black"),
                         axis.ticks = element_line(colour = "black"),
                         axis.text = element_text(size = 7, family="Arial", colour = "black"),
                         axis.line = element_line(colour = "black"),
                         axis.title = element_text(size = 7, family="Arial", colour = "black"),
                         axis.text.y = element_text(size = 7, family="Arial", colour = "black"),
                         axis.text.x = element_text(size = 7, family="Arial", colour = "black", hjust = 0.5, vjust = 0.5),
                         strip.text = element_text(size = 7, family="Arial", colour = "black"),
                         legend.text = element_text(size = 7, family="Arial", colour = "black"),
                         legend.position = "top",
                         legend.box = "vertical")

########################################
## FUNCTIONS 
########################################


## SECTION 1 ---------------------------------------------

dummie_test <- function() {
  
  #######################################
  ## Test that the read counts assigned to a random intron 
  #######################################
  
  
  ## Find in the original split read data whether the counts assigned for that annotated intron
  ## or novel junction whithin the database are the same within the original split read counts data
  
  
  cluster_id <- "control"
  table_name <- paste0(cluster_id, "_", project_id, "_misspliced")
  
  query = paste0("SELECT * FROM '", table_name, "' ")
  df_missplicing <- dbGetQuery(con, query) 
  
  ## Choose a random annotated intron
  junID <- sample(x = df_missplicing$ref_junID, size = 1)
  
  split_read_counts <- readRDS(file = paste0(results_folder, "/", project_id, 
                                             "/base_data/",project_id, "_", cluster_id, "_split_read_counts.rds"))
  
  if ( is.null(names(split_read_counts)) ) { 
    split_read_counts <- split_read_counts %>% 
      as_tibble(rownames = "junID")
  }
  
  #############
  
  query = paste0("SELECT ref_junID, ref_coordinates 
                 FROM 'intron' 
                 WHERE ref_junID == ", junID)
  df_intron <- dbGetQuery(con, query) 
  df_intron
  
  query = paste0("SELECT * FROM '", table_name, "' WHERE ref_junID == ", junID)
  df_missplicing <- dbGetQuery(con, query) 
  df_missplicing
  
  query = paste0("SELECT novel_junID, novel_coordinates 
                 FROM 'novel' 
                 WHERE novel_junID IN (", 
                 paste(df_missplicing$novel_junID, collapse = ","), ")")
  df_novel <- dbGetQuery(con, query) 
  df_novel
  
  df_merged <- merge(x = df_missplicing,
                     y = df_novel,
                     by = "novel_junID")
  
  df_merged <- df_merged %>%
    inner_join(y = df_intron,
               by = c("ref_junID"))
  
  
  original_counts_ref_junID <- split_read_counts %>%
    dplyr::filter(junID == (df_merged$ref_coordinates %>% unique())) %>%
    mutate(sum = rowSums(across(where(is.numeric)))) %>%
    dplyr::select(junID, sum)
  
  if ((df_merged$ref_sum_counts %>% unique()) != original_counts_ref_junID$sum) {
    print(paste0("ERROR! The number of counts stored on IntroVerse for the intron '", 
                 junID, "' differ from the original figure."))
  }
  
  
  ## Choose a random novel junction
  novel_jxn <- sample(df_merged$novel_coordinates, size = 1)
  
  ## Test the read counts with the novel junction table
  original_counts_novel_junID <- split_read_counts %>%
    dplyr::filter(junID == novel_jxn) %>%
    mutate(sum = rowSums(across(where(is.numeric)))) %>%
    dplyr::select(junID, sum)
  original_counts_novel_junID
  
  ## Ideally, this should be tested for all the novel junctions.
  if ((df_merged %>%
       filter(novel_coordinates == novel_jxn) %>% 
       pull(novel_sum_counts)) != original_counts_novel_junID$sum) {
    print(paste0("ERROR! The number of counts stored on IntroVerse for the first novel junction attached to the intron '", 
                 junID, "' differ from the original figure."))
  }
}


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
  

  master_introns %>% head()
  master_introns %>% distinct(ref_junID) %>% nrow()
  master_introns %>%
    dplyr::count(misspliced)
  
  
  ## Collectively, we detected 3,865,268 unique novel junctions, equating to 14 novel junctions per annotated intron. 
  master_novel_junctions %>% head
  master_novel_junctions %>% nrow() 
  
  master_novel_junctions %>% 
    dplyr::count(novel_type)
  
  (master_novel_junctions %>% nrow()) / (master_introns %>% filter(misspliced==1) %>% nrow())
  
  
  master_transcripts %>% distinct(transcript_id) %>% nrow()
  master_gene %>% distinct(gene_id) %>% nrow()

}


#' Title
#' Get an overview of the metadata of the project
#' @return
#' @export
#'
#' @examples
plot_metadata <- function() {
  

  
  # df_metadata$cluster = factor(df_metadata$cluster, 
  #                                   levels = c("control", "AD"))
  ## Num samples
  plot_num_samples <- ggplot(df_metadata %>% 
                               # mutate(cluster = factor(cluster, 
                               #                         levels = c("control", "AD"))) %>%
                               dplyr::count(SRA_project, cluster)  ) +
    geom_bar(aes(x = n, 
                 y = SRA_project,
                 fill = cluster),
             stat = "identity", position = position_dodge()) + 
    theme_light() +
    scale_fill_hue() +
    labs(y = "", x = "Num. samples" ) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1)) +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
    scale_fill_manual(values = c("#bfbfbf","#666666"),
                      breaks = c("control", "AD"),
                      label = c("Control", "AD")) +
    custom_ggtheme
  
  plot_num_samples
  
  
  ## Gender
  plot_gender <- ggplot(df_metadata %>% 
                          mutate(gender = gender %>% as.character()) %>%
                          dplyr::count(cluster, gender)) +
    geom_bar(aes(y = cluster, x = n, group = gender, fill = gender), alpha = 0.8,
             stat = "identity", position = "dodge") + 
    #facet_grid(~cluster)+
    theme_light() +
    labs(y = "", x = "Num. samples" ) +
    scale_fill_manual(labels = c("Male", "Female"), values = c("1", "2"), palette=scales::hue_pal()) +
    guides(fill = guide_legend(title = "Gender: ", ncol = 2, nrow = 1)) +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
    custom_ggtheme
  
  plot_gender
  
  if (all(is.na(df_metadata$rin))) {
    ## AGE
    plot_age <- ggplot(df_metadata ) +
      geom_density(aes(x = age, fill = cluster), alpha = 0.7) + 
      theme_light() +
      labs(x = "AGE" ) +
      scale_fill_hue() +
      guides(fill = guide_legend(title = NULL, ncol = 3, nrow = 1)) +
      scale_fill_manual(values = c("#bfbfbf","#666666"),
                        breaks = c("control", "AD"),
                        label = c("Control", "AD")) +  
      custom_ggtheme
  } else {
    ## RIN
    plot_rin <- ggplot(df_metadata ) +
      geom_density(aes(x = rin, fill = cluster), alpha = 0.7) + 
      theme_light() +
      labs(x = "RIN" ) +
      scale_fill_hue() +
      guides(fill = guide_legend(title = NULL, ncol = 3, nrow = 1)) +
      scale_fill_manual(values = c("#bfbfbf","#666666"),
                        breaks = c("control", "AD"),
                        label = c("Control", "AD")) +  
      custom_ggtheme
  }
  
  
  
  ## Mapped read depth
  plot_read_depth <- ggplot(df_metadata ) +
    geom_density(aes(x = all_mapped_reads, fill = cluster), alpha = 0.6) + 
    theme_light() +
    scale_fill_hue() +
    labs(x = "Mapped Read Count" ) +
    guides(fill = guide_legend(title = NULL, ncol = 3, nrow = 1)) +
    scale_fill_manual(values = c("#bfbfbf","#666666"),
                      breaks = c("control", "AD"),
                      label = c("Control", "AD")) +   
    custom_ggtheme
  
  
  plot_read_depth
  
  if (!exists("plot_rin")) {
    ggpubr::ggarrange(plot_num_samples,
                      plot_gender,
                      plot_age,
                      plot_read_depth,
                      labels = c("a","b","c","d"))
  } else {
    ggpubr::ggarrange(plot_num_samples,
                      plot_gender,
                      plot_rin,
                      plot_read_depth,
                      labels = c("a","b","c","d"))
  }
 
  
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
  
  
  clusters <- df_metadata$cluster %>% unique
  
  ###############################
  ## GET DATA FOR CONTROL
  ###############################
  
  cluster_id <- "Control"
  
  print(paste0(Sys.time(), " - ", project_name, " - ", cluster_id))
  
  query <- paste0("SELECT DISTINCT ref_junID, MSR_D, MSR_A, ref_type, ref_sum_counts, ref_n_individuals 
                  FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
  introns <- dbGetQuery(con, query) %>% as_tibble()
  introns %>% nrow()
  
  query <- paste0("SELECT DISTINCT ref_junID, MSR_D,  MSR_A, ref_type, ref_sum_counts, ref_n_individuals 
                  FROM '", cluster_id, "_", project_id, "_misspliced'")
  introns_control <- dbGetQuery(con, query) %>% as_tibble()
  introns_control %>% nrow()
  
  introns_control <- rbind(introns, introns_control)
 
  
  ###############################
  ## GET DATA FOR CASES
  ###############################
  
  if (project_id == "SRP100948") {
    cluster_id <- "AD"  
  } else {
    cluster_id <- "PD"
  }
 
  
  print(paste0(Sys.time(), " - ", project_name, " - ", cluster_id))
  
  query <- paste0("SELECT DISTINCT ref_junID, MSR_D, MSR_A, ref_type, ref_sum_counts, ref_n_individuals
                  FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
  introns <- dbGetQuery(con, query) %>% as_tibble()
  introns %>% nrow()
  
  query <- paste0("SELECT DISTINCT ref_junID, MSR_D,  MSR_A, ref_type, ref_sum_counts, ref_n_individuals 
                  FROM '", cluster_id, "_", project_id, "_misspliced'")
  introns_case <- dbGetQuery(con, query) %>% as_tibble()
  introns_case %>% nrow()
  
  introns_case <- rbind(introns, introns_case)
  
  
  ##############################################################################
  ## TIDY DATAFRAME - ONLY USE COMMON INTRONS
  ##############################################################################
  
  common_introns <- intersect(introns_case$ref_junID,
                              introns_control$ref_junID) %>% unique()
  common_introns %>% length()
  
  rbind(introns_case %>%
          mutate(sample_type = "AD"), 
        introns_control %>%
          mutate(sample_type = "control"))  %>%
    filter(ref_junID %in% common_introns) %>%
    dplyr::count(sample_type)
  
  ## TIDY DATAFRAME
  df_introns_tidy <- rbind(introns_case %>%
                             mutate(sample_type = "AD"), 
                           introns_control %>%
                             mutate(sample_type = "control"))  %>%
    filter(ref_junID %in% common_introns) %>%
    group_by(sample_type) %>%
    mutate(mean_coverage = (ref_sum_counts/ref_n_individuals) %>% log10()) %>%
    ungroup()
  
  df_introns_tidy
  
  
  df_introns_tidy %>%
    dplyr::count(sample_type)
  
  
  ##########################################################
  ## PLOT COVERAGE OF SHARED INTRONS BEFORE SUBSAMPLING
  ##########################################################
  
  
  plot_BS <- ggplot(data = df_introns_tidy) +
    geom_density(mapping = aes(x = mean_coverage, fill = sample_type), alpha = 0.8) +
    ggtitle("Before subsampling") +
    theme_light() +
    custom_ggtheme +
    xlab("log10 mean expression level") +
    scale_fill_manual(name = "Sample group: ",
                      values = c("#bfbfbf","#666666"),
                      breaks = c("control", "AD"),
                      label = c("Control", "AD")) +
    ylim(c(0, 1))+
    xlim(c(0, 3))
  
  plot_BS
  
  
  #####################################
  #### START SUBSAMPLIING
  #####################################
  
  if ( file.exists(paste0(results_folder, "/common_subsampled_introns.rds")) ) {
    
    subsample <- readRDS(file = paste0(results_folder, "/common_subsampled_introns.rds"))
    
  } else {
    
    set.seed(1000)
    print(paste0(Sys.time(), " - starting subsampling ... "))
    
    ## Subsampling introns to control by similarity in mean read coverage
    m.out <- MatchIt::matchit(sample_type ~ mean_coverage,
                              data = df_introns_tidy %>%
                                mutate(sample_type = sample_type %>% as.factor()),
                              distance = df_introns_tidy$mean_coverage,
                              method = "nearest",
                              caliper = c(mean_coverage = 0.005),
                              std.caliper = FALSE)
    subsample <- MatchIt::match.data(m.out)
    subsample %>% dplyr::count(sample_type)
    
    saveRDS(object = subsample,
            file = paste0(results_folder, "/common_subsampled_introns_seed1000.rds"))
    
  }
  
  subsample %>%
    group_by(sample_type) %>%
    distinct(ref_junID) %>%
    dplyr::count(sample_type)
  
  
  #######################################
  ## PLOT COVERAGE AFTER SUBSAMPLE
  #######################################
  
  plot_AS <- ggplot(data = subsample) +
    geom_density(mapping = aes(x = mean_coverage, fill = sample_type), alpha = 0.8) +
    ggtitle("After subsampling") +
    theme_light() +
    custom_ggtheme +
    xlab("log10 mean expression level") +
    scale_fill_manual(name = "Sample group: ",
                      values = c("#bfbfbf","#666666"),
                      breaks = c("control", "AD"),
                      label = c("Control", "AD")) +
    ylim(c(0, 1))+
    xlim(c(0, 3))
  plot_AS
  
  ggpubr::ggarrange(plot_BS,
                    plot_AS,
                    labels = c("a", "b"),
                    common.legend = T) 
  

  dir.create(path = figures_path, recursive = T,showWarnings = T)
  ggplot2::ggsave(paste0(figures_path, "/subsampling.png"), width = 180, height = 100, units = "mm", dpi = 300)
  
  


}


#' Title
#' Visualise differences in number of unique novel junctions from common introns from AD vs control samples
#' Only using common, subsampled introns
#' @return
#' @export
#'
#' @examples
plot_unique_junctions <- function() {
  
  ## Load common introns
  common_introns_subsample <- readRDS(file = paste0(results_folder, "/common_subsampled_introns_seed1000.rds")) %>%
    dplyr::select(-distance,-weights,-subclass)
  
  common_introns_subsample$ref_junID %>% unique() %>% length()
  common_introns_subsample %>% dplyr::count(sample_type)
  
  df_unique_junctions <- map_df(all_clusters, function(cluster_id) {
    
    print(cluster_id)
    
    ####################
    ## GET THE INTRONS
    ####################
    
    query <- paste0("SELECT DISTINCT ref_junID FROM '", 
                    cluster_id, "_", project_id, "_nevermisspliced'")
    introns <- dbGetQuery(con, query) %>% as_tibble()
    
    
    query <- paste0("SELECT DISTINCT ref_junID FROM '", 
                    cluster_id, "_", project_id, "_misspliced'")
    introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
    
    ## Only common, subsampled introns
    introns <- introns %>%
      filter(ref_junID %in% (common_introns_subsample %>%
                               filter(sample_type == cluster_id) %>%
                               pull(ref_junID)))
    
    
    
    ############################
    ## GET THE NOVEL JUNCTIONS
    ############################
    
    query <- paste0("SELECT *
                    FROM '", cluster_id, "_", project_id, "_misspliced'")
    novel_junctions <- dbGetQuery(con, query) %>% as_tibble() 
    
    

    novel_junctions <- novel_junctions %>%
      inner_join(y = master_novel_junctions %>% 
                   dplyr::select(novel_junID,novel_type),
                by = "novel_junID") %>% 
      as_tibble() %>%
      ## Only common, subsampled introns
      filter(ref_junID %in%  (common_introns_subsample %>%
                                filter(sample_type == cluster_id) %>%
                                pull(ref_junID)))
    
    
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
  
  df_unique_junctions
  
  ######################
  ## BAR PLOT
  ######################
  
  
  df_unique_junctions_tidy <- df_unique_junctions %>%
    dplyr::select(cluster, 
                  donor = donor_prop, 
                  acceptor = acceptor_prop, 
                  annotated_intron = annotated_prop) %>%
    tidyr::gather(key = "type", value = "prop", -cluster ) 
  
  df_unique_junctions_tidy <- df_unique_junctions_tidy %>%
    mutate(prop = prop * 100)
  
  
  ggplot(data = df_unique_junctions_tidy %>%
           mutate(type = factor(type, levels = c("donor", "acceptor", "annotated_intron")),
                  prop_label = paste0(round(x = prop, digits = 2), "%")),
         aes(x = cluster, y = prop, fill =  type)) + 
    geom_bar(stat="identity", position = "stack")+
    geom_text(aes(label=prop_label), vjust=2, position = "stack", 
              color="white", size=2.5, fontface = "bold") +
    theme_light() +
    ylab("% unique junctions") +
    xlab("") +
    theme(axis.line = element_line(colour = "black"), 
          #text = element_text(colour = "black", size = 12),
          legend.position = "top") +
    scale_fill_manual(values = c( "#999999","#64037d", "#35B779FF"),
                      breaks = c("annotated_intron","acceptor", "donor" ),
                      label = c("Annotated","Acceptor", "Donor" )) + 
    guides(fill = guide_legend(title = "", ncol = 3, nrow = 1 )) +
    custom_ggtheme  + 
    theme( plot.margin = margin(t = 5,r = 5,l = 5, b = -5),
           legend.box.margin=margin(b = -10, t = -5))#+
  
  
  ggplot2::ggsave(paste0(figures_path, "/jxn_unique_junctions.png"), 
                  width = 80, height = 70, units = "mm", dpi = 300)
  
  
  
}


#' Title
#' Visualise differences in the cummulative number of novel reads in common introns from AD vs control samples
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
  common_introns_subsample %>% dplyr::count(sample_type)
  
  df_mean_counts <- map_df(all_clusters, function(cluster_id) {
    
    print(cluster_id)
    
    ####################
    ## GET EXPRESSION LEVELS FROM THE ANNOTATED INTRONS
    ####################
    
    query <- paste0("SELECT DISTINCT ref_junID, ref_sum_counts, ref_n_individuals 
                    FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
    introns <- dbGetQuery(con, query) %>% as_tibble()
    
    
    
    query <- paste0("SELECT DISTINCT ref_junID, ref_sum_counts, ref_n_individuals 
                    FROM '", cluster_id, "_", project_id, "_misspliced'")
    introns <- rbind(introns, dbGetQuery(con, query) %>% 
                       as_tibble())

    ## Only common, subsampled introns
    introns <- introns %>%
      filter(ref_junID %in% (common_introns_subsample %>%
                               filter(sample_type == cluster_id) %>%
                               pull(ref_junID)))
    
    ###########################
    ## GET THE NOVEL JUNCTIONS
    ###########################
    
    query <- paste0("SELECT ref_junID, novel_junID, novel_sum_counts, novel_n_individuals 
                    FROM '", cluster_id, "_", project_id, "_misspliced'")
    novel_junctions <- dbGetQuery(con, query) %>% as_tibble() 
    
    
    novel_junctions <- novel_junctions %>%
      inner_join(y = master_novel_junctions %>% dplyr::select(novel_junID, novel_type) %>% as_tibble(),
                 by = "novel_junID") %>% 
      as_tibble() %>%
      ## Only NOVEL READS FROM common, subsampled introns
      filter(ref_junID %in% (common_introns_subsample %>%
                               filter(sample_type == cluster_id) %>%
                               pull(ref_junID)))
    
    ###########################
    ## GET THE PROPORTIONS
    ###########################
    
    
    annotated <- introns %>%
      dplyr::distinct(ref_junID, .keep_all = T) %>%
      pull(ref_sum_counts) %>% 
      sum()
      #rowwise() %>%
      #mutate(exp_levels = ref_sum_counts/ref_n_individuals) %>%
      #pull(exp_levels) %>%
      #sum()
    
    acceptor <- novel_junctions %>%
      filter(novel_type == "novel_acceptor") %>%
      dplyr::distinct(novel_junID, .keep_all = T) %>%
      pull(novel_sum_counts) %>% 
      sum()
      #rowwise() %>%
      #mutate(exp_levels = novel_sum_counts/novel_n_individuals) %>%
      #pull(exp_levels) %>%
      #sum()
    
    donor <- novel_junctions %>%
      filter(novel_type == "novel_donor") %>%
      dplyr::distinct(novel_junID, .keep_all = T) %>%
      pull(novel_sum_counts) %>% 
      sum()
      #rowwise() %>%
      #mutate(exp_levels = novel_sum_counts/novel_n_individuals) %>%
      #pull(exp_levels) %>%
      #sum()
    
    annotated_p = annotated * 100 / (annotated + acceptor + donor)
    acceptor_p = acceptor * 100 / (annotated + acceptor + donor)
    donor_p = donor * 100 / (annotated + acceptor + donor)
    
    return(data.frame(cluster = cluster_id,
                      type = c("annotated","acceptor", "donor"),
                      prop = c(annotated_p, acceptor_p, donor_p),
                      sum_counts = c(annotated, acceptor, donor)))
    
  })
  
  df_mean_counts

  ######################
  ## BAR PLOT
  ######################
  
  p <- ggplot(df_mean_counts %>%
           mutate(type = factor(type, levels = c("donor", "acceptor", "annotated")),
                  prop_label = paste0(round(x = prop, digits = 2), "%")),
         aes(x = cluster, y = prop, fill =  type)) + 
    geom_bar(stat = "identity") +
    ggforce::facet_zoom(ylim = c(95,100)) +
    geom_text(aes(label=prop_label), 
              vjust=1.3, 
              position = "stack", 
              color="white", size=2.5, fontface = "bold") +
    theme_light() +
    ylab("% cumulative split read counts") +
    xlab("") +   
    custom_ggtheme  +
    theme(axis.line = element_line(colour = "black"), 
          legend.position = "top")  +
    scale_fill_manual(values = c( "#999999","#64037d", "#35B779FF"),
                      breaks = c( "annotated","acceptor", "donor"),
                      label = c("Annotated","Acceptor", "Donor")) + 
    guides(fill = guide_legend(title = "", ncol = 3, nrow = 1 )) + 
    theme( plot.margin = margin(t = 5,r = 5,l = 5,
                                b = -5),
           legend.box.margin=margin(b = -10, t = -5))#+
  
  ## Only annotated zoomed graph
  pb <- ggplot_build(p)
  pb$data[[2]][1:6, 'alpha'] <- 0
  pg <- ggplot_gtable(pb)
  
  
  ## save plot
  
  png(filename=paste0(figures_path, "/jxn_cummulative_read_counts.png"), width = 100, height = 70,  units = "mm", res = 300)
  plot(pg)
  dev.off()

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
                  FROM '", cluster_id, "_", project_id, "_misspliced' AS tissue 
                  INNER JOIN 'novel' ON novel.novel_junID = tissue.novel_junID")
    df_misspliced <- dbGetQuery(con, query) %>% 
      as_tibble() %>%
      filter(ref_junID %in% (common_introns_subsample %>%
                               filter(sample_type == cluster_id) %>%
                               pull(ref_junID))) %>%
      mutate(cluster = cluster_id)
  })
  

  ################################
  ## PLOT DISTANCES
  ################################
  
  novel.labs <- c("Novel Donor", "Novel Acceptor")
  names(novel.labs) <- c("novel_donor", "novel_acceptor")
  
  plot_distances <- ggplot(data = df_misspliced %>%
                             mutate(type_colour = paste0(cluster, "_", novel_type)) %>%
                             filter(abs(distance) <= limit_bp) %>%
                             mutate(novel_type = factor(novel_type,
                                                        levels = c("novel_donor", "novel_acceptor"),
                                                        labels = c("Novel Donor", "Novel Acceptor")
                             ))) + 
    geom_histogram(aes(x = distance, fill = type_colour),
                   bins = limit_bp * 2,
                   binwidth = 1,
                   position = "stack") +
    xlab("Distance (in bp)") +
    ylab("Unique novel junctions") +
    facet_grid(novel_type~fct_rev(cluster)) +
    theme_light() +
    scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
                       breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6)) +
    scale_fill_manual(values = c("#9ce2c0","#35B779FF","#e99cfc", "#8c04ae"),
                      breaks = c("control_novel_donor", "AD_novel_donor","control_novel_acceptor", "AD_novel_acceptor"),
                      label = c("Donor - Control", "Donor - AD","Acceptor - Control", "Acceptor - AD")) + 
    guides(fill = guide_legend(title = NULL, ncol = 4, nrow = 1 )) +
    custom_ggtheme +
    theme(legend.position = "none") 
  plot_distances
  
  
  distance_rectangle <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 1, ymax = 100),
              fill = "grey", color = "black") +
    geom_text(aes(x = 15, y = 55),  size = 3, label = "exon") +
    geom_rect(aes(xmin = (limit_bp)*-1, xmax = 0, ymin = 49, ymax = 51),
              fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -15, y = 70),  size = 3, label = "intron") +
    theme_void()
  
  
  plot_distances <- plot_distances / (distance_rectangle + distance_rectangle ) + 
    patchwork::plot_layout(heights = c(8, 1))
  plot_distances
  
  ggplot2::ggsave(paste0(figures_path, "/distances.png"), 
                  width = 180, height = 100, units = "mm", dpi = 300)
  
}


#' Title
#' Visualise differences in modulo3 values in introns from AD vs control samples
#' Only using common, subsampled introns
#' @return
#' @export
#'
#' @examples
plot_modulo <- function() {
  
  common_introns_subsample <- readRDS(file = paste0(results_folder, "/common_subsampled_introns.rds"))
  
  common_introns_subsample %>%
    dplyr::count(sample_type)
  
  common_introns_subsample <- common_introns_subsample %>% distinct(ref_junID)
  
  ## Calculate modulo3 from MANE transcripts in 100bp distance
  df_modulo <- map_df(all_clusters, function(cluster_id) {
      
    # cluster_id <- all_clusters[1]
    print(paste0(Sys.time(), " - ", cluster_id))
    
    #####################################
    ## GET NOVEL JUNCTIONS
    #####################################
    
    query <- paste0("SELECT novel_junID 
                    FROM '", cluster_id, "_", project_id, "_misspliced'")
    novel_jnx <- dbGetQuery(con, query) %>% as_tibble()
      
    
    query <- paste0("SELECT ref_junID, novel_junID, distance
                    FROM 'novel' WHERE novel_junID IN (",
                    paste(novel_jnx$novel_junID, collapse = ","),")")
    novel_jnx <- novel_jnx %>%
      left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                by = "novel_junID") %>% 
      as_tibble() 
    
      
    #####################################
    ## GET INTRON AND TRANSCRIPT INFO
    #####################################
    
    ## Add the transcript and MANE info
    query <- paste0("SELECT intron.ref_junID, intron.protein_coding, transcript.MANE
                    FROM 'intron' 
                    INNER JOIN transcript
                    ON intron.transcript_id = transcript.id
                    WHERE ref_junID IN (", paste(novel_jnx$ref_junID, collapse = ","),")")
    all_jxn <- novel_jnx %>%
      inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
                by = "ref_junID") %>% 
      as_tibble() 
    
    
    df_novel_tidy <- all_jxn %>%
      distinct(novel_junID, .keep_all = T) %>%
      filter(abs(distance) <= 100, 
             MANE == 1,
             ref_junID %in% (common_introns_subsample$ref_junID)) %>% 
      mutate(modulo = abs(distance) %% 3) 
      
    df_novel_tidy %>% distinct(ref_junID) %>% nrow %>% print
    
    df_novel_tidy <- df_novel_tidy %>% 
      group_by(modulo) %>%
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>%
      mutate(cluster = cluster_id)
    
    return(df_novel_tidy)
    
  })
  
  df_modulo
  
  df_modulo$modulo = factor(df_modulo$modulo, 
                            levels = c( "0", "1", "2"))
  
  
  ################
  ## BAR PLOT
  ################

  df_modulo_tidy <- df_modulo %>%
    mutate(freq = freq * 100)
  
  
  ggplot(data = df_modulo_tidy %>%
           mutate(freq_label = paste0(round(x = freq, digits = 2), "%"),
                  modulo = factor(modulo, levels = c("2", "1", "0"))),
         aes(x = cluster, y = freq, fill = modulo)) +
    geom_bar(stat = "identity") +
    xlab("") +
    ylab("% of junctions") +
    geom_text(aes(label=freq_label), vjust=2, position = "stack", 
              color="white", size=2.5, fontface = "bold") +
    ggsci::scale_fill_npg(breaks=c('0', '1', '2')) +
    theme_light() +
    custom_ggtheme +
    guides(fill = guide_legend(title = "Modulo: ", ncol = 3, nrow = 1 )) +
    theme( plot.margin = margin(t = 5, r = 5, l = 5, b = -3),
           legend.box.margin=margin(b = -10, t = -5)) 
  
  
  ggplot2::ggsave(paste0(figures_path, "/modulo.png"), 
                  width = 60, height = 70, units = "mm", dpi = 300)
  
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
  common_junID <- common_introns_subsample %>%
    dplyr::count(ref_junID) %>%
    filter(n == 2) %>%
    pull(ref_junID)
  
  
  common_introns_subsample_tidy <- common_introns_subsample %>%
    filter(ref_junID %in% common_junID)
  
  common_introns_subsample_tidy_MSRD <- common_introns_subsample_tidy %>%
    dplyr::select(ref_junID, MSR_D, sample_type) %>%
    spread(sample_type, value = MSR_D)
  
  common_introns_subsample_tidy_MSRA <- common_introns_subsample_tidy %>%
    dplyr::select(ref_junID, MSR_A, sample_type) %>%
    spread(sample_type, value = MSR_A)
  
  #########################################
  ## TEST DONOR - AFTER SUBSAMPLING
  #########################################
  
  
  common_introns_subsample_tidy %>% filter(sample_type == "AD") %>% pull(MSR_D) %>% summary
  common_introns_subsample_tidy %>% filter(sample_type == "control") %>% pull(MSR_D) %>% summary
  
  
  wilcox.test(x = common_introns_subsample_tidy %>% filter(sample_type == "AD") %>% pull(MSR_D),
              y = common_introns_subsample_tidy %>% filter(sample_type == "control") %>% pull(MSR_D),
              alternative = "greater",
              paired = T,
              correct = T)
  
  wilcox.test(x = common_introns_subsample_tidy_MSRD %>% pull(AD),
              y = common_introns_subsample_tidy_MSRD %>% pull(control),
              alternative = "greater",
              paired = T,
              correct = T)
  
  rstatix::wilcox_effsize(data = common_introns_subsample_tidy_MSRD %>% 
                            gather(key = ref_junID, value = type) %>%
                            dplyr::rename(type = ref_junID, MSR_D = type) %>%
                            mutate(type = type %>% as.factor()),
                          formula = MSR_D ~ type,
                          paired = T)
  
  ########################################
  ## TEST ACCEPTOR
  ########################################
  
  
  wilcox.test(x = common_introns_subsample_tidy %>% filter(sample_type == "AD") %>% pull(MSR_A),
              y = common_introns_subsample_tidy %>% filter(sample_type == "control") %>% pull(MSR_A),
              alternative = "greater",
              paired = T,
              correct = T)
  
  wilcox.test(x = common_introns_subsample_tidy_MSRA %>% pull(AD),
              y = common_introns_subsample_tidy_MSRA %>% pull(control),
              alternative = "greater",
              paired = T,
              correct = T)
  
  
  common_introns_subsample_tidy %>% filter(sample_type == "AD") %>% pull(MSR_A) %>% summary
  common_introns_subsample_tidy %>% filter(sample_type == "control") %>% pull(MSR_A) %>% summary
  
  rstatix::wilcox_effsize(data = common_introns_subsample_tidy_MSRA %>% 
                            gather(key = ref_junID, value = type) %>%
                            dplyr::rename(type = ref_junID, MSR_A = type) %>%
                            mutate(type = type %>% as.factor()),
                          formula = MSR_A ~ type,
                          paired = T)
  
  
  #######################################
  ## EXACT PAIRING
  #######################################
  
  common_introns_subsample_column <- common_introns_subsample_tidy %>% 
    filter(sample_type == "AD") %>% 
    dplyr::select(ref_junID, MSR_A_AD = MSR_A) %>%
    inner_join(y = common_introns_subsample_tidy %>% 
                 filter(sample_type == "control") %>% 
                 dplyr::select(ref_junID, MSR_A_control = MSR_A) )
  
  wilcox.test(x = common_introns_subsample_column$MSR_A_AD,
              y = common_introns_subsample_column$MSR_A_control,
              paired = T,
              alternative = "greater",
              correct = T)
  
}



plot_MSR_by_biotype <- function() {
  
  common_introns_subsample <- readRDS(file = paste0(results_folder, "/common_subsampled_introns.rds"))
  common_junID <- common_introns_subsample %>%
    dplyr::count(ref_junID) %>%
    filter(n == 2) %>%
    pull(ref_junID)
  
  
  ################################
  ## TIDY DATA BEFORE PLOTTING
  ################################
  
  df_introns_biotype <- common_introns_subsample %>%
    filter(ref_junID %in% common_junID) %>%
    inner_join(y = master_introns %>% dplyr::select(ref_junID, protein_coding),
               by = "ref_junID") %>%
    filter(protein_coding %in% c(0,100)) %>%
    mutate(biotype = ifelse(protein_coding == 100, "PC", "non PC"))
  


  df_introns_biotype$biotype = factor(df_introns_biotype$biotype, 
                                      levels = c("non PC","PC"))
  
  df_introns_biotype <- df_introns_biotype %>%
    dplyr::select(ref_junID, biotype, MSR_D, MSR_A) %>%
    gather(key = "MSR_type", value = "MSR", -biotype, -ref_junID)
  
  
  df_introns_biotype$MSR_type = factor(df_introns_biotype$MSR_type, 
                                       levels = c("MSR_A", "MSR_D"))
  
  
  print(paste0(Sys.time(), " - ", df_introns_biotype %>% nrow(), " - introns after tidying!"))
  
  
  df_introns_biotype %>% 
    filter(biotype == "PC") %>%
    pull(MSR) %>%
    summary()
  
  df_introns_biotype %>% 
    filter(biotype == "non PC") %>%
    pull(MSR) %>%
    summary()
  
  
  ###########################
  ## PLOT BY MSR GROUP
  ###########################
  
  # scales::show_col(scales::viridis_pal(option = "B")(40))
  # scales::show_col(scales::viridis_pal(option = "C")(40))
  # scales::show_col(scales::viridis_pal(option = "D")(40))
  # scales::show_col(scales::viridis_pal(option = "E")(40))
  # scales::show_col(scales::viridis_pal(option = "F")(40))
  # scales::show_col(scales::viridis_pal(option = "G")(40))
  # scales::show_col(scales::viridis_pal(option = "H")(40))
  
  
  df_introns_biotype_tidy <- df_introns_biotype %>%
    group_by(biotype, MSR_type) %>%
    distinct(ref_junID, .keep_all = T) %>%
    mutate(percentile_group = case_when(MSR == 0 ~ "0",
                                        MSR > 0 & MSR <= 0.2 ~ "(0,.2]",
                                        MSR > 0.2 & MSR <= 0.4 ~ "(.2,.4]",
                                        MSR > 0.4 & MSR <= 0.6 ~ "(.4,.6]",
                                        MSR > 0.6 & MSR <= 0.8 ~ "(.6,.8]",
                                        MSR > 0.8 & MSR <= 1 ~ "(.8,1]"))
  
  df_introns_biotype_tidy$MSR_type = factor(df_introns_biotype_tidy$MSR_type, 
                                            levels = c("MSR_D","MSR_A"))
  df_introns_biotype_tidy$percentile_group = factor(df_introns_biotype_tidy$percentile_group, 
                                                    levels = c("0","(0,.2]","(.2,.4]","(.4,.6]","(.6,.8]","(.8,1]"))
  
  df_introns_biotype_tidy$biotype = factor(df_introns_biotype_tidy$biotype, 
                                           levels = c("PC","non PC"))
  
  plotMSR_donor <- ggplot(data = df_introns_biotype_tidy ) + 
    geom_bar(aes(x = percentile_group, fill = biotype),
             position = "dodge", linewidth = .5, color = "#333333")+
    #ggtitle("MSR Donor") +
    xlab("Mis-splicing ratio value group") +
    ylab("Number of annotated introns") +
    #ylim(c(0,23000))+
    ggforce::facet_row(~MSR_type) +
    scale_fill_manual(values = c("#333333","#999999"),
                      breaks = c("PC","non PC"),
                      labels = c("Protein-coding   ","Non-protein-coding ")) +
    theme_light() +
    custom_ggtheme +
    guides(fill = guide_legend(title = "",
                               ncol = 2, nrow = 1))
  plotMSR_donor
  
  
  ## Non-protein-coding
  plotMSR_acceptor <- ggplot(data = df_introns_biotype_tidy %>% 
                               filter(MSR_type == "MSR_A")) + 
    geom_bar(aes(x = percentile_group, fill = biotype),
             position = "dodge", linewidth = 0.5, color = "#333333")+
    ggtitle("MSR Acceptor") +
    xlab("Mis-splicing ratio value group") +
    ylab("Number of annotated introns") +
    #ylim(c(0,23000))+
    
    scale_fill_manual(values = c("#333333","#999999"),
                      breaks = c("PC","non PC"),
                      labels = c("Protein-coding ","Non-protein-coding ")) +
    theme_light() +
    custom_ggtheme +
    guides(fill = guide_legend(ncol = 2, nrow = 1))
  
  
  ggpubr::ggarrange(plotMSR_donor,
                    plotMSR_acceptor,
                    nrow = 1,
                    ncol = 2,
                    common.legend = T)
  
  
  
  
  print(paste0(Sys.time(), " - saving plot..."))
  file_name <- paste0(figures_folder, "/panel4bc")
  ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 90, units = "mm", dpi = 300)
}


#' Title
#' GO, KEGG and REACTOME ENRICHMENT analysis of the introns showing increasing MSR values in AD compared to control samples
#' Only using common, subsampled introns
#' @return
#' @export
#'
#' @examples
plot_GO_KEGG__REACTOME_terms <- function() {
  
  ## Load common introns
  common_introns_subsample <- readRDS(file = paste0(results_folder, "/common_subsampled_introns.rds")) %>%
    dplyr::select(-distance,-weights,-subclass)
  
  common_junID <- common_introns_subsample %>%
    dplyr::count(ref_junID) %>%
    filter(n == 2) %>%
    pull(ref_junID)
  
  common_introns_subsample_tidy <- common_introns_subsample %>%
    filter(ref_junID %in% common_junID)
  
  
  common_introns_subsample_column <- common_introns_subsample_tidy %>% 
    filter(sample_type == "AD") %>% 
    dplyr::select(ref_junID, MSR_D_AD = MSR_D, MSR_A_AD = MSR_A) %>% 
    distinct(ref_junID, .keep_all = T) %>%
    inner_join(y = common_introns_subsample_tidy %>% 
                 filter(sample_type == "control")  %>% 
                 dplyr::select(ref_junID, MSR_D_control = MSR_D, MSR_A_control = MSR_A) %>% 
                 distinct(ref_junID, .keep_all = T),
               by = "ref_junID")
  
  
  ## GET MASTER INTRON DATA
  df_intron <- master_introns %>% 
    dplyr::select(ref_junID, transcript_id) %>%
    inner_join(y = master_transcripts %>% dplyr::select(id, gene_id),
               by =c("transcript_id" = "id")) %>%
    inner_join(y = master_gene %>% dplyr::select(id, gene_name, gene_id),
               by =c("gene_id" = "id")) 
  
  
  ## Add gene name to introns with increasing MSR_A values in AD compared to control
  introns_increasing_msr <- common_introns_subsample_column %>%
    filter((MSR_A_AD > MSR_A_control) | (MSR_D_AD > MSR_D_control)) %>%
    inner_join(y = df_intron,
               by = c("ref_junID")) %>%
    distinct(ref_junID, .keep_all = T)
  
  
  introns_increasing_msr
  # introns_descending_msra <- common_introns_subsample_column %>%
  #   filter(MSR_A_AD < MSR_A_control) %>%
  #   inner_join(y = df_intron,
  #              by =c("ref_junID")) 
  
  
  ## Get gene background data
  bg_genes <- common_introns_subsample_column %>%
    inner_join(y = df_intron,
               by =c("ref_junID")) 
  
  ################################
  ## GO ENRICHMENT
  ################################
  
  ego_MSRA <- clusterProfiler::enrichGO(
    gene          = introns_increasing_msr$gene_name %>% unique,
    universe      = bg_genes$gene_name %>% unique,
    keyType       = "SYMBOL",
    OrgDb         = "org.Hs.eg.db", ##Genome wide annotation for Human, primarily based on mapping using Entrez Gene identifiers.
    ont           = "ALL",
    pAdjustMethod = "fdr",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  
  category_terms <-25
  clusterProfiler::dotplot(ego_MSRA %>%
                             filter(ONTOLOGY == "CC") %>%
                             filter(str_detect(string = Description,negate = T,
                                               pattern = "with bulged")) %>%
                             filter(str_detect(string = Description,negate = T,
                                               pattern = "ubiquitin-dependent")), 
                           x = "GeneRatio", 
                           showCategory = category_terms, 
                           split="ONTOLOGY") +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 60)) +
    xlab("Gene Ratio") +
    ggforce::facet_col(ONTOLOGY~., scales = "free_y", space = "free") +
    custom_ggtheme+
    theme(#axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = "top",
          legend.box="horizontal",
          plot.margin = margin(t = -5,b = 0,l = 2,r = 2),
          legend.margin=margin(t = -5,b = -5,0,0),
          legend.box.margin=margin(t = -5,r = -5,b = -5,l = -5)) + 
    scale_size(range = c(1, 5)) +
    guides(size = guide_legend(title = "Gene Count: "),
           colour = guide_legend(title = "q: "))               
  
  
  ggplot2::ggsave(paste0(figures_path, "/AD_Control_go_enrichment_MSR_",category_terms,".png"), 
                  width = 180, height = 75, units = "mm", dpi = 300)
  
 
  
  ################################
  ## KEGG ENRICHMENT
  ################################
  
  library('org.Hs.eg.db')
  
  
  category_terms <-20

  
  # mapIds(org.Hs.eg.db, introns_increasing_msra, 'ENTREZID', 'SYMBOL')
  
  ekegg_MSR <- clusterProfiler::enrichKEGG(
    gene          = mapIds(x = org.Hs.eg.db, 
                           keys = introns_increasing_msr %>% drop_na %>% pull(gene_name) %>% unique, 
                           column = 'ENTREZID', 
                           keytype = 'SYMBOL')%>% unlist %>% unname ,
    organism      = "hsa",
    keyType       = "kegg",
    universe      = mapIds(x = org.Hs.eg.db, 
                           keys = bg_genes %>% drop_na %>% pull(gene_name) %>% unique(), 
                           column = 'ENTREZID', 
                           keytype = 'SYMBOL') %>% unlist %>% unname , 
    pAdjustMethod = "fdr")
  
  
  plotKEGG <-  clusterProfiler::dotplot(ekegg_MSR %>% 
                                          mutate(ONTOLOGY = "KEGG") %>%
                                          filter(!(str_detect(string = Description, pattern = "Chemical carcino"))) %>%
                                          filter(!(str_detect(string = Description, pattern = "Protein processing in"))) %>%
                                          filter(!(str_detect(string = Description, pattern = "Amino sugar"))), 
                                        showCategory = category_terms, 
                                        split="ONTOLOGY")
  
  plotKEGG +
    scale_y_discrete(labels = 
                       function(x) stringr::str_wrap(x, width = 40)) +
    xlab("Gene Ratio") +
    ggforce::facet_row(ONTOLOGY~., scales = "free", space = "free") +
    #Coord_flip() +
    custom_ggtheme +
    theme(text = element_text(colour = "black", size = 7),
          legend.position = "top",
          legend.box="vertical",
          plot.margin = margin(t = 0, b =0, 0, 0),
          legend.margin = margin(t = -5, b = -5, r = 1,l = 1),
          legend.box.margin=margin(t = -10,b = -5, r = 0,l = 0)) + 
    scale_size(range = c(1, 5))+
    guides(colour = guide_legend(title = "q: "),
           size = guide_legend(title = "Gene count: "),
           nrow = 2, ncol = 2 ) 
  
  ggplot2::ggsave(paste0(figures_path, "/AD_Control_kegg_enrichment_MSR_",category_terms,".png"), 
                  width = 120, height = 70, units = "mm", dpi = 300)
  
  ################################
  ## REACTOME ENRICHMENT
  ################################
  
  
  library(ReactomePA)
  
  reactome_MSR <- ReactomePA::enrichPathway(gene = mapIds(x = org.Hs.eg.db, 
                                                          keys = introns_increasing_msr$gene_name %>% unique,
                                                          column = 'ENTREZID', 
                                                          keytype = 'SYMBOL') %>% unlist %>% unname, 
                                         pvalueCutoff = 0.05,
                                         pAdjustMethod = "BH",
                                         organism = "human",
                                         universe = mapIds(x = org.Hs.eg.db, 
                                                           keys = bg_genes$gene_name %>% unique, 
                                                           column = 'ENTREZID',
                                                           keytype = 'SYMBOL'),
                                         readable=TRUE)
  
  
  plotREAC <-  clusterProfiler::dotplot(reactome_MSR %>% mutate(ONTOLOGY = "REACTOME"), 
                                        showCategory = category_terms, 
                                        split="ONTOLOGY")
  
  plotREAC +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
    xlab("Gene Ratio") +
    ggforce::facet_row(ONTOLOGY~., scales = "free", space = "free") +
    custom_ggtheme +
    theme(legend.position = "top",
          legend.box="horizontal") + 
    scale_size(range = c(1, 5))+
    guides(colour = guide_legend(title = "q: "),
           size = guide_legend(title = "Gene count: "))
  
  ggplot2::ggsave(paste0(figures_path, "/AD_Control_reactome_enrichment_MSR_",category_terms,".png"), 
                  width = 180, height = 180, units = "mm", dpi = 300)
}


##################################
## CALLS
##################################
get_common_subsample_introns
  