library(tidyverse)
library(GenomicRanges)
library(DBI)

####################################################
## CONNECT TO THE SPLICING DATABASE ################
####################################################


## CONNECT TO THE DATABASE ------------------------------


gtf_version <- 105
main_project <- "ad_control"
setwd(normalizePath(paste0("/home/soniagr/PROJECTS/", main_project, "/database_generation/")))
database_path <- paste0(getwd(), "/database/v", gtf_version, "/", main_project,
                          "/", main_project, ".sqlite")


con <- dbConnect(RSQLite::SQLite(), database_path)
dbListTables(con)

## GET FROM MASTER TABLE
query = paste0("SELECT * FROM 'master'")
df_metadata <- dbGetQuery(con, query) %>% as_tibble()

all_projects <- df_metadata$SRA_project %>% unique
# all_projects %>% length() %>% print()


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



get_MSR <- function()  {
  
  project_id <- "SRP100948" #"SRP058181"
  
  ###############################
  ## GET DATA FOR CONTROL
  ###############################
  
  #query <- paste0("SELECT intron.ref_junID, gene.gene_name FROM 'intron'
  #                INNER JOIN 'transcript' ON transcript.id = intron.transcript_id
  #                INNER JOIN 'gene' ON gene.id = transcript.gene_id")
  #master_introns <- dbGetQuery(con, query) %>% as_tibble()
  
  
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
  
  #project_id <- "SRP058181"
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
  
  common_introns <- intersect(introns_pd$ref_junID,introns_control$ref_junID) %>% unique()
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
  
  
  #########################################
  ## TEST DONOR
  #########################################
  
  subsample_tidy_D <- df_introns_tidy %>% filter(sample_type == "AD") %>% select(ref_junID, MSR_D_AD = MSR_D) %>%
    inner_join(y = df_introns_tidy %>% filter(sample_type == "control") %>% select(ref_junID, MSR_D_control = MSR_D),
               by = "ref_junID")# %>%
  #inner_join(y = master_introns,
  #           by = "ref_junID")
  
  wilcox.test(x = subsample_tidy_D$MSR_D_AD,
              y = subsample_tidy_D$MSR_D_control,
              alternative = "greater",
              correct = T)
  
  subsample_tidy_D$MSR_D_AD %>% summary
  subsample_tidy_D$MSR_D_control %>% summary
  
  
  ########################################
  ## TEST ACCEPTOR
  ########################################
  
  subsample_tidy_A <- df_introns_tidy %>% filter(sample_type == "AD") %>% select(ref_junID, MSR_A_AD = MSR_A) %>%
    inner_join(y = df_introns_tidy %>% filter(sample_type == "control") %>% select(ref_junID, MSR_A_control = MSR_A),
               by = "ref_junID") #%>%
  #inner_join(y = master_introns,
  #           by = "ref_junID")
  
  #subsample_tidy_A %>%
  #  filter(MSR_A_PD > MSR_A_control)
  
  
  wilcox.test(x = subsample_tidy_A$MSR_A_AD,
              y = subsample_tidy_A$MSR_A_control,
              alternative = "greater",
              paired = T,
              correct = T)
  
  
  subsample_tidy_A$MSR_A_AD %>% summary
  subsample_tidy_A$MSR_A_control %>% summary
  
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
  
  saveRDS(object = subsample,
          file = "~/PROJECTS/ad_control/database_generation/results/SRP100948/v105/ad_control/results/MSR_intron_subsample_readcoverage.rds")
  
  #########################################
  ## TEST DONOR - AFTER SUBSAMPLING
  #########################################
  
  subsample_tidy_D <- subsample %>% filter(sample_type == "AD") %>% select(ref_junID, MSR_D_AD = MSR_D) %>%
    inner_join(y = subsample %>% filter(sample_type == "control") %>% select(ref_junID, MSR_D_control = MSR_D),
               by = "ref_junID")# %>%
  #inner_join(y = master_introns,
  #           by = "ref_junID")
  
  wilcox.test(x = subsample_tidy_D$MSR_D_AD,
              y = subsample_tidy_D$MSR_D_control,
              alternative = "greater",
              correct = T)
  
  subsample_tidy_D$MSR_D_AD %>% summary
  subsample_tidy_D$MSR_D_control %>% summary
  
  
  ########################################
  ## TEST ACCEPTOR
  ########################################
  
  subsample_tidy_A <- subsample %>% filter(sample_type == "AD") %>% select(ref_junID, MSR_A_AD = MSR_A) %>%
    inner_join(y = subsample %>% filter(sample_type == "control") %>% select(ref_junID, MSR_A_control = MSR_A),
               by = "ref_junID") #%>%
  #inner_join(y = master_introns,
  #           by = "ref_junID")
  
  #subsample_tidy_A %>%
  #  filter(MSR_A_PD > MSR_A_control)
  
  
  wilcox.test(x = subsample_tidy_A$MSR_A_AD,
              y = subsample_tidy_A$MSR_A_control,
              alternative = "greater",
              paired = T,
              correct = T)
  
  
  subsample_tidy_A$MSR_A_AD %>% summary
  subsample_tidy_A$MSR_A_control %>% summary
  
  master_introns
}


 

##################################
## CALLS
##################################

get_MSR()

# compare_tissues_somatic_mutations(project_id1 = "BRAIN",
#                                   cluster_id1 = "Brain - Nucleus accumbens (basal ganglia)",
#                                   project_id2 = "BLOOD",
#                                   cluster_id2 = "Whole Blood")
  