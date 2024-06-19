library(tidyverse)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(DBI)
library(dplyr)
library(here)
library(doParallel)
library(Biostrings)
library(tidyverse)
library(protr)
library(optparse)

# source("/home/sruiz/PROJECTS/recount3-database-project/init_age.R")

base_folder <- here::here()
dependencies_folder <- paste0(base_folder, "/dependencies/")

#####################################
## LOAD SOURCE SCRIPTS
#####################################

scripts_folder <- paste0(base_folder, "/scripts/")
setwd(file.path(scripts_folder))
files.sources = list.files()
sapply(files.sources, source)
setwd(file.path(base_folder))


#####################################
## SET MAIN VARIABLES
## AND PATHS
#####################################


gtf_version <- 111

supportive_reads <- 1
replace <- T

data_source <- "data_sources/gtex"


project_name <- paste0("GTEX_", supportive_reads, "read_subsampleFALSE")
database_name <- paste0(project_name, "_age")


levelqc1_folder <- paste0(base_folder, "/database/", project_name, "/", gtf_version, "/")
database_folder <- paste0(base_folder, "/database/", database_name, "/", gtf_version, "/")
results_folder <- paste0(base_folder, "/results/", project_name, "/", gtf_version, "/")
tpm_folder <- paste0(base_folder, "/results/tpm/")

age_projects <- readRDS(file = paste0(results_folder, "/all_final_projects_used.rds"))


##################################
## CALLS - SQL GENERATION
##################################


if ( !file.exists(paste0(results_folder, "/age_init_projects.rds")) ) {
  project_init <- age_stratification_init_data(projects.id = age_projects,
                                               gtf.version = gtf_version,
                                               project.name = project_name,
                                               results.folder = results_folder) %>%
    as_tibble()

  saveRDS(object = project_init,
          file = paste0(results_folder, "/age_init_projects.rds"))
} else {
  project_init <- readRDS(file = paste0(results_folder, "/age_init_projects.rds"))
}

age_projects <- project_init$project %>% unique()

# for (project_id in age_projects) {
# 
#   # project_id <- age_projects[1]
#   # project_id <- "BRAIN"
# 
#   print(paste0(Sys.time(), " --> ", project_id))
#   project_init_local <- project_init %>% filter(project == project_id)
#   project_results_folder <- file.path(results_folder, project_id)
# 
#   age_stratification_annotate(age.groups = project_init_local$age_group %>% unique(),
#                               project.id = project_id,
#                               gtf.version = gtf_version,
#                               project.name = project_name,
#                               age.samples.clusters = project_init_local,
#                               results.folder = project_results_folder)
# 
#   age_stratification_junction_pairing(age.groups = project_init_local$age_group %>% unique(),
#                                       project.id = project_id,
#                                       gtf.version = gtf_version,
#                                       project.name = project_name,
#                                       age.samples.clusters = project_init_local,
#                                       results.folder = project_results_folder,
#                                       replace)
# 
# }


# get_all_annotated_split_reads(recount3.project.IDs = age_projects,
#                               all.clusters = project_init$age_group %>% unique(),
#                               database.folder = database_folder,
#                               results.folder = results_folder,
#                               num.cores = 8)
# 
# 
# 
# get_all_raw_jxn_pairings(recount3.project.IDs = age_projects,
#                          all.clusters = project_init$age_group %>% unique(),
#                          database.folder = database_folder,
#                          results.folder = results_folder,
#                          num.cores = 8)


age_projects <- readRDS(file = paste0(results_folder, "/all_final_projects_used.rds"))


# generate_recount3_tpm(recount3.project.IDs = age_projects,
#                       data.source = data_source,
#                       tpm.folder = tpm_folder,
#                       results.folder = results_folder)



# tidy_data_pior_sql(recount3.project.IDs = age_projects,
#                    levelqc1.folder = levelqc1_folder,
#                    all.clusters = project_init$age_group %>% unique(),
#                    database.folder = database_folder,
#                    results.folder = results_folder)

# if ( !file.exists(paste0(database_folder, "/all_split_reads_qc_level2_PC_biotype.rds")) ) {
#   generate_transcript_biotype_percentage(gtf.version = gtf_version,
#                                          database.folder = paste0(base_folder, "/database/", 
#                                                                   project_name, "/", gtf_version, "/"),
#                                          results.folder = database_folder)
# }

database_path <- paste0(database_folder,  "/", database_name, ".sqlite")

sql_database_generation(database.path = database_path,
                        recount3.project.IDs = age_projects,
                        remove.all = T,
                        database.folder = database_folder,
                        results.folder = results_folder,
                        gtf.version = gtf_version)

