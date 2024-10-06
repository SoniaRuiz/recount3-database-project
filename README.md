[![DOI](https://zenodo.org/badge/663495574.svg)](https://zenodo.org/doi/10.5281/zenodo.12794619)

# Splicing accuracy varies across human introns, tissues, age and disease

*Sonia Garcia-Ruiz, [David Zhang](https://github.com/dzhang32), [Emil K Gustavsson](https://github.com/egustavsson), Guillermo Rocamora-Perez, [Melissa Grant-Peters](https://github.com/mgrantpeters), Aine Fairbrother-Browne, Regina H Reynolds, Jonathan W Brenton, Ana L Gil-Martinez, Zhongbo Chen, Donald C Rio, Juan A Botia, Sebastian Guelfi, [Leonardo Collado-Torres](https://lcolladotor.github.io/), Mina Ryten*

bioRxiv 2023.03.29.534370;
doi: [https://doi.org/10.1101/2023.03.29.534370](https://doi.org/10.1101/2023.03.29.534370)


# Overview 

The recount3-database-project repository contains the code used to generate all the databased produced for the manuscript [**Splicing accuracy varies across human introns, tissues and age**](https://doi.org/10.1101/2023.03.29.534370).
It contains R scripts designed to generate three different SQLite databases from junction data. The scripts utilize publicly available datasets to facilitate biological data analysis.

## Table of Contents
1. [Installation](#installation)
2. [Usage](#usage)
   - [init.R](#initr)
   - [init_age.R](#init_ager)
   - [init_ENCODE.R](#init_encoder)
3. [Function Descriptions](#function-descriptions)
4. [License](#license)

---

## Installation
To use the scripts in this repository, ensure you have R installed. You can install the required packages by running:

From R:
```R
install.packages(c("dplyr", "DBI", "RSQLite", "recount3"))
```

From Bash command-line:

```bash
git clone https://github.com/SoniaRuiz/recount3-database-project.git
cd recount3-database-project
```

## Usage

### [init.R](https://github.com/SoniaRuiz/recount3-database-project/blob/main/init.R)

The init.R script downloads publicly available junction data from any recount3 project and builds a junction database. 

```R
## Please, remember to update the variable `recount3_project_IDs`, to indicate the recount3 project ID that you would like to download and database.
source("init.R")
```

### [init_age.R](https://github.com/SoniaRuiz/recount3-database-project/blob/main/init_age.R)

The init_age.R script utilizes previously downloaded junction data from the GTEx project, stratifying samples by age before constructing a junction database. 

It starts the age sample clustering using the age groups "20-39", "40-59" and "60-79" years-old. Then, using the previously downloaded exon-exon junction data used for the creation of the Splicing database (init.R), this script clusters exon-exon split reads and count matrices across the samples of each age category. It pairs the split reads from the annotated category with the split reads from the novel donor and novel acceptor junctions across the samples of each age cluster at the tissue level. Finally, it creates the "Age Stratification" intron database.

```R
## To run the age stratification of the GTEx samples, it is necessary to have downloaded, processed and databased the GTEx v8 junctions using the `init.R` script.
source("init_age.R")
```

### [init_ENCODE.R](https://github.com/SoniaRuiz/recount3-database-project/blob/main/init_ENCODE.R)

The init_ENCODE.R script downloads BAM files from the ENCODE platform, extracts junctions, and generates a junction database.
BAM files downloaded correspond to the RNA-binding proteins involved in post-transcriptional processes published by [Van Nostrand et at. 2020](https://www.nature.com/articles/s41586-020-2077-3).

```R
source("init_ENCODE.R")
```

## Function Descriptions

### init.R:

**Data download and processing:**
* [*download_recount3_data()*](https://github.com/SoniaRuiz/recount3-database-project/blob/main/scripts/01_download_recount3_data.R): downloads, process and annotates, using a given Ensembl annotation, the exon-exon junction split read data from the recount3 project indicated.
* [*prepare_recount3_data()*](https://github.com/SoniaRuiz/recount3-database-project/blob/main/scripts/02_prepare_recount3_data.R): groups the processed split-read data by sample cluster for a given recount3 project.

**Junction pairing:**
* [*junction_pairing()*](https://github.com/SoniaRuiz/recount3-database-project/blob/main/scripts/04_junction_pairing.R): pairs novel junctions and annotated introns across the samples of each sample cluster.

**Junction processing prior databasing:**
* [*get_all_annotated_split_reads()*](https://github.com/SoniaRuiz/recount3-database-project/blob/main/scripts/08_get_all_annotated_split_reads.R): loops through the samples clusters from the current recount3 project and obtains all unique split reads found across their samples.
* [*get_all_raw_jxn_pairings()*](https://github.com/SoniaRuiz/recount3-database-project/blob/main/scripts/09_get_all_raw_jxn_pairings.R): loops through the samples clusters from the current recount3 project and obtains all junction pairings.
* [*tidy_data_pior_sql()*](https://github.com/SoniaRuiz/recount3-database-project/blob/main/scripts/12_tidy_data_pior_sql.R): discards ambiguous junctions and prepares the data prior generation of the SQL database.
* [*generate_transcript_biotype_percentage()*](https://github.com/SoniaRuiz/recount3-database-project/blob/main/scripts/10_generate_transcript_biotype_percentage.R): calculates the percentage of protein-coding transcripts in which a given junction may appear.
* [*generate_recount3_tpm()*](https://github.com/SoniaRuiz/recount3-database-project/blob/main/scripts/11_generate_recount3_tpm.R): obtains and transforms (scaled by library size) raw counts data. Calculates the median gene TPM across all samples from each sample cluster.

**SQLITE database generation:**
* [*sql_database_generation()*](https://github.com/SoniaRuiz/recount3-database-project/blob/main/scripts/14_sql_database_generation.R): creates the different Database tables sequentially (see [IntroVerse: a comprehensive database of introns across human tissues](https://doi.org/10.1093/nar/gkac1056))

### init_age.R:
* [*age_stratification_init_data()*](https://github.com/SoniaRuiz/recount3-database-project/blob/main/scripts/23_age_stratification_init_data.R): clusters the GTEx v8 samples by age supergroup, i.e. "20-39", "40-59" and "60-79" years-old.
* [*age_stratification_annotate()*](https://github.com/SoniaRuiz/recount3-database-project/blob/main/scripts/25_age_stratification_annotate.R): creates the split reads annotation files per age group.
* [*age_stratification_junction_pairing()*](https://github.com/SoniaRuiz/recount3-database-project/blob/main/scripts/26_age_stratification_junction_pairing.R): performs the split read junction pairing between novel junctions and annotated introns across the samples of each age group.
* get_all_annotated_split_reads()
* get_all_raw_jxn_pairings()
* generate_recount3_tpm()
* tidy_data_pior_sql()
* generate_transcript_biotype_percentage()
* sql_database_generation()

### init_ENCODE.R:
* [*ENCODE_download_metadata()*](https://github.com/SoniaRuiz/recount3-database-project/blob/main/scripts/28_ENCODE_download_metadata.R): downloads metadata from the gene-silencing knockdown experiments of RBPs from the ENCODE platform.
* [*ENCODE_download_bams()*](https://github.com/SoniaRuiz/recount3-database-project/blob/main/scripts/29_ENCODE_download_bams.R): downloads the BAM files corresponding to the gene-silencing knockdown experiments of RBPs from the ENCODE platform.
* [*prepare_encode_data()*](https://github.com/SoniaRuiz/recount3-database-project/blob/main/scripts/30_ENCODE_prepare_encode_data.R): extracts the exon-exon junction split read data from each knockdown ENCODE experiments, processes and annotates them.
* junction_pairing()
* get_all_annotated_split_reads()
* get_all_raw_jxn_pairings()
* tidy_data_pior_sql()
* generate_transcript_biotype_percentage()
* sql_database_generation()

## Environments
The code included within this repository has been successfully tested on:
* Ubuntu version "16.04.7 LTS (Xenial Xerus)"
* Ubuntu version "22.04.2 LTS (Jammy Jellyfish)"

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

## Acknowledgments

* All [RytenLab](https://rytenlab.com/) members
* [Aligning Science Across Parkinson's (ASAP)](https://parkinsonsroadmap.org/#)
* [UCL Great Ormond Street Institute Of Child Health](https://www.ucl.ac.uk/child-health/great-ormond-street-institute-child-health-0)

![Aligning Science Across Parkinson's (ASAP)](https://parkinsonsroadmap.org/wp-content/uploads/2020/10/cropped-ASAP_Logo_FullColor.png)
