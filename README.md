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

### init.R

The init.R script downloads publicly available junction data from any recount3 project and builds a junction database. 

```R
## Please, remember to update the variable `recount3_project_IDs`, to indicate the recount3 project ID that you would like to download and database.
source("init.R")
```

### init_age.R

The init_age.R script utilizes previously downloaded junction data from the GTEx project, stratifying samples by age before constructing a junction database. 

It starts the age sample clustering using the age groups "20-39", "40-59" and "60-79" years-old. Then, using the previously downloaded exon-exon junction data used for the creation of the Splicing database (init.R), this script clusters exon-exon split reads and count matrices across the samples of each age category. It pairs the split reads from the annotated category with the split reads from the novel donor and novel acceptor junctions across the samples of each age cluster at the tissue level. Finally, it creates the "Age Stratification" intron database.

```R
## To run the age stratification of the GTEx samples, it is necessary to have downloaded, processed and databased the GTEx v8 junctions using the `init.R` script.
source("init_age.R")
```

### init_ENCODE.R

The init_ENCODE.R script downloads BAM files from the ENCODE platform, extracts junctions, and generates a junction database.

```R
source("init_ENCODE.R")
```

## Function Descriptions

### init.R:

* download_recount3_data()
* prepare_recount3_data()
* junction_pairing()
* get_all_annotated_split_reads()
* get_all_raw_jxn_pairings()
* tidy_data_pior_sql()
* generate_transcript_biotype_percentage()
* generate_recount3_tpm()
* sql_database_generation()

### init_age.R:
* age_stratification_annotate()
* age_stratification_junction_pairing()
* get_all_annotated_split_reads()
* get_all_raw_jxn_pairings()
* generate_recount3_tpm()
* tidy_data_pior_sql()
* generate_transcript_biotype_percentage()
* sql_database_generation()

### init_ENCODE.R:
* ENCODE_download_metadata()
* ENCODE_download_bams()
* prepare_encode_data()
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
