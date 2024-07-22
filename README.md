## Splicing accuracy varies across human introns, tissues, age and disease

*Sonia Garcia-Ruiz, [David Zhang](https://github.com/dzhang32), [Emil K Gustavsson](https://github.com/egustavsson), Guillermo Rocamora-Perez, [Melissa Grant-Peters](https://github.com/mgrantpeters), Aine Fairbrother-Browne, Regina H Reynolds, Jonathan W Brenton, Ana L Gil-Martinez, Zhongbo Chen, Donald C Rio, Juan A Botia, Sebastian Guelfi, [Leonardo Collado-Torres](https://lcolladotor.github.io/), Mina Ryten*

bioRxiv 2023.03.29.534370;
doi: [https://doi.org/10.1101/2023.03.29.534370](https://doi.org/10.1101/2023.03.29.534370)


## Repository 

This 'Breadcrumbsrecount3-database-project' repository contains the code used to generate all the databased produced for the manuscript [**Splicing accuracy varies across human introns, tissues and age**](https://doi.org/10.1101/2023.03.29.534370).


### 1. *"Splicing"* Database Generation

To produce the *"Splicing"* intron database, please follow the file pipeline indicated below:

*init.R*. Main file. It starts the exon-exon junction data download from recount3. It annotates the split reads in three categories: annotated intron, novel donor and novel acceptor junction. Once the split reads are annotated, this script pairs the split reads from the annotated category with those from the novel donor and novel acceptor junction groups across the samples of each tissue. Finally, after accessing multiple auxiliary functions, it creates the "Splicing" intron database using SQL commands ('DBI' R package).

Important auxiliary functions called from the *init.R* file:

1. download_recount3_data()
2. prepare_recount3_data()
3. junction_pairing()
4. get_all_annotated_split_reads()
5. get_all_raw_jxn_pairings()
6. tidy_data_pior_sql()
7. generate_transcript_biotype_percentage()
8. generate_recount3_tpm()
9. sql_database_generation()
  

### 2. *"Age Stratification"* Database Generation

To produce the *"Age Stratification"* intron database, please follow the file pipeline indicated below:

*init_age.R*. Main file. It starts the GTEx v8 sample clustering by age within the age groups "20-39", "40-59" and "60-79" years-old. Then, using the previously downloaded and QC'ed exon-exon junction data used for the creation of the Splicing database, this script groups the exon-exon split reads and extracts read count metrics across the samples of each age category. Then, it pairs the split reads from the annotated category with the split reads from the novel donor and novel acceptor junctions across the samples of each age cluster at the tissue level. Finally, after some intermediary calls to multiple auxiliary functions, it starts the SQL commands and creates the "Age Stratification" intron database.

Important auxiliary functions called from the *init_age.R* file:

1. age_stratification_init_data()
2. age_stratification_annotate()
3. age_stratification_junction_pairing()
4. get_all_annotated_split_reads()
5. get_all_raw_jxn_pairings()
6. tidy_data_pior_sql()
7. generate_transcript_biotype_percentage()
8. generate_recount3_tpm()
9. sql_database_generation()

## Environments

The code included within this repository has been successfully tested on:
* Ubuntu version "16.04.7 LTS (Xenial Xerus)"
* Ubuntu version "22.04.2 LTS (Jammy Jellyfish)"

## Acknowledgments

* All [RytenLab](https://rytenlab.com/) members
* [Aligning Science Across Parkinson's (ASAP)](https://parkinsonsroadmap.org/#)
* [UCL Great Ormond Street Institute Of Child Health](https://www.ucl.ac.uk/child-health/great-ormond-street-institute-child-health-0)

![Aligning Science Across Parkinson's (ASAP)](https://parkinsonsroadmap.org/wp-content/uploads/2020/10/cropped-ASAP_Logo_FullColor.png)

