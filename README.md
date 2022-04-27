# Longitudinal COVID-19 Multi-omics Repository

Repository for the paper: "Multi-omics before, during and after COVID-19 reveals severity-dependent trajectories and persistent pro-thrombotic transcriptome after recovery" by Gisby and Norzawani _et al._ doi: XX

In this repository, we have stored:
 - A processed version of the dataset used in the paper (`data/`)
 - Scripts that demonstrate the pre-processing applied to generate the datasets (`data/preprocessing_scripts/`)
 - RMD scripts replicating our major analyses (`notebooks/`)
 - Helper functions to process the data and perform basic analyses (`R/`)

## Description of the Dataset

RNA sequencing, SomaLogic proteomics and flow cytometry data were generated for two cohorts of end-stage kidney disease patients with COVID-19. For a full description of the cohort and methods used to generate these data, see the corresponding publication: XX

We have made the raw RNA sequencing reads available on the European Phenome-Genome archive: XX. In this repository, we have also included a counts matrix (`data/htseq_counts.csv`) generated for these reads. The steps used to generate this matrix are described in the corresponding publication (XX) and in the file `data/preprocessing_scripts/README.md`. The markdown documents in this repository (`notebooks/`) demonstrate the use of the `get_summarized_experiment` function to read and analyse these data.

We have made the SomaLogic proteomics and flow cytometry data available by means of Mendeley Data: XX. The proteomics data is also included in this repository by means of three files: `sample_technical_meta.csv`, `feature_meta.csv` and `soma_abundance.csv`. The first two files contain metadata columns for the samples and protein features, respectively. The final file includes the protein abundance data. The markdown documents in this repository (`notebooks/`) demonstrate the use of the `get_soma_data` function to read and analyse these data.

Finally, clinical metadata is available for the two cohorts described in this study. The Wave 1 cohort consists of samples collected from patients during the first wave of COVID-19 in early 2020 while samples were collected from the Wave 2 cohort the following year. A full description of this cohort is available in the corresponding publication: XX. This data is available through Mendeley Data (XX) and in this repository (`data/w1_metadata.csv`) and (`w2_metadata.csv`).

The features in the clinical metadata include:
Column Name | Data Type | Description
| :---: | :---: | :---:
sample_id | Character | Unique identifier for samples
individual_id | Character | Unique identifier for individuals
ethnicity | Character | The individual's ethnicity (asian, white, black or other)
sex | Character | The individual's sex (M or F)
calc_age | Integer | Age in years
ihd	| Character | Information on coronary heart disease
previous_vte | Character | Whether individuals have had venous thromboembolism
copd | Character | Whether individuals have chronic obstructive pulmonary disease
diabetes | Character | Whether individuals have diabetes, and, if so, the type of diabetes
smoking | Character | Smoking status
cause_eskd | Character | Cause of ESKD
WHO_severity | Character | The peak (WHO) severity for the patient over the disease course
WHO_temp_severity | Character | The (WHO) severity at time of sampling
fatal_disease | Logical | Whether the disease was fatal
case_control | Character | Whether the individual was COVID-19 `POSITIVE` or `NEGATIVE` at time of sampling. Convalescent patients are denoted by the label `RECOVERY`
radiology_evidence_covid | Character | Evidence of COVID-19 from radiology
time_from_first_symptoms | Integer | The number of days since the individual first experienced COVID symptoms at time of sampling
time_from_first_positive_swab | Integer | The number of days since the individual's first positive swab was taken at time of sampling

## Dependencies

Note that if you get an error (e.g. that a function cannot be found) it is likely that you need to load a particular package. Among the packages you might need to load are:
 - ggplot2
 - SummarizedExperiment
 - data.table
 - variancePartition
 - edgeR
 - BiocParallel
 - biomaRt
 - topGO
 - GO.db
 - splines
 - ReactomePA
 - org.Hs.eg.db

# License

This work is licensed under a XX
