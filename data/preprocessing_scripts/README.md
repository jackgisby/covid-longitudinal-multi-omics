# Data processing steps

We intend to deposit the raw RNA sequencing data in the European Phenome-Genome Archive. These reads were processed using the nf-core RNA-seq v3.2 pipeline, as described in our preprint: https://doi.org/10.1101/2022.04.29.22274267. Finally, htseq-count was used to generate a raw counts matrix. The script used to do this can be found in this repository `data/preprocessing_scripts/run_htseq`. This script was run on the Imperial High Performance Computing cluster. The resulting file is available in this manuscript: `data/htseq_counts.csv`. 

The proteomics data is available via Zenodo: https://doi.org/10.5281/zenodo.6497251. The original normalised data was received from SomaLogic as a .adat file. We converted these data from a .adat file to three .csv files (`sample_technical_meta.csv`, `feature_meta.csv` and `soma_abundance.csv`) using the script in `data/preprocessing_scripts/load_soma_from_adat.R`.

# Normalising the data

We provide RMD files that demonstrate the reading and analysis of the RNA sequencing and SomaLogic proteomics data in `notebooks/`. However, for a quick guide to loading these data, see below.

The functions for loading in the datasets are in `R/preprocessing.R`. The main functions have descriptions above the function definition. A minimal workflow would look something like this:

```
# Load the original raw counts
htseq_counts <- data.frame(fread("data/20210714_htseq_counts/Unnormalized_counts.csv"))

# Subset the counts to the wave 1 data by using this function and specifying the wave 1 metadata
wave1_counts <- get_summarized_experiment(htseq_counts, "scripts/20210713_salmon_counts_initial_look/w1_metadata.csv")

# Normalise the wave1 data in preparation for downstream processes (in this case we get logCPM)
wave1_normalised_counts <- normalize_se(wave1_counts, normalisation_level = "logcpm")

# Extract the original raw counts from the object (note: these are NOT logcpm transformed, they are the original HTSeq counts)
wave1_raw_htseq_counts <- assay(wave1_normalised_counts)
print(wave1_raw_htseq_counts[1:5,1:5])

# Extract the normalised counts from the object (note: you MUST specify 2 as the second argument to assay to get these)
wave1_normalised_count_matrix <- assay(wave1_normalised_counts, 2)
print(wave1_normalised_count_matrix[1:5,1:5])

# See the information for the genes
head(rowData(wave1_normalised_counts))

# See the sample metadata
head(colData(wave1_normalised_counts))
```

The main things to notice in the above script are:
 1. You can get wave1 or wave2 counts based on which metadata files you pass to the `get_summarized_experiment` function.
 2. You can choose the type of normalisation you want to perform on the counts by specifying the `normalisation_level` argument to the `normalize_se` function.
 3. The `normalisation_level` `SummarizedExperiment` object contains both raw and normalised counts. Make sure to pick the right assay when doing your analysis.

Further details and options can be found in the `R/preprocessing.R` file.

The process for loading the proteomics data is even simpler:
```

# When the argument ret_wide == TRUE, gets a list of three dataframes
# When the argument ret_wide == FALSE (default), returns a single dataframe in "long form"
wave1_soma <- get_soma_data(
  soma_abundance = "../data/soma_abundance.csv",
  sample_meta = "../data/w1_metadata.csv",
  sample_technical_meta = "../data/sample_technical_meta.csv",
  feature_meta = "../data/feature_meta.csv",
  ret_wide = TRUE
)

# View a sample of the normalised proteomics data
print(wave1_soma$soma_abundance[1:5,1:5])

# View the sample metadata
head(wave1_soma$sample_meta)

# View the feature metadata
head(wave1_soma$feature_meta)

```

The `get_soma_data` function performs a rank-based inverse normal transformation to the proteomics data, as described in the preprint: https://doi.org/10.1101/2022.04.29.22274267 
