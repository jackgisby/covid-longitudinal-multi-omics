# covid-longitudinal-rnaseq

Functions for longitudinal RNA-seq project

## Normalising the counts data

The two main functions to worry about are in `R/preprocessing.R`. Each of these have descriptions above the function definition. A minimal workflow would look something like this:

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
