Differential Expression Analysis
================
Jack Gisby
2022-04-28

-   [Load data and code](#load-data-and-code)
    -   [Load packages and custom
        functions](#load-packages-and-custom-functions)
    -   [Load RNA-seq data](#load-rna-seq-data)
    -   [Load SOMA protein data](#load-soma-protein-data)
-   [Differential expression analyses of COVID-19
    status](#differential-expression-analyses-of-covid-19-status)
    -   [RNA-seq analysis using dream](#rna-seq-analysis-using-dream)
    -   [Proteomic analysis using
        lmerTest](#proteomic-analysis-using-lmertest)
-   [Differential expression analyses of COVID-19 severity at time of
    sampling](#differential-expression-analyses-of-covid-19-severity-at-time-of-sampling)
    -   [RNA-seq analysis using dream](#rna-seq-analysis-using-dream-1)
    -   [Proteomics](#proteomics)
-   [Session info](#session-info)

# Load data and code

## Load packages and custom functions

This chunk loads packages required for the notebook in addition to
custom functions in the `../R` directory. These include functions for
loading data, running differential expression and enrichment analyses
and making plots. The main functions and their arguments are described
by comments within these files.

``` r
# load relevant packages
library(data.table)
library(SummarizedExperiment)
library(variancePartition)
library(BiocParallel)
library(edgeR)
library(lmerTest)
library(GSVA)
library(biomaRt)

# set up plotting settings
library(ggplot2)
library(ggpubr)

theme_set(theme_pubr(border = TRUE))
theme_update(text = element_text(size = 8))

# load custom functions from the R/ directory
source("../R/preprocessing.R")
source("../R/de.R")

set.seed(1)
```

## Load RNA-seq data

The raw RNA-seq counts are loaded, transformed into
`SummarizedExperiment` objects, and split into the “Wave 1” and “Wave 2”
cohorts. `SnowParam` is also initialised for the differential expression
analysis using the `variancePartition` package.

``` r
# get the raw RNA-seq counts
htseq_counts <- data.frame(fread("../data/htseq_counts.csv"))

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl", host="uswest.ensembl.org"))

# get the wave 1 and wave 2 counts and their associated metadata
wave1_counts <- get_summarized_experiment(
  htseq_counts,
  "../data/w1_metadata.csv",
  mart = mart
)

wave2_counts <- get_summarized_experiment(
  htseq_counts,
  "../data/w2_metadata.csv",
  mart = mart
)

# prepare multiprocessing for RNA-seq differential expression
num_cores <- 5
param <- SnowParam(num_cores, "SOCK", progressbar = TRUE)
register(param)
```

## Load SOMA protein data

The SOMA data is loaded into long format and split into the “Wave 1” and
“Wave 2” cohorts by the `get_soma_data` function.

``` r
# get the data for wave 1
wave1_soma <- get_soma_data(
  soma_abundance = "../data/soma_abundance.csv",
  sample_meta = "../data/w1_metadata.csv",
  sample_technical_meta = "../data/sample_technical_meta.csv",
  feature_meta = "../data/feature_meta.csv",
  ret_wide = FALSE
)

# get the data for wave 2
wave2_soma <- get_soma_data(
  soma_abundance = "../data/soma_abundance.csv",
  sample_meta = "../data/w2_metadata.csv",
  sample_technical_meta = "../data/sample_technical_meta.csv",
  feature_meta = "../data/feature_meta.csv",
  ret_wide = FALSE
)
```

# Differential expression analyses of COVID-19 status

The following code chunks apply linear mixed models to genes and
proteins for the analysis of COVID-19 status. The formulae applied is
described in the Methods section of the associated publication (XX)
under the heading “Differential expression analysis: COVID-19 positive
versus negative”.

## RNA-seq analysis using dream

In order to analyse the effect of COVID-19 status on the expression of
individual genes, we have used the `dream` function from the
`variancePartition` package.

The function `normalize_se` applies normalisation using edgeR. The
`normalisation_level` argument is set to “pre_dream” in order to
retrieve the data in the correct format for dream. The `group` argument
is set to COVID-19 status, which is passed to the `group` argument of
`edgeR::filterByExpr`.

The `run_dream` function is used to apply dream to compare the
expression of genes based on COVID-19 status.

The analysis is performed for each of the Wave 1 and Wave 2 cohorts.

### Wave 1 transcriptomic COVID-19 infected vs. non-infected comparison

``` r
# apply normalisation to counts (returns edgeR DGE object)
wave1_case_control_dge <- normalize_se(
  wave1_counts,
  normalisation_level = "pre_dream",
  group = "case_control"
)

# run dream differential expression analysis (from variancePartition package)
wave1_case_control_de <- run_dream(
  f = ~ case_control + sex + ethnicity + calc_age + (1 | individual_id),
  dge = wave1_case_control_dge,
  L = list(case_control = "case_controlPOSITIVE"),
  L_is_contrasts = FALSE
)
```

``` r
write.csv(
    wave1_case_control_de$case_control_tt,
    "../results/1_differential_expression/wave_1/gene_negative_vs_positive.csv",
    row.names = FALSE
)

head(wave1_case_control_de$case_control_tt)
```

    ##           gencode_id      ensembl_id gene_id    logFC  AveExpr         t
    ## 1 ENSG00000148773.14 ENSG00000148773   MKI67 1.934942 5.014080 10.400552
    ## 2 ENSG00000131747.15 ENSG00000131747   TOP2A 1.773565 3.883617 10.229445
    ## 3  ENSG00000276043.5 ENSG00000276043   UHRF1 1.629814 2.634389 10.030059
    ## 4 ENSG00000175063.17 ENSG00000175063   UBE2C 2.054433 1.138486  9.735528
    ## 5 ENSG00000145386.11 ENSG00000145386   CCNA2 1.974922 2.360058  9.689831
    ## 6 ENSG00000126787.13 ENSG00000126787  DLGAP5 2.314382 1.510518  9.641550
    ##        P.Value    adj.P.Val    z.std
    ## 1 5.082284e-19 7.828368e-15 8.910457
    ## 2 8.063416e-19 7.828368e-15 8.859140
    ## 3 2.512331e-18 1.626064e-14 8.731538
    ## 4 9.193332e-18 4.253827e-14 8.583620
    ## 5 1.095387e-17 4.253827e-14 8.563451
    ## 6 1.501622e-17 4.859499e-14 8.527023

### Wave 2 transcriptomic COVID-19 infected vs. pre-infection comparison

``` r
# apply normalisation to counts (returns edgeR DGE object)
# only NEGATIVE (pre-infection) and POSITIVE samples are compared in this analysis
wave2_case_control_dge <- normalize_se(
  wave2_counts[, wave2_counts$case_control != "RECOVERY"],
  normalisation_level = "pre_dream",
  group = "case_control"
)

# run dream differential expression analysis (from variancePartition package)
wave2_case_control_de <- run_dream(
  f = ~ case_control + sex + ethnicity + calc_age + (1 | individual_id),
  dge = wave2_case_control_dge,
  L = list(case_control = "case_controlPOSITIVE"),
  L_is_contrasts = FALSE
)
```

``` r
write.csv(
    wave2_case_control_de$case_control_tt, 
    "../results/1_differential_expression/wave_2/gene_negative_vs_positive.csv", 
    row.names = FALSE
)

print(head(wave2_case_control_de$case_control_tt))
```

    ##           gencode_id      ensembl_id         gene_id      logFC   AveExpr
    ## 1 ENSG00000165949.12 ENSG00000165949           IFI27  5.3186178 7.1994212
    ## 2  ENSG00000227081.5 ENSG00000227081 ENSG00000227081 -0.8327455 1.1359637
    ## 3  ENSG00000169385.3 ENSG00000169385          RNASE2  2.2437637 6.0478778
    ## 4  ENSG00000243199.1 ENSG00000243199 ENSG00000243199 -0.9640504 0.7854245
    ## 5  ENSG00000234664.1 ENSG00000234664         HMGN2P5  0.9511312 2.5923348
    ## 6 ENSG00000088827.13 ENSG00000088827         SIGLEC1  3.3545174 8.3205458
    ##            t      P.Value    adj.P.Val     z.std
    ## 1  11.024030 2.365192e-18 4.571206e-14  8.738360
    ## 2 -10.455787 3.763076e-17 3.636449e-13 -8.420055
    ## 3   9.891359 5.464663e-16 3.520518e-12  8.100694
    ## 4  -9.479492 4.022283e-15 1.943466e-11 -7.854232
    ## 5   8.729579 1.437078e-13 4.806420e-10  7.392856
    ## 6   8.703161 1.492136e-13 4.806420e-10  7.387857

### Wave 2 transcriptomic COVID-19 recovery vs. pre-infection comparison

In addition to comparing COVID-19 positive and negative individuals, we
were able to compare pre-infection samples to samples taken after
recovery from COVID-19 for twelve of the patients of the Wave 2 cohort.

``` r
# only NEGATIVE (pre-infection) and RECOVERY samples are compared in this analysis
wave2_recov_counts <- wave2_counts[, wave2_counts$case_control != "POSITIVE"]

# some patients do not have a RECOVERY sample
# since this is a paired analysis, these individual's pre-infection timepoints are removed
individuals_with_recovery_sample <- wave2_recov_counts$individual_id[wave2_recov_counts$case_control == "RECOVERY"]
wave2_recov_counts <- wave2_recov_counts[,  wave2_recov_counts$individual_id %in% individuals_with_recovery_sample]

# apply normalisation to counts (returns edgeR DGE object)
wave2_recovery_dge <- normalize_se(
  wave2_recov_counts, 
  normalisation_level = "pre_dream", 
  group = "case_control"
)

# run dream differential expression analysis (from variancePartition package)
wave2_recovery_de <- run_dream(
  f = ~ case_control + sex + ethnicity + calc_age + (1 | individual_id),
  dge = wave2_recovery_dge,
  L = list(case_control = "case_controlRECOVERY"),
  L_is_contrasts = FALSE
)
```

``` r
write.csv(
    wave2_recovery_de$case_control_tt, 
    "../results/1_differential_expression/wave_2/gene_negative_vs_recovery.csv", 
    row.names = FALSE
)

print(head(wave2_recovery_de$case_control_tt))
```

    ##           gencode_id      ensembl_id         gene_id     logFC   AveExpr
    ## 1  ENSG00000236304.1 ENSG00000236304 ENSG00000236304 1.8376567 2.7579686
    ## 2 ENSG00000163430.12 ENSG00000163430           FSTL1 1.1722049 1.0489600
    ## 3 ENSG00000095303.17 ENSG00000095303           PTGS1 0.6416105 6.1308448
    ## 4 ENSG00000134668.12 ENSG00000134668          SPOCD1 1.3122343 0.7431222
    ## 5  ENSG00000163735.7 ENSG00000163735           CXCL5 1.8402389 1.5466860
    ## 6 ENSG00000108839.12 ENSG00000108839          ALOX12 1.3330727 2.5465500
    ##           t      P.Value    adj.P.Val    z.std
    ## 1 12.930396 5.372831e-08 0.0009198286 5.438509
    ## 2 11.627480 1.597729e-07 0.0013676560 5.240926
    ## 3 10.718000 3.694085e-07 0.0021080913 5.084082
    ## 4 10.200070 6.411830e-07 0.0025292432 4.978381
    ## 5  9.871564 7.386808e-07 0.0025292432 4.950908
    ## 6  9.612284 1.237765e-06 0.0035317571 4.849494

## Proteomic analysis using lmerTest

We analysed the effect of COVID-19 status on individuals proteins using
linear mixed models applied using the `lmerTest` package. This is a
similar analytical strategy to that applied to individual genes, however
the `variancePartition::dream` function applies parameter sharing.

The analysis is performed for each of the Wave 1 and Wave 2 cohorts.

### Wave 1 proteomic COVID-19 infected vs. non-infected comparison

``` r
# fit linear mixed models for proteins
wave1_case_control_models <- soma_de(
  soma_long = wave1_soma,
  formula_string = "RFU ~ case_control + sex + ethnicity + calc_age + (1 | individual_id)"
)

# takes a list of models and applies anova to each one
# the anova results are concatenated, written as a CSV, and returned as a data.frame
models_to_df <- function(m, file_path) {
    
    model_df <- data.frame()
    
    for (i in 1:length(m)) {
        
        anova_df <- data.frame(anova(m[[i]]))
        colnames(anova_df) <- c("sum_sq", "mean_sq", "num_df", "den_df", "f_value", "p_value")
        
        anova_df$feature <- names(m)[i]
        anova_df$term <- rownames(anova_df)

        model_df <- rbind(model_df, anova_df)
    }
    
    rownames(model_df) <- NULL
    model_df$adjusted_p_value <- NA
    
    for (term in unique(model_df$term)) {
        
        model_df$adjusted_p_value[model_df$term == term] <- p.adjust(model_df$p_value[model_df$term == term], method = "BH")
    }
    
    write.csv(model_df, file_path, row.names = FALSE)
    
    return(model_df)
}

wave1_case_control_df <- models_to_df(
    wave1_case_control_models, 
    "../results/1_differential_expression/wave_1/protein_negative_vs_positive.csv"
)

head(wave1_case_control_df[wave1_case_control_df$term == "case_control" ,])
```

    ##         sum_sq     mean_sq num_df   den_df   f_value    p_value      feature
    ## 1  1.684838230 1.684838230      1 79.46688 6.2934101 0.01415361 seq.10000.28
    ## 5  0.180439368 0.180439368      1 80.64892 0.6948487 0.40698207  seq.10001.7
    ## 9  1.381962637 1.381962637      1 77.54401 6.6497796 0.01180947 seq.10003.15
    ## 13 2.357238638 2.357238638      1 83.51676 8.6209086 0.00429230 seq.10006.25
    ## 17 0.003197357 0.003197357      1 81.88581 0.0153449 0.90171819 seq.10008.43
    ## 21 0.074366091 0.074366091      1 76.11538 0.3532911 0.55401867 seq.10010.10
    ##            term adjusted_p_value
    ## 1  case_control       0.04895944
    ## 5  case_control       0.54679107
    ## 9  case_control       0.04263130
    ## 13 case_control       0.01926612
    ## 17 case_control       0.93517257
    ## 21 case_control       0.67413364

### Wave 2 proteomic COVID-19 infected vs. non-infected comparison

``` r
# fit linear mixed models for proteins
wave2_case_control_models <- soma_de(
  soma_long = wave2_soma,
  formula_string = "RFU ~ case_control + sex + ethnicity + calc_age + (1 | individual_id)"
)

wave2_case_control_df <- models_to_df(
    wave2_case_control_models, 
    "../results/1_differential_expression/wave_2/protein_negative_vs_positive.csv"
)

head(wave2_case_control_df[wave2_case_control_df$term == "case_control" ,])
```

    ##       sum_sq   mean_sq num_df   den_df  f_value      p_value      feature
    ## 1   5.904916  5.904916      1 100.1205 11.17815 1.165218e-03 seq.10000.28
    ## 5  41.378870 41.378870      1 100.2807 86.82358 3.025937e-15  seq.10001.7
    ## 9  22.035503 22.035503      1 100.1466 57.86034 1.567094e-11 seq.10003.15
    ## 13 10.426370 10.426370      1 100.3637 32.45096 1.225230e-07 seq.10006.25
    ## 17 14.312679 14.312679      1 100.2665 29.73308 3.572577e-07 seq.10008.43
    ## 21 27.224366 27.224366      1 100.4696 75.68207 6.729278e-14 seq.10010.10
    ##            term adjusted_p_value
    ## 1  case_control     1.524427e-03
    ## 5  case_control     1.318158e-14
    ## 9  case_control     3.616755e-11
    ## 13 case_control     2.055298e-07
    ## 17 case_control     5.813769e-07
    ## 21 case_control     2.154798e-13

# Differential expression analyses of COVID-19 severity at time of sampling

In addition to exploring the effect of COVID-19 status on gene
expression and protein abundance, we modelled COVID-19 severity at time
of sampling. Severity (`WHO_temp_severity`) is encoded as an ordered
factor. By default, polynomial contrasts (`?contr.poly`) are generated
for ordered factors, meaning that the linear mixed models consider
severity to be an ordinal variable.

## RNA-seq analysis using dream

Like for the analysis of COVID-19 status, we explored the expression of
individual genes using `variancePartition::dream`. The analysis was run
for each cohort separately. Only COVID-19 positive samples were modelled
in this analysis.

### Wave 1 transcriptomic COVID-19 severity analysis

``` r
# convert contemporaneous severity to an ordered factor
wave1_counts$WHO_temp_severity <- ordered(wave1_counts$WHO_temp_severity, c("mild", "moderate", "severe", "critical"))

# get the normalised dataset for COVID-19 cases
wave1_severity_dge <- normalize_se(
  wave1_counts[, wave1_counts$case_control == "POSITIVE" & !is.na(wave1_counts$WHO_temp_severity)],
  normalisation_level = "pre_dream",
  group = "WHO_temp_severity"
)

# run dream differential expression analysis (from variancePartition package)
wave1_severity_de <- run_dream(
  f = ~ WHO_temp_severity + sex + ethnicity + calc_age + (1 | individual_id),
  dge = wave1_severity_dge,
  L = list("WHO_temp_severity" = c("WHO_temp_severity.L", "WHO_temp_severity.Q", "WHO_temp_severity.C")),
  L_is_contrasts = FALSE
)
```

``` r
write.csv(
    wave1_severity_de$WHO_temp_severity_tt, 
    "../results/1_differential_expression/wave_1/gene_severity.csv", 
    row.names = FALSE
)

print(head(wave1_severity_de$WHO_temp_severity_tt))
```

    ##           gencode_id      ensembl_id  gene_id WHO_temp_severity.L
    ## 1 ENSG00000078699.22 ENSG00000078699  CBFA2T2          -0.3263431
    ## 2 ENSG00000084073.10 ENSG00000084073 ZMPSTE24           0.3727049
    ## 3  ENSG00000211891.6 ENSG00000211891     IGHE           2.1268354
    ## 4 ENSG00000112787.14 ENSG00000112787   FBRSL1          -0.3184568
    ## 5 ENSG00000198431.16 ENSG00000198431   TXNRD1           0.4431768
    ## 6  ENSG00000211637.2 ENSG00000211637 IGLV4-69           1.2719832
    ##   WHO_temp_severity.Q WHO_temp_severity.C    AveExpr        F      P.Value
    ## 1         0.014977459         0.001258513  5.7266858 24.27277 6.474893e-13
    ## 2        -0.102005051        -0.042670446  5.3542315 22.61247 3.494240e-12
    ## 3         0.449762357         0.312872877 -0.2412975 21.12926 1.624748e-11
    ## 4         0.008961688         0.015378854  7.1838742 20.55258 2.977162e-11
    ## 5        -0.067352645        -0.012781824  6.1270927 19.90168 5.930376e-11
    ## 6        -1.447146303         0.597863069  2.6456058 19.59481 8.223856e-11
    ##      adj.P.Val    F.std
    ## 1 1.355648e-08 19.93449
    ## 2 3.657945e-08 18.79161
    ## 3 1.133912e-07 17.74869
    ## 4 1.558321e-07 17.33741
    ## 5 2.483285e-07 16.86921
    ## 6 2.869715e-07 16.64699

### Wave 2 transcriptomic COVID-19 severity analysis

``` r
# convert contemporaneous severity to an ordered factor
wave2_counts$WHO_temp_severity <- ordered(wave2_counts$WHO_temp_severity, c("mild", "moderate", "severe", "critical"))

# get the normalised dataset for COVID-19 cases
wave2_severity_dge <- normalize_se(
  wave2_counts[, wave2_counts$case_control == "POSITIVE" & !is.na(wave2_counts$WHO_temp_severity)],
  normalisation_level = "pre_dream",
  group = "WHO_temp_severity"
)

# run dream differential expression analysis (from variancePartition package)
wave2_severity_de <- run_dream(
  f = ~ WHO_temp_severity + sex + ethnicity + calc_age + (1 | individual_id),
  dge = wave2_severity_dge,
  L = list(WHO_temp_severity = c("WHO_temp_severity.L", "WHO_temp_severity.Q", "WHO_temp_severity.C")),
  L_is_contrasts = FALSE
)
```

``` r
write.csv(
    wave2_severity_de$WHO_temp_severity_tt, 
    "../results/1_differential_expression/wave_2/gene_severity.csv", 
    row.names = FALSE
)

print(head(wave2_severity_de$WHO_temp_severity_tt))
```

    ##           gencode_id      ensembl_id gene_id WHO_temp_severity.L
    ## 1 ENSG00000154258.17 ENSG00000154258   ABCA9          -1.2497979
    ## 2 ENSG00000003436.16 ENSG00000003436    TFPI           1.4948912
    ## 3  ENSG00000119632.4 ENSG00000119632 IFI27L2           0.2767223
    ## 4  ENSG00000124635.9 ENSG00000124635  H2BC11           0.8399268
    ## 5  ENSG00000270276.2 ENSG00000270276   H4C15           0.8269051
    ## 6  ENSG00000141744.4 ENSG00000141744    PNMT           0.4813776
    ##   WHO_temp_severity.Q WHO_temp_severity.C   AveExpr        F      P.Value
    ## 1           0.7352337         -0.18088670 -0.983901 21.19820 3.860432e-10
    ## 2          -0.1214807         -0.44856405 -0.286022 19.48483 1.596296e-09
    ## 3          -0.3512347         -0.05790458  3.939837 17.23648 1.114017e-08
    ## 4          -0.6593023         -0.18485987  2.130155 17.14604 1.206997e-08
    ## 5          -1.3196204         -0.92314151 -4.380679 15.64057 4.695367e-08
    ## 6          -1.2354299         -0.82861937 -1.788681 15.55184 5.094144e-08
    ##      adj.P.Val    F.std
    ## 1 7.989550e-06 15.59520
    ## 2 1.651847e-05 14.62842
    ## 3 6.245001e-05 13.30293
    ## 4 6.245001e-05 13.24818
    ## 5 1.412503e-04 12.31950
    ## 6 1.412503e-04 12.26372

## Proteomics

We then analysed the effect of 4-level COVID-19 WHO severity on the
abundance of individual proteins using linear mixed models. As before,
the analysis was applied to each cohort separately.

### Wave 1 proteomic COVID-19 severity analysis

``` r
# convert contemporaneous severity to an ordered factor
wave1_soma$WHO_temp_severity <- ordered(wave1_soma$WHO_temp_severity, c("mild", "moderate", "severe", "critical"))

# get dataset for COVID-19 cases
wave1_soma_positive <- wave1_soma[wave1_soma$case_control == "POSITIVE" & !is.na(wave1_soma$WHO_temp_severity) ,]

# fit linear mixed models for proteins
wave1_severity_models <- soma_de(
  soma_long = wave1_soma_positive,
  formula_string = "RFU ~ WHO_temp_severity + sex + ethnicity + calc_age + (1 | individual_id)"
)

wave1_severity_df <- models_to_df(
    wave1_severity_models, 
    "../results/1_differential_expression/wave_1/protein_severity.csv"
)

head(wave1_severity_df[wave1_severity_df$term == "WHO_temp_severity" ,])
```

    ##        sum_sq    mean_sq num_df   den_df   f_value   p_value      feature
    ## 1  0.09979642 0.03326547      3 68.75061 0.1172321 0.9497044 seq.10000.28
    ## 5  0.36065428 0.12021809      3 70.98834 0.4428789 0.7230729  seq.10001.7
    ## 9  0.28694214 0.09564738      3 67.58893 0.4359281 0.7279933 seq.10003.15
    ## 13 0.99714151 0.33238050      3 59.99785 1.1634360 0.3312199 seq.10006.25
    ## 17 0.12707917 0.04235972      3 62.60548 0.1897272 0.9030149 seq.10008.43
    ## 21 0.47171063 0.15723688      3 66.27665 0.7205405 0.5432580 seq.10010.10
    ##                 term adjusted_p_value
    ## 1  WHO_temp_severity        0.9790994
    ## 5  WHO_temp_severity        0.8818120
    ## 9  WHO_temp_severity        0.8834515
    ## 13 WHO_temp_severity        0.6837844
    ## 17 WHO_temp_severity        0.9632869
    ## 21 WHO_temp_severity        0.8026711

### Wave 2 proteomic COVID-19 severity analysis

``` r
# convert contemporaneous severity to an ordered factor
wave2_soma$WHO_temp_severity <- ordered(wave2_soma$WHO_temp_severity, c("mild", "moderate", "severe", "critical"))

# get dataset for COVID-19 cases
wave2_soma_severity <- wave2_soma[wave2_soma$case_control == "POSITIVE" & !is.na(wave2_soma$WHO_temp_severity) ,]

# fit linear mixed models for proteins
wave2_severity_models <- soma_de(
  soma_long = wave2_soma_severity,
  formula_string = "RFU ~ WHO_temp_severity + sex + ethnicity + calc_age + (1 | individual_id)"
)

wave2_severity_df <- models_to_df(
    wave2_severity_models, 
    "../results/1_differential_expression/wave_2/protein_severity.csv"
)

head(wave2_severity_df[wave2_severity_df$term == "WHO_temp_severity" ,])
```

    ##       sum_sq   mean_sq num_df   den_df   f_value     p_value      feature
    ## 1  8.2081314 2.7360438      3 88.07090 5.1503855 0.002505862 seq.10000.28
    ## 5  4.8600287 1.6200096      3 73.14561 3.5747943 0.017941590  seq.10001.7
    ## 9  0.9786665 0.3262222      3 89.40482 0.8836756 0.452777505 seq.10003.15
    ## 13 1.8353434 0.6117811      3 80.72684 2.0032240 0.120076186 seq.10006.25
    ## 17 2.3253638 0.7751213      3 78.96428 1.6541174 0.183674658 seq.10008.43
    ## 21 2.3037738 0.7679246      3 70.65678 2.0387521 0.116230452 seq.10010.10
    ##                 term adjusted_p_value
    ## 1  WHO_temp_severity       0.01018666
    ## 5  WHO_temp_severity       0.04963588
    ## 9  WHO_temp_severity       0.56016903
    ## 13 WHO_temp_severity       0.21752394
    ## 17 WHO_temp_severity       0.29742888
    ## 21 WHO_temp_severity       0.21266580

# Session info

Session information collected at the point this markdown document was
knit.

``` r
sessionInfo()
```

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19044)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United Kingdom.1252 
    ## [2] LC_CTYPE=English_United Kingdom.1252   
    ## [3] LC_MONETARY=English_United Kingdom.1252
    ## [4] LC_NUMERIC=C                           
    ## [5] LC_TIME=English_United Kingdom.1252    
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] ggpubr_0.4.0                biomaRt_2.50.3             
    ##  [3] GSVA_1.42.0                 lmerTest_3.1-3             
    ##  [5] lme4_1.1-28                 Matrix_1.3-4               
    ##  [7] edgeR_3.36.0                variancePartition_1.24.0   
    ##  [9] BiocParallel_1.28.3         limma_3.50.1               
    ## [11] ggplot2_3.3.5               SummarizedExperiment_1.24.0
    ## [13] Biobase_2.54.0              GenomicRanges_1.46.1       
    ## [15] GenomeInfoDb_1.30.1         IRanges_2.28.0             
    ## [17] S4Vectors_0.32.3            BiocGenerics_0.40.0        
    ## [19] MatrixGenerics_1.6.0        matrixStats_0.61.0         
    ## [21] data.table_1.14.2          
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] minqa_1.2.4                 colorspace_2.0-2           
    ##   [3] ggsignif_0.6.3              ellipsis_0.3.2             
    ##   [5] XVector_0.34.0              rstudioapi_0.13            
    ##   [7] bit64_4.0.5                 AnnotationDbi_1.56.2       
    ##   [9] fansi_1.0.2                 xml2_1.3.3                 
    ##  [11] codetools_0.2-18            splines_4.1.2              
    ##  [13] sparseMatrixStats_1.6.0     doParallel_1.0.17          
    ##  [15] cachem_1.0.6                knitr_1.37                 
    ##  [17] nloptr_2.0.0                pbkrtest_0.5.1             
    ##  [19] broom_0.7.12                annotate_1.72.0            
    ##  [21] dbplyr_2.1.1                png_0.1-7                  
    ##  [23] graph_1.72.0                HDF5Array_1.22.1           
    ##  [25] compiler_4.1.2              httr_1.4.2                 
    ##  [27] backports_1.4.1             assertthat_0.2.1           
    ##  [29] fastmap_1.1.0               cli_3.1.1                  
    ##  [31] BiocSingular_1.10.0         htmltools_0.5.2            
    ##  [33] prettyunits_1.1.1           tools_4.1.2                
    ##  [35] rsvd_1.0.5                  gtable_0.3.0               
    ##  [37] glue_1.6.1                  GenomeInfoDbData_1.2.7     
    ##  [39] reshape2_1.4.4              dplyr_1.0.8                
    ##  [41] rappdirs_0.3.3              Rcpp_1.0.8                 
    ##  [43] carData_3.0-5               vctrs_0.3.8                
    ##  [45] Biostrings_2.62.0           rhdf5filters_1.6.0         
    ##  [47] nlme_3.1-153                DelayedMatrixStats_1.16.0  
    ##  [49] iterators_1.0.14            xfun_0.29                  
    ##  [51] stringr_1.4.0               beachmat_2.10.0            
    ##  [53] irlba_2.3.5                 lifecycle_1.0.1            
    ##  [55] gtools_3.9.2                rstatix_0.7.0              
    ##  [57] XML_3.99-0.8                zlibbioc_1.40.0            
    ##  [59] MASS_7.3-54                 scales_1.1.1               
    ##  [61] hms_1.1.1                   parallel_4.1.2             
    ##  [63] rhdf5_2.38.0                curl_4.3.2                 
    ##  [65] SingleCellExperiment_1.16.0 yaml_2.2.2                 
    ##  [67] memoise_2.0.1               stringi_1.7.6              
    ##  [69] RSQLite_2.2.10              ScaledMatrix_1.2.0         
    ##  [71] foreach_1.5.2               filelock_1.0.2             
    ##  [73] caTools_1.18.2              boot_1.3-28                
    ##  [75] rlang_1.0.1                 pkgconfig_2.0.3            
    ##  [77] bitops_1.0-7                evaluate_0.15              
    ##  [79] lattice_0.20-45             Rhdf5lib_1.16.0            
    ##  [81] purrr_0.3.4                 bit_4.0.4                  
    ##  [83] tidyselect_1.1.1            GSEABase_1.56.0            
    ##  [85] plyr_1.8.6                  magrittr_2.0.2             
    ##  [87] R6_2.5.1                    snow_0.4-4                 
    ##  [89] gplots_3.1.1                generics_0.1.2             
    ##  [91] DelayedArray_0.20.0         DBI_1.1.2                  
    ##  [93] pillar_1.7.0                withr_2.4.3                
    ##  [95] abind_1.4-5                 KEGGREST_1.34.0            
    ##  [97] RCurl_1.98-1.6              tibble_3.1.6               
    ##  [99] car_3.0-12                  crayon_1.5.0               
    ## [101] KernSmooth_2.23-20          utf8_1.2.2                 
    ## [103] BiocFileCache_2.2.1         rmarkdown_2.11             
    ## [105] progress_1.2.2              locfit_1.5-9.4             
    ## [107] grid_4.1.2                  blob_1.2.2                 
    ## [109] digest_0.6.29               xtable_1.8-4               
    ## [111] tidyr_1.2.0                 numDeriv_2016.8-1.1        
    ## [113] munsell_0.5.0
