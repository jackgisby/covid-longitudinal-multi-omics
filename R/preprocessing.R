#' Get the wave 1 or wave 2 RNA-seq data as a `SummarizedExperiment`
#' 
#' This should be used to generate a `SummarizedExperiment` object to be
#' passed to `normalize_se` for normalisation.
#' 
#' @param counts
#' The HTSeq counts matrix. This could be loaded as in the following command:
#' `data.frame(fread("data/20210714_htseq_counts/Unnormalized_counts.csv"))`
#' 
#' @param metadata_path
#' The path to either the wave 1 or wave 2 metadata. Depending on which metadata
#' you pass to this function, you will get the corresponding wave 1 or wave 2 
#' samples in the returned `SummarizedExperiment.`
#' 
#' @param mart
#' The line `mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))`
#' is used to get gene names from ensembl, however this requires an internet
#' connection, is somewhat time-consuming and puts needless strain on the 
#' API. Therefore, you can save the "mart" object using the `saveRDS` function
#' for subsequent uses. You can later use `readRDS` to get the mart object
#' and give it to the `mart` argument of this function. This will then be
#' used instead of querying the API to generate the `mart` object each time.
#' 
#' @return 
#' Returns a `SummarizedExperiment` object containing the samples in the 
#' file given by `metadata_path`. This is loaded with `rowData`, including
#' the `gencode_id`, `ensembl_id` and `gene_id`. The metadata given as input
#' to this function is loaded into the `colData`. 

get_summarized_experiment <- function(
    counts, 
    metadata_path="scripts/20210713_salmon_counts_initial_look/w1_metadata.csv",
    mart=NULL
) {
    
    # get ensembl ids from gencode ids that come with the count matrix
    gencode_ids <- counts[,1]
    ensembl_ids <- sapply(strsplit(gencode_ids, ".", fixed=T), function(x) x[1])
    counts <- counts[,-1]
    
    # remove duplicate ensembl IDs
    counts <- counts[!duplicated(ensembl_ids),]
    gencode_ids <- gencode_ids[!duplicated(ensembl_ids)]
    ensembl_ids <- ensembl_ids[!duplicated(ensembl_ids)]
    rownames(counts) <- ensembl_ids
    
    # get gene ids from ensembl ids via biomart
    if (is.null(mart)) {
        mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    }
    
    ensembl_to_geneid <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=ensembl_ids, mart=mart)
    ensembl_to_geneid <- ensembl_to_geneid[!duplicated(ensembl_to_geneid$ensembl_gene_id),]
    
    count_row_names_joined <- dplyr::left_join(data.frame(gencode_id=gencode_ids, ensembl_id=ensembl_ids), ensembl_to_geneid, by = c("ensembl_id"="ensembl_gene_id"))
    colnames(count_row_names_joined) <- c("gencode_id", "ensembl_id", "gene_id")
    
    count_row_names_joined$gene_id[count_row_names_joined$gene_id == ""] <- count_row_names_joined$ensembl_id[count_row_names_joined$gene_id == ""]
    
    # incorrect sample naming
    colnames(counts) <- gsub("C133","C105", gsub("\\.", "_", gsub("CC", "C", gsub("HC", "", colnames(counts)))))
    
    # setup the metadata
    counts_coldata <- data.frame(fread(metadata_path))
    
    # subset the counts by the metadata
    counts <- counts[,colnames(counts) %in% counts_coldata$sample_id]
    counts_coldata <- counts_coldata[counts_coldata$sample_id %in% colnames(counts),]
    
    stopifnot(all(colnames(counts) %in% counts_coldata$sample_id))
    
    counts_coldata_joined <- dplyr::left_join(
        data.frame(sample_id = colnames(counts)),
        counts_coldata,
        by = c("sample_id" = "sample_id")
    )
    
    # create the final SE object
    se_object <- SummarizedExperiment(counts, colData = counts_coldata_joined, rowData = count_row_names_joined)
    assayNames(se_object) <- "counts"
    
    # create a new time column
    colnames(colData(se_object)) <- gsub("date_positive_swab", "date_first_positive_swab", colnames(colData(se_object)))
    
    se_object$time_from_first_positive_swab <- as.numeric(
        as.Date(se_object$sample_date, format = "%d/%m/%Y") -
        as.Date(se_object$date_first_positive_swab, format = "%d/%m/%Y")
    )
    
    se_object$time_from_first_symptoms <- as.numeric(
        as.Date(se_object$sample_date, format = "%d/%m/%Y") -
        as.Date(se_object$date_first_symptoms, format = "%d/%m/%Y")
    )
    
    se_object$date_first_x <- se_object$date_first_symptoms
    se_object$time_from_first_x <- se_object$time_from_first_symptoms
    
    for (i in 1:length(se_object$time_from_first_symptoms)) {
        
        if (is.na(se_object$time_from_first_symptoms[i])) {
            
            max_time <- "swab"
            
        } else if (is.na(is.na(se_object$time_from_first_positive_swab[i]))) {
            
            max_time <- "symptoms"
            
        } else {
            
            if (se_object$time_from_first_symptoms[i] > se_object$time_from_first_positive_swab[i]) {
                max_time <- "symptoms"
            } else {
                max_time <- "swab"
            }
        }
        
        if (max_time == "swab") {
            
            se_object$date_first_x[i] <- se_object$date_first_positive_swab[i]
            se_object$time_from_first_x[i] <- se_object$time_from_first_positive_swab[i]
            
        } else {
            se_object$date_first_x[i] <- se_object$date_first_symptoms[i]
            se_object$time_from_first_x[i] <- se_object$time_from_first_symptoms[i]
        }
    }
    
    # check we are using the maximum time
    stopifnot(!any(se_object$time_from_first_x < se_object$time_from_first_symptoms, na.rm = TRUE))
    stopifnot(!any(se_object$time_from_first_x < se_object$time_from_first_positive_swab, na.rm = TRUE))
    
    stopifnot(!any(is.na(se_object$time_from_first_x[!is.na(se_object$time_from_first_symptoms)])))
    stopifnot(!any(is.na(se_object$time_from_first_x[!is.na(se_object$time_from_first_positive_swab)])))
    
    return(se_object)
}

#' Normalise the HTSeq counts matrix
#' 
#' This can be used to do different levels of normalisation to the counts data.
#' First you must convert the counts data into a `SummarizedExperiment` object
#' via the `get_summarized_experiment` function.
#' 
#' To use this function, you will usually want to normalise your `SummarizedExperiment`
#' object using one of the `normalisation_level` arguments with other parameters as default.
#' However, there are additional arguments to this function that allow for further customisation
#' of the normalisation procedure.
#' 
#' @param se
#' This should be a SummarizedExperiment object containing row and column data
#' as generated by the `get_summarized_experiment`.
#' 
#' @param normalisation_level
#' This parameter is the main determinant of the normalisation that is carried out.
#' The other parameters to this function can usually be ignored, but care must 
#' be taken to specify this one correctly. It can take on the following values:
#' 
#'  - `pre_dream`         : This string should be given if you wish to carry out the normalisation
#'                          required prior to running the limma-voom-dream functions. Note that
#'                          this will lead to the generation of an `edgeR::DGEList` rather than
#'                          a `SummarizedExperiment`.
#'                          
#'  - `logcpm`            : If "logcpm" is specified, then the function will return the `SummarizedExperiment`
#'                          object after running the `edgeR::cpm function`. Note that if you also specify
#'                          `log_transform = FALSE`, then this takes precendence and the cpm transformed matrix 
#'                          will not actually be log transformed.
#'                          
#'  - `variance_filtered` : This option will take the cpm transformed counts and apply variance filtering
#'                          on top of this. By default, the top 1000 variable genes will be returned.
#'                          
#' @param group
#' Note that lowly expressed counts are removed prior to all of the above normalisation_level 
#' configurations (unless `filter_by_expr = FALSE` is specified). The `group` parameter will
#' be sent to the `group` argument of the `edgeR::filterByExpr` function. This is useful when
#' carrying out an analysis, such as differential expression, that is looking at a particular group.
#' For instance, when doing an analysis comparing the `case_control` variable, you should set
#' `group` to "case_control".
#' 
#' @param log_transform
#' By default, both the `logcpm` and `variance_filtered` values of `normalisation_level` will
#' lead to the generation of log transformed counts. However, if you set `log_transform` to `FALSE`
#' then counts will not be log transformed (i.e. this variable takes precedence). Note that
#' this parameter has no effect if you specify `normalisation_level` as `pre_dream`.
#' 
#' @param variance_filter_cutoff
#' This parameter allows you to choose how many genes you want to be retained after
#' filtering the top variable genes. If you give a value of `variance_filter_cutoff`
#' between 0 and 1, the top `variance_filter_cutoff` * 10 percentile is kept. Else,
#' you can supply a number above 1 to specify an absolute number of genes to be
#' retained (default 1000). Note that this parameter only has an effect if 
#' you specify `normalisation_level` as `variance_filtered`
#' 
#' @param variance_filter_method
#' This allows for the use of different variance filtering methods. This is passed
#' to the `method` argument of the `M3C::featurefilter` function.
#' 
#' @param filter_by_expr
#' For all of the `normalisation_level` options, lowly expressed genes are removed
#' by default using the `filterByExpr` function. However, if you set this parameter
#' to `FALSE`, then this step will be skipped entirely regardless of `normalisation_level`.
#' 
#' @return 
#' If `normalisation_level` is set to `pre_dream`, then an `edgeR::DGEList` object is
#' returned. Else, a `SummarizedExperiment` object is returned which contains two assays. 
#' Assay 1 is the original count data whereas assay 2 is the count data normalized
#' by whichever procedure you have specified.
#' 
#' Therefore, you should do `assay(se, 2)` to extract your raw counts data - NOT `assay(se)`

normalize_se <- function(
    se, 
    normalisation_level="pre_dream", 
    group = NULL, 
    log_transform = TRUE,
    variance_filter_cutoff=1000,
    variance_filter_method="MAD",
    filter_by_expr = TRUE
) {
    # keep highly expressed genes
    if (!is.null(group)) {
        group <- se[[group]]
        
        stopifnot(!is.null(group))
    }
    
    if (filter_by_expr) {
        fbe_keep <- filterByExpr(se, group=group)
        se <- se[fbe_keep,]
        
        message(paste0(sum(fbe_keep), " out of ", length(fbe_keep), " genes kept after filterByExpr"))
    }
    
    # prepare edgeR normalisation
    dge <- calcNormFactors(se)
    
    # this object is required as input to voomWithDreamWeights
    if (normalisation_level == "pre_dream") {
        return(dge)
    }
    
    # apply cpm and log transformation - does not overwrite the original counts
    assay(se, 2) <- cpm(dge, log = log_transform)
    assayNames(se) <- c("counts", "cpm")
    
    if (log_transform) {
        message("counts have been log transformed")
    }
    
    if (normalisation_level == "logcpm") {
        return(se)
    }
    
    # apply variance filter - if you specify a number rather than a percentige, this code converts it into a percentage
    if (variance_filter_cutoff > 1) {
        variance_filter_cutoff <- (variance_filter_cutoff / nrow(assay(se, 2))) * 100
    } else {
        variance_filter_cutoff <- variance_filter_cutoff * 100
    }
    
    # figure out which genes to keep based on variability (default is mean absolute deviation)
    ff_keep <- rownames(M3C::featurefilter(assay(se, 2), percentile = variance_filter_cutoff, method=variance_filter_method, topN = 5)[[1]])
    ff_keep <- rownames(se) %in% ff_keep
    se <- se[ff_keep,]
    
    message(paste0(sum(ff_keep), " out of ", length(ff_keep), " genes kept after featurefilter"))
    
    if (normalisation_level == "variance_filtered") {
        return(se)
    }
}

#' load soma data

get_soma_data <- function(soma_abundance, sample_meta, sample_technical_meta, feature_meta, ret_wide = FALSE) {
    
    soma_abundance <- read.table(soma_abundance, row.names = 1, sep = ",", header = TRUE)
    feature_meta <- read.csv(feature_meta)
    
    w_metadata <- data.frame(fread(sample_meta))
    sample_technical_meta <- read.csv(sample_technical_meta)
    
    w_metadata <- w_metadata[w_metadata$sample_id %in% sample_technical_meta$sample_id,]
    sample_technical_meta <- sample_technical_meta[sample_technical_meta$sample_id %in% w_metadata$sample_id,]
    
    combined_meta <- dplyr::left_join(sample_technical_meta, w_metadata, by = c("sample_id" = "sample_id", "individual_id" = "individual_id"))
    
    soma_abundance <- soma_abundance[rownames(soma_abundance) %in% combined_meta$sample_id,]
    
    rank_norm <- function(u, k = 0.375) {
        
        n <- length(u)
        r <- rank(u)
        
        return(qnorm((r - k) / (n - 2 * k + 1)))
    }
    
    for (i in 1:ncol(soma_abundance)) {
        
        soma_abundance[,i] <- rank_norm(soma_abundance[,i])
    }
    
    colnames(combined_meta) <- gsub("date_positive_swab", "date_first_positive_swab", colnames(combined_meta))
    
    combined_meta$time_from_first_positive_swab <- as.numeric(
        as.Date(combined_meta$sample_date, format = "%d/%m/%Y") -
            as.Date(combined_meta$date_first_positive_swab, format = "%d/%m/%Y")
    )
    
    combined_meta$time_from_first_symptoms <- as.numeric(
        as.Date(combined_meta$sample_date, format = "%d/%m/%Y") -
            as.Date(combined_meta$date_first_symptoms, format = "%d/%m/%Y")
    )
    
    combined_meta$date_first_x <- combined_meta$date_first_symptoms
    combined_meta$time_from_first_x <- combined_meta$time_from_first_symptoms
    
    for (i in 1:length(combined_meta$time_from_first_symptoms)) {
        
        if (is.na(combined_meta$time_from_first_symptoms[i])) {
            
            max_time <- "swab"
            
        } else if (is.na(is.na(combined_meta$time_from_first_positive_swab[i]))) {
            
            max_time <- "symptoms"
            
        } else {
            
            if (combined_meta$time_from_first_symptoms[i] > combined_meta$time_from_first_positive_swab[i]) {
                max_time <- "symptoms"
            } else {
                max_time <- "swab"
            }
        }
        
        if (max_time == "swab") {
            
            combined_meta$date_first_x[i] <- combined_meta$date_first_positive_swab[i]
            combined_meta$time_from_first_x[i] <- combined_meta$time_from_first_positive_swab[i]
            
        } else {
            combined_meta$date_first_x[i] <- combined_meta$date_first_symptoms[i]
            combined_meta$time_from_first_x[i] <- combined_meta$time_from_first_symptoms[i]
        }
    }

    # check we are using the maximum time
    stopifnot(!any(combined_meta$time_from_first_x < combined_meta$time_from_first_symptoms, na.rm = TRUE))
    stopifnot(!any(combined_meta$time_from_first_x < combined_meta$time_from_first_positive_swab, na.rm = TRUE))

    stopifnot(!any(is.na(combined_meta$time_from_first_x[!is.na(combined_meta$time_from_first_symptoms)])))
    stopifnot(!any(is.na(combined_meta$time_from_first_x[!is.na(combined_meta$time_from_first_positive_swab)])))
    
    if (ret_wide) {
        return(list(
            soma_abundance = soma_abundance,
            sample_meta = combined_meta,
            feature_meta = feature_meta
        ))
    }
    
    soma_abundance$sample_id <- rownames(soma_abundance)
    
    soma_abundance <- tidyr::pivot_longer(soma_abundance, -sample_id, names_to = "converted_seq_ids", values_to = "RFU")
    
    prev_nrow <- nrow(soma_abundance)
    soma_abundance <- dplyr::left_join(soma_abundance, combined_meta, by = c("sample_id" = "sample_id"))
    stopifnot(nrow(soma_abundance) == prev_nrow)
    
    prev_nrow <- nrow(soma_abundance)
    soma_abundance <- dplyr::left_join(soma_abundance, feature_meta, by = c("converted_seq_ids" = "converted_seq_ids"))
    stopifnot(nrow(soma_abundance) == prev_nrow)
    
    return(soma_abundance)
}
