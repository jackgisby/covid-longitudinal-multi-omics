# split up data into summarized experiments for wave1 and wave2
get_summarized_experiment <- function(
    counts, 
    metadata_path="scripts/20210713_salmon_counts_initial_look/w1_metadata.csv"
) {
    
    # get ensembl ids from gencode ids that come with the count matrix
    gencode_ids <- counts[,1]
    ensembl_ids <- sapply(strsplit(gencode_ids, ".", fixed=T), function(x) x[1])
    counts <- counts[,-1]
    rownames(counts) <- gencode_ids
    
    # get gene ids from ensembl ids via biomart
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    
    ensembl_to_geneid <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=ensembl_ids, mart=mart)
    ensembl_to_geneid <- ensembl_to_geneid[!duplicated(ensembl_to_geneid$ensembl_gene_id),]
    
    count_row_names_joined <- dplyr::left_join(data.frame(gencode_id=gencode_ids, ensembl_id=ensembl_ids), ensembl_to_geneid, by = c("ensembl_id"="ensembl_gene_id"))
    colnames(count_row_names_joined) <- c("gencode_id", "ensembl_id", "gene_id")
    
    # incorrect sample naming
    colnames(counts) <- gsub("C133","C105", gsub("\\.", "_", gsub("CC", "C", gsub("HC", "", colnames(counts)))))
    
    # setup the metadata
    counts_coldata <- data.frame(fread(metadata_path))
    
    # subset the counts by the metadata
    counts <- counts[,colnames(counts) %in% w1_metadata$sample_id]
    
    # create the final SE object
    se_object <- SummarizedExperiment(counts, colData = counts_coldata, rowData = count_row_names_joined)
    assayNames(se_object) <- "counts"
    
    return(se_object)
}

normalize_se <- function(
    se, 
    normalisation_level="pre_dream", 
    group = NULL, 
    log_transform = TRUE,
    variance_filter_cutoff=1000,
    variance_filter_method="MAD"
) {
    
    # keep highly expressed genes
    fbe_keep <- filterByExpr(assay(se), group=group)
    se <- se[fbe_keep,]
    
    message(paste0(sum(fbe_keep), " out of ", length(fbe_keep), " genes kept after filterByExpr"))
    
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
    
    # apply variance filter - if you specify a number rather than a percentige, this code converts it into a percentage
    if (variance_filter_cutoff > 1) {
        variance_filter_cutoff <- (variance_filter_cutoff / nrow(assay(se, 2))) * 100
    }
    
    # figure out which genes to keep based on variability (default is mean absolute deviation)
    ff_keep <- rownames(M3C::featurefilter(assay(se, 2), percentile = variance_filter_cutoff, method=variance_filter_method, topN = 5)[[1]])
    ff_keep <- rownames(se) %in% ff_keep
    se <- se[ff_keep,]
    
    message(paste0(sum(ff_keep), " out of ", length(ff_keep), " genes kept after featurefilter"))
    
    return(se)
}
