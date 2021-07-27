get_summarised_experiments <- function(counts) {
    
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    
    # get ensembl ids from gencode ids that come with the count matrix
    gencode_ids <- counts[,1]
    ensembl_ids <- sapply(strsplit(gencode_ids, ".", fixed=T), function(x) x[1])
    counts <- counts[,-1]
    rownames(counts) <- gencode_ids
    
    # get gene ids from ensembl ids via biomart
    ensembl_to_geneid <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=ensembl_ids, mart=mart)
    ensembl_to_geneid <- ensembl_to_geneid[!duplicated(ensembl_to_geneid$ensembl_gene_id),]
    
    count_row_names_joined <- dplyr::left_join(data.frame(gencode_id=gencode_ids, ensembl_id=ensembl_ids), ensembl_to_geneid, by = c("ensembl_id"="ensembl_gene_id"))
    colnames(count_row_names_joined) <- c("gencode_id", "ensembl_id", "gene_id")
    
    # incorrect sample naming
    colnames(counts) <- gsub("C133","C105", gsub("\\.", "_", gsub("CC", "C", gsub("HC", "", colnames(counts)))))
    
    # setup the wave1 counts
    wave1_counts_coldata <- w1_metadata
    wave1_counts <- counts[,colnames(counts) %in% w1_metadata$sample_id]
    
    # setup the wave2 counts
    wave2_counts_coldata <- w2_metadata
    wave2_counts <- counts[,colnames(counts) %in% w2_metadata$sample_id]
    
    # set up the summarized experiments
    wave1_counts_se <- SummarizedExperiment(wave1_counts, colData = wave1_counts_coldata, rowData = count_row_names_joined)
    wave2_counts_se <- SummarizedExperiment(wave2_counts, colData = wave2_counts_coldata, rowData = count_row_names_joined)
    
    # return as a list
    return(list(wave1=wave1_counts_se, wave2=wave2_counts_se))
}