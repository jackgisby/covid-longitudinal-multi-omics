library(data.table)
library(SomaDataIO)
library(org.Hs.eg.db)

# use SomaDataIO to read normalised SOMAScan data
soma_wide <- read_adat("data/20210812_SOMA/ICL-216787_v4.1_SLEplasma_20210728/SS-216787_v4.1_SLEplasma.hybNorm.medNormInt.plateScale.calibrate.medNormSMP.adat")

# non-COVID samples were run in the same batch, so we remove these
soma_wide <- soma_wide[soma_wide$SampleType == "Sample",]
soma_wide <- soma_wide[substr(soma_wide$SampleId, 1, 1) == "C",]

# 242 covid samples pre filtering

print(sum(soma_wide$RowCheck == "FLAG"))  # 1
soma_wide <- soma_wide[soma_wide$RowCheck == "PASS",]

# 241 covid samples post filtering

names(attributes(soma_wide))
attr(soma_wide, "Col.Meta")
attr(soma_wide, "row.meta")

# extract metadata for samples
sample_meta <- data.frame(soma_wide[, !grepl("seq.", colnames(soma_wide))])
sample_meta <- sample_meta[, apply(sample_meta, 2, function(col) {return(!all(is.na(col)))})]

# extract and write sample metadata
sample_meta_select <- dplyr::select(sample_meta, PlateId, PlateRunDate, ScannerID, PlatePosition, SlideId, Subarray, SampleId, SubjectID, RowCheck)
colnames(sample_meta_select) <- c("plate_id", "plate_run_date", "scanner_id", "plate_position", "slide_id", "sub_array", "sample_id", "individual_id", "row_check")

write.csv(sample_meta_select, "data/20210819_soma_processed/sample_meta.csv", row.names = FALSE)

# extract feature metadata
feature_meta <- attr(soma_wide, "Col.Meta")

feature_meta <- feature_meta[feature_meta$Organism == "Human",]
feature_meta <- feature_meta[feature_meta$Type == "Protein",]

feature_meta$UniProtSingle <- gsub(".*\\|", "", feature_meta$UniProt)
feature_meta$EntrezGeneIDSingle <- gsub(".*\\|", "", feature_meta$EntrezGeneID)

# get gene IDs
ensembl_hg19 <- biomaRt::useMart(biomart= "ENSEMBL_MART_ENSEMBL",
                                 dataset="hsapiens_gene_ensembl")

uniprot_to_hgnc <- biomaRt::getBM(attributes = c('uniprotswissprot', 'hgnc_symbol'), 
                                  filters = 'uniprotswissprot', 
                                  values = unique(feature_meta$UniProtSingle), 
                                  mart = ensembl_hg19)

uniprot_to_hgnc <- uniprot_to_hgnc[!duplicated(uniprot_to_hgnc$uniprotswissprot),]
colnames(uniprot_to_hgnc) <- c("UniProtSingle", "gene_id")

feature_meta <- dplyr::left_join(feature_meta, uniprot_to_hgnc, by=c("UniProtSingle"="UniProtSingle"))

hs <- org.Hs.eg.db

entrez_to_gene_id <- select(hs, 
                            keys = feature_meta$EntrezGeneIDSingle,
                            columns = c("ENTREZID", "SYMBOL"),
                            keytype = "ENTREZID")

colnames(entrez_to_gene_id) <- c("EntrezGeneIDSingle", "gene_id")
entrez_to_gene_id <- entrez_to_gene_id[!duplicated(entrez_to_gene_id$EntrezGeneIDSingle),]

for (i in 1:nrow(feature_meta)) {
    
    if (is.na(feature_meta$gene_id[i]) | feature_meta$gene_id[i] == "") {
        
        if (feature_meta$EntrezGeneIDSingle[i] %in% entrez_to_gene_id$EntrezGeneIDSingle) {
            feature_meta$gene_id[i] <- entrez_to_gene_id$gene_id[entrez_to_gene_id$EntrezGeneIDSingle == feature_meta$EntrezGeneIDSingle[i]]
        }
    }
}

feature_meta$gene_id[feature_meta$UniProtSingle == "Q9BYB0"] <- "SHANK3"
feature_meta$gene_id[feature_meta$UniProtSingle == "Q6F5E7"] <- "TXNRD3NB"
feature_meta$gene_id[feature_meta$UniProtSingle == "P01880"] <- "IGHD"

feature_meta$converted_seq_ids <- paste0("seq.", gsub("-", ".", feature_meta$SeqId))

# 7240 pass colcheck
feature_meta <- feature_meta[feature_meta$ColCheck == "PASS",]

rownames(soma_wide) <- soma_wide$SampleId
soma_wide_matrix <- data.frame(soma_wide[, colnames(soma_wide) %in% feature_meta$converted_seq_ids])

stopifnot(all(colnames(soma_wide_matrix) %in% feature_meta$converted_seq_ids))
stopifnot(all(feature_meta$converted_seq_ids %in% colnames(soma_wide_matrix)))

# write abundances
write.csv(soma_wide_matrix, "data/20210819_soma_processed/soma_abundances.csv", row.names = TRUE)

# write feature metadata
feature_meta_select <- dplyr::select(feature_meta, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, Dilution, ColCheck, UniProtSingle, EntrezGeneIDSingle, gene_id, converted_seq_ids)
colnames(feature_meta_select) <- c("soma_id", "target_full_name", "target", "uniprot", "entrez_gene_id", "entrez_gene_symbol", "dilution", "col_check", "uniprot_single", "entrez_gene_id_single", "gene_id", "converted_seq_ids")

write.csv(feature_meta_select, "data/20210819_soma_processed/feature_meta.csv", row.names = FALSE)
