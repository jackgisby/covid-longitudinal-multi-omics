# run differential expression with dream
run_dream <- function(f, se, L, group="case_control", debug=FALSE) {
    
    # subset to expressed genes and calcNormFactors for downstream
    dge <- DGEList(assay(se))
    
    # keep highly expressed genes
    keep <- filterByExpr(dge, group = se[[group]])
    dge <- dge[keep,,keep.lib.sizes=FALSE]
    
    # prepare edgeR normalisation
    dge <- calcNormFactors(dge)
    
    if (debug & nrow(dge) > 500) {
        dge <- dge[1:500,]
    }
    
    # transform to log cpm, estimate mean variance relationship, prepare data for mixed modelling
    vobjDream <- voomWithDreamWeights(dge, f, data.frame(colData(se)), plot=TRUE)
    
    if (typeof(L) == "list") {
        
        coefs_to_compare <- L
        coef_names <- c()
        
        for (L_name in names(coefs_to_compare)) {
            
            coefs_to_compare[[L_name]] <- getContrast(vobjDream, f, data.frame(colData(se)), coefs_to_compare[[L_name]])
            
            if (typeof(L) == "list") {
                L <- coefs_to_compare[[L_name]]
            } else {
                L <- cbind(L, coefs_to_compare[[L_name]])
            }
            
            coef_names <- c(coef_names, L_name)
        }
        
        colnames(L) <- coef_names
        contrast <- NULL
        
        fitmm <- dream(vobjDream, f, data.frame(colData(se)), L)
        model_and_table <- list("fit" = fitmm)
        
        for (contrast in colnames(L)) {
            
            tt <- topTable(fitmm, coef=contrast, number = nrow(dge))
            tt$gencode_id <- rownames(tt)
            tt <- dplyr::left_join(tt, data.frame(rowData(se)))
            
            model_and_table[[paste0(contrast, "_tt")]] <- tt
        }
        
    } else {
        contrast <- L
        fitmm <- dream(vobjDream, f, data.frame(colData(se)))
        model_and_table <- list("fit" = fitmm)
        
        tt <- topTable(fitmm, coef=contrast, number = nrow(dge))
        tt$gencode_id <- rownames(tt)
        tt <- dplyr::left_join(tt, data.frame(rowData(se)))
        
        model_and_table[["tt"]] <- tt
    }
    
    return(model_and_table)
}

# run GSEA based on differential expression results
run_topgo <- function(tt, de_p_cutoff=0.01, go_p_cutoff=0.05, desc="case_control") {
    
    # this named vector is required to create a topGO object
    geneList <- tt$P.Value
    names(geneList) <- tt$ensembl_id
    
    # get raw p value cutoff from adjusted p values
    p_val_cutoff <- max(tt$P.Value[tt$adj.P.Val <= de_p_cutoff])
    
    # this selects the significant genes in the topgo object
    get_sig_genes <- function (allScore) {
        return(allScore < p_val_cutoff)
    }
    
    # create the topgo object
    w1_case_control_tg <- new(
        "topGOdata", 
        description = desc,
        ontology = "BP",
        allGenes = geneList,
        annotationFun = annFUN.org,
        ID = "ensembl",
        geneSelectionFun = get_sig_genes,
        mapping = "org.Hs.eg.db"
    )
    
    # run topGO GSEA
    w1_case_control_terms <- runTest(w1_case_control_tg, "weight01", "fisher")
    
    # extract topGO results
    w1_case_control_df <- data.frame(w1_case_control_terms@score[w1_case_control_terms@score < go_p_cutoff])
    colnames(w1_case_control_df) <- "p_value"
    
    w1_case_control_df$go_term <- rownames(w1_case_control_df)
    w1_case_control_df$term_description <- Term(w1_case_control_df$go_term)
    
    return(w1_case_control_df)
}

# create volcano plots from differential expression data
tt_volcano_plot <- function(tt, fc_col="logFC", de_p_cutoff=0.01, n_pos_fc=8, n_neg_fc=n_pos_fc, n_pos_pval=n_pos_fc, n_neg_pval=n_pos_pval) {
    
    # where there is no gene ID use ensembl ID
    tt$gene_id[tt$gene_id == ""] <- tt$ensembl_id[tt$gene_id == ""]
    
    # prepare colouring
    tt$de_col <- mapply(tt$adj.P.Val, tt[[fc_col]], FUN = function(adj_p, fc) {
        if (adj_p < de_p_cutoff) {
            if (fc > 0) {
                return("upreg")
            } else {
                return("downreg")
            }
        } else {
            return("non_sig")
        }
    })
    
    # label the top up and downregulated genes
    top_pos_fc <- Rfast::nth(tt[[fc_col]], n_pos_fc, descending = TRUE)
    top_neg_fc <- Rfast::nth(tt[[fc_col]], n_neg_fc, descending = FALSE)
    top_pos_pval <- Rfast::nth(tt$P.Value[tt[[fc_col]] > 0], n_pos_pval, descending = FALSE)
    top_neg_pval <- Rfast::nth(tt$P.Value[tt[[fc_col]] < 0], n_neg_pval, descending = FALSE)
    
    tt_to_label <- tt$adj.P.Val <= de_p_cutoff & (
        tt[[fc_col]] >= top_pos_fc |
            tt[[fc_col]] <= top_neg_fc |
            (tt$P.Value <= top_pos_pval & tt[[fc_col]] > 0) |
            (tt$P.Value <= top_neg_pval & tt[[fc_col]] < 0))
    
    tt$logp <- -log10(tt$P.Value)
    
    # make the ggplot
    volcano_plot <- ggplot(tt, aes_string(fc_col, "logp", col="de_col")) +
        geom_point() +
        scale_color_manual(values=c("#2C7BB6", "black", "#D7191C")) +
        geom_text_repel(aes(label=gene_id), size = 2.25, color = "black", data = subset(tt, tt_to_label)) +
        theme(legend.position = "none")
    
    return(volcano_plot)
}

single_lmer <- function(data, formula_string) {
    
    out.model <- tryCatch(
        lmerTest::lmer(
            as.formula(formula_string),
            data=data,
            REML=TRUE,
            control = lmerControl(check.conv.singular = "ignore"
        )),
        warning = function(w){
            return(lmerTest::lmer(
                as.formula(formula_string),
                data=data,
                REML=TRUE,
                control=lmerControl(optimizer = "Nelder_Mead", check.conv.singular = "ignore")
            ))
        }
    )
    
    
    if(class(out.model) == "lmerModLmerTest") {
        return(out.model)
        
    } else {
        stop("Convergence issue not caught by single_lmer")
    }
}
