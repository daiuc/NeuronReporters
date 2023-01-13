# Rerun network analysis.
# This code should be roughly an replica of the jupyter notebook I used back
# in 2019-2020 that generated the network interactions.
# This is the next step of `compute_interactions.R`
# This combine interactions from regulator -> target and target -> regulator
# procedures found interactions, and later integrate RNA-seq reads.


suppressMessages(library(stringr))
suppressMessages(library(purrr))
suppressMessages(library(readxl))
suppressMessages(library(data.table))

#-------------------------------------------------------------------------------
# INPUTS, OUTPUTS
#-------------------------------------------------------------------------------


# inputfiles = list(
#     interactions_raw = "results/reviews/interactions/Interactions_AllTF_intergenic_raw_20230110.xlsx",
#     interactions_norm = "results/reviews/interactions/Interactions_AllTF_intergenic_norm_20230110.xlsx",
#     interaction_sheets = c("Regulator_to_Target"),
#     dge = "resources/TFscreen/RNA-seq_timepoint_deseq_result_20200102.xlsx",
#     dgetimepoints = c("H15_vs_ES", "D1_vs_ES", "D4_vs_ES"),
#     rawcount = "rawCounts",
#     hitlist = "resources/TFscreen/Hitlist_20191230.csv", # has header
#     tflist = "resources/TFscreen/TFlist_20191230.csv" # has header
#     # RNAsamples = "resources/TFscreen/RNASeqSampleNames.csv"
# )

inputfiles = list(
    interactions_raw = snakemake@input[["interaction_raw"]],
    interactions_norm = snakemake@input[["interaction_norm"]],
    interaction_sheets = snakemake@params[["interaction_sheets"]],
    dge = snakemake@input[["dge"]],
    dgetimepoints = snakemake@params[["dge_timepoints_sheets"]],
    rawcount = snakemake@params[["dge_rawcount_sheet"]],
    hitlist = snakemake@input[["hitlist"]],
    tflist = snakemake@input[["tflist"]]
)


outputfiles = list(
    csv = snakemake@output[[1]],
    readme = snakemake@output[[2]]
)

cat("input files:\n")
print(inputfiles)
cat("parameters:\n")
print(snakemake@params)
cat("\noutput files:\n")
print(outputfiles)


#-------------------------------------------------------------------------------
# CUSTOM FUNCS
#-------------------------------------------------------------------------------

myflatten = function(l, use.prefix=TRUE) {
    #' @param l: is a nested list

    # use variable name as a string, use it as prefix
    prefix = deparse(substitute(l))
    l = flatten(l)
    if (use.prefix) names(l) = paste(prefix, names(l), sep=".")
    return(l)
}


readInteraction = function(excel_file, sheets) {
    #' @param excel_file: xlsx file of the interactions
    #' @param sheets: vector, sheet names in the xslx file
    
    #' @return a data.table of interactions, with columns:
    #'         interaction, regulator_gene, target_gene, ES, H1, H4, H16, D1, D4

    
    names(sheets) = sheets
    xlsx = map(sheets, ~read_excel(excel_file, sheet = .x) %>% as.data.table)

    # add an interaction column
    walk(xlsx, ~.x[, interaction := paste(regulator_gene, target_gene, sep="_")])

    if (length(xlsx) == 1) {
        interactions = xlsx[[1]]
        interactions = interactions[,
        .(interaction, regulator_gene, target_gene, ES, H1, H4, H16, D1, D4)]
    }
    
    return(interactions) # a data.table
}

#-------------------------------------------------------------------------------
# READ IN INTERACTIONS
#-------------------------------------------------------------------------------

# interactions computed with raw ATAC counts
inter_raw = readInteraction(inputfiles$interactions_raw, inputfiles$interaction_sheets)

# interactions computed with normalized ATAC counts
inter_norm = readInteraction(inputfiles$interactions_norm, inputfiles$interaction_sheets)



#-------------------------------------------------------------------------------
# READ IN RNA-SEQ REASULTS - READS & FDR & L2FC
#-------------------------------------------------------------------------------

# DESeq2 results that has l2fc and FDR info
sheets2 = inputfiles$dgetimepoints
if (all(sheets2 == c("H15_vs_ES", "D1_vs_ES", "D4_vs_ES"))) {
    dge = map(sheets2,
              ~read_excel(inputfiles$dge, sheet = .x) %>%
              as.data.table)
} else {
    stop(cat("Excel sheets not in correct order!"))
}

dge = map(dge,
~.x[!is.na(gene_name), .(gene_id, gene_name, log2FoldChange, pvalue, padj)])
names(dge) = factor(c("H16", "D1", "D4"), levels = c("H16", "D1", "D4"))


# gene expression raw counts
dge.rawCounts = read_excel(inputfiles$dge, sheet = inputfiles$rawcount) %>%
    as.data.table

# columns for gene_id and gene_name
idcols = dge.rawCounts[, .(gene_id, gene_name)]

# remove genes correspond to multiple or NA gene_ids
# same set for both rawCounts and normCounts
remove.gene.idx = idcols[
    , .(N_geneid = length(unique(gene_id))),
    by = .(gene_name)
    ][N_geneid > 1, gene_name]

# normalize RNA-seq counts - sqrt of reads per million
dge.normCounts = dge.rawCounts[, -c("gene_id", "gene_name")] %>%
    apply(2, function(x) sqrt(x * 1e6 / sum(x)) ) %>%
    as.data.table

# replicate names，each represents ES, H16, D1, D4
rna.rep1 <- c("S01_B1", "S15_B1", "S17_B1", "S07_B1")
rna.rep2 <- c("S01_B2", "S15_B2", "S17_B2", "S07_B2")
rna.rep3 <- c("S02_B1", "S16_B1", "S18_B1", "S08_B1")
rna.rep4 <- c("S02_B2", "S16_B2", "S18_B2", "S08_B2")

# combine 4 replicates into 1, using average
dge.normCounts <- pmap(
    list(dge.normCounts[, rna.rep1, with = F],
         dge.normCounts[, rna.rep2, with = F],
         dge.normCounts[, rna.rep3, with = F],
         dge.normCounts[, rna.rep4, with = F]),
    ~ pmap_dbl(list(..1, ..2, ..3, ..4), ~ mean(c(..1, ..2, ..3, ..4)))
    ) %>%
    as.data.table

# add back gene_id and gene_name columns, and rename the time points
names(dge.normCounts) = c("ES", "H16", "D1", "D4")
dge.normCounts[, `:=`(gene_id = idcols$gene_id, gene_name = idcols$gene_name)]
dge.normCounts = dge.normCounts[, .(gene_id, gene_name, ES, H16, D1, D4)]

# remove some genes with NA name or multiple gene_ids
dge.normCounts = dge.normCounts[
    !gene_name %in% remove.gene.idx & !is.na(gene_name)]

# nest list column
dge.normCounts = dge.normCounts[
    , .(ge = list(list(ES = ES, H16 = H16, D1 = D1, D4 = D4))),
    by = .(gene_id, gene_name)]

# to unnest list column run this
# dge.normCounts[, flatten(ge), by=.(gene_id, gene_name)]


# Do the same to the raw counts table
# combine 4 replicates into 1, using average
dge.rawCounts = pmap(
    list(dge.rawCounts[, rna.rep1, with = F],
         dge.rawCounts[, rna.rep2, with = F],
         dge.rawCounts[, rna.rep3, with = F],
         dge.rawCounts[, rna.rep4, with = F]),
    ~ pmap_dbl(list(..1, ..2, ..3, ..4), ~ round(mean(c(..1, ..2, ..3, ..4))))
    ) %>%
    as.data.table
names(dge.rawCounts) = c("ES", "H16", "D1", "D4")
dge.rawCounts[, `:=`(gene_id = idcols$gene_id, gene_name = idcols$gene_name)]
dge.rawCounts = dge.rawCounts[, .(gene_id, gene_name, ES, H16, D1, D4)]

# remove some genes with NA name or multiple gene_ids
dge.rawCounts = dge.rawCounts[
    !gene_name %in% remove.gene.idx & !is.na(gene_name)]

# nest list column
dge.rawCounts = dge.rawCounts[
    , .(cnt = list(list(ES=ES, H16=H16, D1=D1, D4=D4))),
    by = .(gene_id, gene_name)]


#-------------------------------------------------------------------------------
# INTEGRATE TARGET GENE EXPRESSION (NORM & RAW)
#-------------------------------------------------------------------------------


# nest interaction ATAC-seq normalized reads
inter_raw = inter_raw[
, .(atac.cnt = list(list(ES=round(ES), H16=round(H16), D1=round(D1), D4=round(D4)))),
    by = .(regulator_gene, target_gene)
]

# nest interaction ATAC-seq normalized reads
inter_norm = inter_norm[
    , .(atac = list(list(ES=ES, H16=H16, D1=D1, D4=D4))),
    by = .(regulator_gene, target_gene)
]

# combine raw and norm interaction reads, makesure all rows the same
if (all(inter_raw$regulator_gene == inter_norm$regulator_gene &
        inter_raw$target_gene == inter_norm$target_gene)) {
    inter_corr = cbind(inter_norm, inter_raw[, .(atac.cnt)])
} else {
    stop(cat("ID rows don't match in normalized vs. raw interactions!"))
}

# innner join with normalized counts to get targets' normalized gene expression
inter_corr = inter_corr[
    dge.normCounts[, .(target_gene = gene_name, tar.ge = ge)],
    on = "target_gene", nomatch = NULL
    ][order(regulator_gene, target_gene)]

# inner join with raw conts to get targets' raw gene exression, termed as `cnt`
inter_corr = inter_corr[
    dge.rawCounts[, .(target_gene = gene_name, tar.cnt = cnt)],
    on = "target_gene", nomatch = NULL
][order(regulator_gene, target_gene)]


#-------------------------------------------------------------------------------
# COR.TEST OF TARGET EXPRESSION (NORM) AND TARGET ATAC (NORM)
#-------------------------------------------------------------------------------

# Run correlation test on targets' normalized gene expression
# and normalized interaction reads on target genes --- this is
# essentially a motif's reads on the target gene
tictoc::tic() # about 1 minute
suppressWarnings(
    cor_test <- map2(inter_corr$atac, inter_corr$tar.ge,
                     ~ cor.test(unlist(.x), unlist(.y), method = "p")))
tictoc::toc()

# Extract pearson correlation and p-values
# store along with read counts for each interactions
inter_corr[, `:=`(
    corr = map_dbl(cor_test, ~.x$estimate),
    cor_p = map_dbl(cor_test, ~.x$p.value)
)]


#-------------------------------------------------------------------------------
# DESEQ RESULTS - LOG2 FC & FDR
# this step REMOVES interactions where regulators motifs are dimers
#-------------------------------------------------------------------------------

# first remove genes we don't want, same as raw and norm counts above
dge = map(dge, ~.x[!gene_name %in% remove.gene.idx & !is.na(gene_name)])

# extract log2 fold change (l2fc)
l2fc = imap(dge, ~.x[, .(gene_name, log2FoldChange, 
                     time = factor(.y, c("H16", "D1", "D4")))]) %>%
    rbindlist %>%
    dcast(gene_name ~ time, value.var = "log2FoldChange")
# nest l2fc column
l2fc = l2fc[, .(l2fc = list(list(H16 = H16, D1 = D1, D4 = D4))), by = gene_name]

# extract FDR
fdr = imap(dge, ~.x[, .(gene_name, fdr = padj,
                     time = factor(.y, c("H16", "D1", "D4")))]) %>%
    rbindlist  %>%
    dcast(gene_name ~ time, value.var = "fdr")
# nest fdr column
fdr = fdr[, .(fdr = list(list(H16 = H16, D1 = D1, D4 = D4))), by = gene_name]

# integrate target genes' l2fc and fdr
inter_corr = inter_corr[
    l2fc[, .(target_gene = gene_name, tar.l2fc = l2fc)],
    on = "target_gene", nomatch = NULL]

inter_corr = inter_corr[
    fdr[, .(target_gene = gene_name, tar.fdr = fdr)],
    on = "target_gene", nomatch = NULL]

# integrate regulator genes' l2fc and fdr
# Because some regulators are dimers such as ARNT::HIF1A,
# the result of an inner_join would remove these dimer (interactions)
inter_corr = inter_corr[
    l2fc[, .(regulator_gene = gene_name, reg.l2fc = l2fc)],
    on = "regulator_gene", nomatch = NULL]

inter_corr = inter_corr[
    fdr[, .(regulator_gene = gene_name, reg.fdr = fdr)],
    on = "regulator_gene", nomatch = NULL]


#-------------------------------------------------------------------------------
# INTEGRATE REGULATORS' GENE EXPRESSION (NORM & RAW)
#-------------------------------------------------------------------------------

inter_corr = inter_corr[
    dge.normCounts[, .(regulator_gene = gene_name, reg.ge = ge)],
    on = "regulator_gene", nomatch = NULL]

inter_corr = inter_corr[
    dge.rawCounts[, .(regulator_gene = gene_name, reg.cnt = cnt)],
    on = "regulator_gene", nomatch = NULL]



#-------------------------------------------------------------------------------
# LABEL ISHIT / ISTF
#-------------------------------------------------------------------------------

hitlist = fread(inputfiles$hitlist)
tflist = fread(inputfiles$tflist)

inter_corr[, `:=`(
    reg.isHIT = regulator_gene %in% hitlist$gene_name,
    reg.isTF = regulator_gene %in% tflist$gene_name,
    tar.isHIT = target_gene %in% hitlist$gene_name,
    tar.isTF = target_gene %in% tflist$gene_name
)]


#-------------------------------------------------------------------------------
# UNNEST LIST COLUMNS AND WRITE OUT
#-------------------------------------------------------------------------------

cat("Integrated all metrics for interactions, now unnest list columns and write results...\n")

tictoc::tic() # about 2 minutes

outDT = inter_corr[
    , c(list(corr = corr, cor_p = cor_p),
        # unnest list columns
        myflatten(atac),  myflatten(atac.cnt),
        myflatten(tar.ge), myflatten(tar.cnt),
        myflatten(tar.l2fc), myflatten(tar.fdr),
        myflatten(reg.ge), myflatten(reg.cnt),
        myflatten(reg.l2fc), myflatten(reg.fdr),
        # other non-list columns
        list(reg.isHit = reg.isHIT, reg.isTF = reg.isTF,
             tar.isHIT = tar.isHIT, tar.isTF = tar.isTF)
        ),
    by = .(regulator_gene, target_gene)]

tictoc::toc()

colsDescription = c(
    regulator_gene = "gene name of regulator",
    target_gene = "gene name of target",
    corr = "correlation",
    cor_p = "pvalue of correlation",
    atac = "sqrt(rpm) normalized reads in target gene's ATAC-seq peaks that match regulator motifs. suffix after 'atac.' indicates time point.",
    atac.cnt = "same as `atac`, except it is raw reads. suffix indicates time points",
    tar.ge = "target gene's expression level, sqrt(rpm) normalized",
    tar.cnt = "same as `tar.ge`, except raw counts",
    tar.l2fc = "target's log2FoldChange stated timepoint vs. ES",
    tar.fdr = "target's FDR based on DESeq2, stated timepoints vs. ES",
    reg.ge = "regulator gene's expression, sqrt(rpm) normalized",
    reg.cnt = "regulator gene's expression in raw counts",
    reg.l2fc = "regulator's log2FoldChange stated timepoint vs. ES",
    reg.fdr = "target's FDR based on DESeq2, stated time points vs. ES",
    reg.isHIT = "regulator is in crispr screen hits if TRUE",
    reg.isTF = "regulator is a TF if TRUE",
    tar.isHIT = "target is in crispr screen hits if TRUE",
    tar.isTF = "target is a TF if TRUE"
) %>% tibble::enframe(name = "column", value = "description")


fwrite(outDT, outputfiles$csv, sep = ",")
fwrite(colsDescription, outputfiles$readme, sep = " ")



#---- Saving previous code for having two sheets in excel file ---

# readInteraction = function(excel_file, sheets) {
#     #' @param excel_file: xlsx file of the interactions
#     #' @param sheets: vector, sheet names in the xslx file
    
#     #' @return a data.table of interactions, with columns:
#     #'         interaction, regulator_gene, target_gene, ES, H1, H4, H16, D1, D4

    
#     names(sheets) = sheets
#     xlsx = map(sheets, ~read_excel(excel_file, sheet = .x) %>% as.data.table)

#     # add an interaction column
#     walk(xlsx, ~.x[, interaction := paste(regulator_gene, target_gene, sep="_")])

#     # get keys, which is the interaction column
#     keys = map(xlsx, ~.x[, interaction])
#     shared.keys = intersect(keys$Regulator_to_Target, keys$Target_to_Regulator)
#     diff.keys = map(keys, ~ setdiff(.x, shared.keys))

#     # Combine the two datasets so that overlapping interactions are stored once
#     xlsx = map2(xlsx, diff.keys, ~.x[
#         interaction %in% c(.y, shared.keys),
#         .(interaction, regulator_gene, target_gene, ES, H1, H4, H16, D1, D4)
#     ])

#     interactions = rbindlist(xlsx)

#     # Because many interactions are repeated from both methods
#     # taking a mean by interaction
#     interactions = interactions[, .(
#         ES = mean(ES), H1 = mean(H1), H4 = mean(H4), 
#         H16 = mean(H16), D1 = mean(D1), D4 = mean(D4)),
#         by = .(interaction, regulator_gene, target_gene)
#     ]

#     return(interactions) # a data.table
# }







