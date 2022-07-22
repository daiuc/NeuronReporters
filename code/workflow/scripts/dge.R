# differential gene expression analysis
# replicating the analysis previously done using jupyter notebooks

library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(WriteXLS)
library(cowplot)
library(gridExtra)
library(EnhancedVolcano)

RUN_MODE = 2

if (RUN_MODE == 1) {
    input.file = snakemake@input[[0]]
    col.meta = snakemake@params[["COL_DATA"]]
    gene.lookup = snakemake@params[["GENE_LOOKUP"]]
    fig.prefix = snakemake@params[["FIG_PREFIX"]]
} else {
    input.file = "Results/RNAseq/RSEM/combined_RSEM_counts.txt"
    col.meta = "config/rna.tsv"
    gene.lookup = 'Data/Annotations/gencode_refseq_partial_modified_geneNames_ID_20191230.csv'
    fig.prefix = 'Results/Figs'
}




# Read in all RSEM counts -------------------------------------------------


# all counts
cnts = read_tsv(input.file)
cnts = column_to_rownames(cnts, "gene_id")

# all sample col meta
colMeta = read_tsv(col.meta)
colMeta = mutate(colMeta, key = paste(treatment, timepoint, rep, batch, sep = "_"))

# gene-gene-id lookup
lookup = read_csv(gene.lookup, col_types = 'cc') %>% 
    mutate_at('gene_id', ~str_extract(.x, "ENSG[0-9]+"))



# compare rep - batch consistancy -----------------------------------------

# for each timepoint, compare rep and batch

# compare rep1 and rep2 in batch A

plotScatterByRepBatch = function(cnts, colNames, imgType=".pdf", prefix = fig.prefix) {
    df = cnts[, colNames]
    X = colNames[1]
    Ys = colNames[2:3]
    TITLEs = map_chr(Ys, ~paste0(X, " VS ", .x))
    FN = paste0(prefix, "/dge_replicateQC_scatter_", X, imgType)
    
    plotBase = function(X, Y) {
        gg = ggplot(cnts) + geom_point(aes_string(as.symbol(X), as.symbol(Y)), alpha = .4) +
            geom_abline(slope = 1, color='firebrick') 
        return(gg)
    }
    
    addstyle = function(gg, TITLE) {
        gg = gg + labs(title = TITLE) +
            xlim(0, 5e4) + ylim(0, 5e4)
            theme_cowplot()
        return(gg)
    }
    
    ggs = map(1:2, ~plotBase(X, Ys[.x]) %>% addstyle(TITLE = TITLEs[.x]))
    
    ggsave(FN, grid.arrange(grobs= ggs, nrow=1), width = 10, height = 8)

}


# plot all comparisons
Scatter_cols = filter(colMeta, rep== 1 & batch == 'B') %>% pull(key)
Scatter_cols = map(Scatter_cols, 
                   ~c(X=.x, 
                         Y1=str_replace(.x, "_[1]_", "_2_"), 
                         Y2=str_replace(.x, "_[1]_B", "_2_A")))

# save replicate and batch comparing scatter plots
walk(Scatter_cols, ~plotScatterByRepBatch(cnts, .x, prefix = fig.prefix))





# Differential analysis on wt ony timepoints ------------------------------


# prep DESeq2 model -------------------------------------------------------



# first subset counts matrix and corresponding column meta data to only wild types
cnts.wtOnly = select(cnts, starts_with("wt"))
colMeta.wtOnly = left_join(data.frame(key = colnames(cnts.wtOnly)), colMeta, by = "key") %>% 
    select(-sample, -fastq)


# construct DESeq2 object with timepoint as contrast
dds.wtOnly = DESeqDataSetFromMatrix(cnts.wtOnly, colMeta.wtOnly, design = ~timepoint)

# add in gene name
rowData(dds.wtOnly) = left_join(data.frame(gene_id = rownames(dds.wtOnly)),
                       lookup, by = "gene_id")

# remove low counts genes
keep.rows = rowSums(counts(dds.wtOnly)) >= 10
dds.wtOnly = dds.wtOnly[keep.rows, ]

# build model - this may take a minute
model.wtOnly = DESeq(dds.wtOnly)

# check all comparisons available
resultsNames(model.wtOnly)





# Perform DGE and plot volcano plots --------------------------------------

# build funciton, with option to plot volcano and/or return dge dataframe
performDGE = function(dds, model, contrast = c('timepoint', 't1', 't0'), 
                      time_labels = c('H15', 'ES'), treatment = "WT", 
                      prefix, imgType=".pdf", returnDF=T, plotfig=F) {
    
    results = results(model, contrast = contrast)
    volcano.plot.fn = paste0(prefix, '/dge_vocano_', time_labels[1], '_vs_', time_labels[2], imgType)
    volcano.title = paste0(treatment, ": ", time_labels[1], '_vs_', time_labels[2])
    
    df = cbind(rowData(dds), data.frame(results@listData))
    EnhancedVolcano(df, 
                    lab = df$gene_name, 
                    x = 'log2FoldChange', 
                    y = 'padj', 
                    xlim=c(-8,8), 
                    ylim=c(0,150),
                    xlab = bquote(~Log[2]~ 'fold change'), 
                    ylab= bquote(~-Log[10]~adjusted~italic(P)),
                    legendPosition='bottom',
                    legendLabels=c('NS',bquote(~Log[2]~ 'FC'),'Adj. P', bquote('Adj P & '~Log[2]~ 'FC')),
                    legendLabSize = 8,
                    legendIconSize = 3.0, 
                    title = volcano.title, 
                    subtitle="", 
                    pCutoff = 0.05, 
                    FCcutoff = 2,
                    colAlpha = .4,
                    pointSize = 1.5,
                    labSize = 3.0
                    )
    ggsave(volcano.plot.fn, width = 10, height = 8)
    
    if (returnDF) {
        return(df)
    }
}


# run for each contrast

# t1 vs t0
performDGE(dds.wtOnly, model.wtOnly, 
           c("timepoint", "t1", "t0"), c("H15", "ES"), "WT", 
           fig.prefix, returnDF = F, plotfig = T)

# t2 vs t0
performDGE(dds.wtOnly, model.wtOnly, 
           c("timepoint", "t2", "t0"), c("Day1", "ES"), "WT", 
           fig.prefix, returnDF = F, plotfig = T)

# t3 vs t0
performDGE(dds.wtOnly, model.wtOnly, 
           c("timepoint", "t3", "t0"), c("Day2", "ES"), "WT", 
           fig.prefix, returnDF = F, plotfig = T)

# t4 vs t0
performDGE(dds.wtOnly, model.wtOnly, 
           c("timepoint", "t4", "t0"), c("Day4", "ES"), "WT", 
           fig.prefix, returnDF = F, plotfig = T)

# t5 vs t0
performDGE(dds.wtOnly, model.wtOnly, 
           c("timepoint", "t5", "t0"), c("Day7", "ES"), "WT", 
           fig.prefix, returnDF = F, plotfig = T)
