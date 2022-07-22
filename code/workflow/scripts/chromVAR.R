library(tidyverse)
library(chromVAR)
library(GenomicRanges)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(motifmatchr)
set.seed(2019)
setwd("/gpfs/commons/groups/sanjana_lab/cdai/TFscreen/atac/")
getwd()

# Samples
samplesheet <- read.csv('samplesheet2.csv', stringsAsFactors = F)

# Peaks
peak.file <- "/c/groups/sanjana_lab/cdai/TFscreen/atac/diffbind/diffbind_consensu_min2overlap.bed"
peaks <- getPeaks(peak.file, sort_peaks = T)
peaks <- resize(peaks, width = 500, fix = "center")

# ATACseq counts
my_counts_matrix <- read.table("diffbind/diffbind_consensu_min2overlap_readcounts.txt", header = T) %>% as.matrix
depth <- apply(my_counts_matrix, 2, sum) # read depth per sample library
colData <- column_to_rownames(samplesheet, 'SampleID')
colData <- add_column(colData, depth) # add read depth per sample

# convert raw counts to SummarizedExperiment object
fragment_counts <- SummarizedExperiment(assays = list(counts = my_counts_matrix),
                                        rowRanges = peaks, colData = colData)

fragment_counts.addGC <- addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38)

# each peak must have at least 10 reads across all 12 samples to be included
fragment_counts.filtered <- filterPeaks(fragment_counts.addGC, min_fragments_per_peak = 10, 
                                        non_overlapping = T)

#make a grange object of normalized atac-seq counts
my_counts_matrix.norm <- apply(my_counts_matrix, MARGIN = 2, function(x) x * 1000000/sum(x))
my_counts_matrix.norm <- peaks %>% as.data.frame %>% cbind(my_counts_matrix.norm) %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)



# Jaspar & annotations -------------------------------------------------------------------

jaspar_2020 <- readJASPARMatrix("JASPAR2020_combined_matrices_20191030.txt", matrixClass = "PFM")

#120 hit list
hitlist <- read.csv('/c/groups/sanjana_lab/cdai/TFscreen/Hitlist_20191230.csv', 
                    header = T, stringsAsFactors = F) %>% pull(hgnc_symbol)
#TF list
tflist <- read.csv('/c/groups/sanjana_lab/cdai/TFscreen/TFlist_20191230.csv', 
                   stringsAsFactors = F) %>% pull(hgnc_symbol)

# construct base motif_lookup table
motif_lookup <- TFBSTools::name(jaspar_2020) %>% 
    data.frame(stringsAsFactors = F) %>% 
    rownames_to_column("motif") %>% 
    mutate(gene_name=str_extract(`.`, "[a-zA-Z:0-9\\-]+")) %>%
    select(motif, gene_name)

# function to check if a tf name is in the substring of jaspar name
mystrfunc <- function(stringA, target) {
    # string A is a possible gene name, scalar
    # target is a vector of proper gene names
    # match if any target is a substring in string A, if so return True
    stringA <- rep(stringA, length(target))
    test <- unlist(map2(stringA, target, ~ str_detect(.x, .y)))
    return(any(test))
}

#motif_lookup <- motif_lookup  %>% mutate(is_hit = map(gene_name, ~ if_else(mystrfunc(.x, hitlist), "Yes", "No"))) %>%
#mutate(is_tf = map(gene_name, ~ if_else(mystrfunc(.x, tflist), "Yes", "No"))) 

is_hit <- map(motif_lookup$gene_name, ~ if_else(mystrfunc(.x, hitlist), "Yes", "No")) %>% unlist
is_tf <- map(motif_lookup$gene_name, ~ if_else(mystrfunc(.x, tflist), "Yes", "No")) %>% unlist

motif_lookup <- motif_lookup %>% add_column(is_hit = is_hit, is_tf = is_tf)


# Compute deviations ------------------------------------------------------

# get motif matches with peaks
motif_ix <- matchMotifs(jaspar_2020, fragment_counts.filtered, genome = BSgenome.Hsapiens.UCSC.hg38)

# computing deviations
dev <- computeDeviations(object = fragment_counts.filtered, annotations = motif_ix)

# compute differentialDiviations
difdev = differentialDeviations(dev, groups = "Condition")

# Add motif's gene names and whethere it's a hit / tf or not
# some motifs has 2 or 3 variants with the same gene name, thus taking the mean of all motif scores for the same gene
# Note deviationScores return the Zscore of deviation, while deviations return deviations
dev.scores <- deviationScores(dev) %>% 
    as.data.frame %>% 
    rownames_to_column('motif') %>% 
    left_join(motif_lookup, by = 'motif') %>%
    select(A1:is_tf) %>% 
    group_by(gene_name, is_hit, is_tf) %>% 
    summarise_all(mean) %>%
    ungroup()


# compute differential deviation --------------------------------

# 1H vs. ESC
NEUROG.motifs = c("NEUROG1","NEUROG2", "NEUROG2(var.2)")
dev[NEUROG.motifs, dev$Condition %in% c("ESC", "1H")] %>% 
    differentialDeviations(groups = "Condition")

# 4H vs. ESC
NEUROG.motifs = c("NEUROG1","NEUROG2", "NEUROG2(var.2)")
dev[NEUROG.motifs, dev$Condition %in% c("ESC", "4H")] %>% 
    differentialDeviations(groups = "Condition")

# 16H vs. ESC
NEUROG.motifs = c("NEUROG1","NEUROG2", "NEUROG2(var.2)")
dev[NEUROG.motifs, dev$Condition %in% c("ESC", "16H")] %>% 
    differentialDeviations(groups = "Condition")

# 24H vs. ESC
NEUROG.motifs = c("NEUROG1","NEUROG2", "NEUROG2(var.2)")
dev[NEUROG.motifs, dev$Condition %in% c("ESC", "24H")] %>% 
    differentialDeviations(groups = "Condition")

# 5D vs. ESC
NEUROG.motifs = c("NEUROG1","NEUROG2", "NEUROG2(var.2)")
dev[NEUROG.motifs, dev$Condition %in% c("ESC", "5D")] %>% 
    differentialDeviations(groups = "Condition")
