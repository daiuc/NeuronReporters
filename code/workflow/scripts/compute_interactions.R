# Rerun network analysis.
# This code should be roughly an replica of the jupyter notebook I used back
# in 2019-2020 that generated the network interactions


suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(TFBSTools))
suppressMessages(library(chromVAR))
suppressMessages(library(motifmatchr))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(BiocParallel))
suppressMessages(library(parallel))


NCore = min(parallel::detectCores(), snakemake@threads[[1]])
cat(paste0("Register ", NCore, " cores.\n"))
register(MulticoreParam(NCore))

# register(MulticoreParam(16))

# if (getwd() != "/gpfs/commons/groups/sanjana_lab/cdai/NeuronReporters/code") {
#   setwd("/gpfs/commons/groups/sanjana_lab/cdai/NeuronReporters/code")
# }

set.seed(2019)


#-------------------------------------------------------------------------------
# INPUT SETUP
#------------------------------------------------------------------------------


# inputfiles = list(
#   samplesheet = "resources/TFscreen/atac/samplesheet2.csv", # has header
#   raw_counts = "resources/TFscreen/atac/diffbind/diffbind_consensu_min2overlap_NO_RECENTER.txt", 
#   hitlist = "resources/TFscreen/Hitlist_20191230.csv", # has header
#   tflist = "resources/TFscreen/TFlist_20191230.csv", # has header
#   jaspar_2020 = "resources/TFscreen/atac/JASPAR2020_combined_matrices_20191030.txt",
#   gene_region = "resources/TFscreen/atac/Protein_coding_genes_Up_100k_20230111.bed",
#   normalize = "norm" # either `norm` or `raw`
# )

inputfiles = list(
  samplesheet = snakemake@input[["samplesheet"]], # has header
  raw_counts = snakemake@input[["rawcounts"]], # has header
  hitlist = snakemake@input[["hitlist"]], # has header
  tflist = snakemake@input[["tflist"]], # has header
  jaspar_2020 = snakemake@input[["jaspar2020"]],
  gene_region = snakemake@input[["generegion"]],
  normalize = snakemake@wildcards[["Normalize"]] # either `norm` or `raw`
)

outputfiles = list(
  interactions = snakemake@output[[1]]
)

norm_flag = inputfiles$normalize == "norm"

cat("input files:\n")
print(inputfiles)
cat("\noutput files:\n")
print(outputfiles)

#-------------------------------------------------------------------------------
# CUSTOM FUNCTIONS
#-------------------------------------------------------------------------------

mystrfunc = function(motifnames, target_match) {
  #' @param motifnames: a character vector of motif names
  #' @param target_match: string vector, e.g. hitlist, or tflist
  #' @returns : equal length vector of TRUE/FALSE indicating whether or not
  #'            it's in target

  x = str_match_all(motifnames, "[A-Z0-9\\-]{2,13}")
  in_list = map_lgl(x, function(m) {
    # each m is a character matrix, col1 is the full match
    motifs = as.vector(m[,1]) # clean gene names of TF motifs
    is.in.list = map_lgl(motifs, ~ .x %in% target_match)
    # usually 1 gene, but in the case of dimers, multiple gene name exists
    # return TRUE as long as one gene is in target list
    return(any(is.in.list))
  })
  return(if_else(in_list, "Yes", "No"))
}


getTargetsFromMotif <- function(motif_name,
                                motif.matches,
                                genomic.features,
                                min.overlap = 10,
                                genelist) {
    #' Given a motif name (JASPAR), find its binding targets with read counts,
    #' observed in ATAC-seq
    
    #' @param motif_name: string, this is the jaspar motif name, not gene name
    #' @param motif.matches: matchMotif object, e.g. motifMatch() result of peaks
    #' and jaspar 2020 motifs
    #' @param genomic.features: GenomicRanges object, e.g. a bed formated promoter
    #' region coordinates, converted into GRanges object
    #' @param min.overlap: integer, minimum overlap in base pairs
    #' @param genelist: vector, a list of gene_names. Note these gene names are some
    #' times not the same as motif names, hence the need of motif_lookup
    
    #' --------------- OUTPUT -----------------
    #' @return overlapped.target.readcount: dataframe, containing read counts of
    #' genes that are targets of a given regulator. Read counts sum of all peaks
    #' that match to the regulator's motif.

    #' Get motif matches matrix (rows are peak coordianes, columns are each motif
    #' name, values are logical values indicating match)
    
    match.matrix <- motifMatches(motif.matches)
    
    # Get peaks that have binding sites matching given motif,
    # each GRange also include read count columns
    GRanges.index <- which(match.matrix[ , motif_name]) 
    matched.target.GRanges <- rowRanges(motif.matches)[GRanges.index]
    
    # Find overlaps between:
    # 1. peaks that match binding sites of a given motif, and
    # 2. coordinates of gene annotation
    target.intersect.feature <- GenomicRanges::findOverlaps(
      query = matched.target.GRanges,
      subject = genomic.features,
      minoverlap = min.overlap,
      ignore.strand = T)

    # Get the genes (from annoation) that overlap with peaks 
    # (that have binding sites of given motif)
    overlapped.gene_names <- genomic.features[
      subjectHits(target.intersect.feature), ] %>%
      as.data.frame %>% pull(gene_name)
    
    # Summarize the read counts of matching peaks into read counts
    # per gene for a given motif
    overlapped.target.readcount <- matched.target.GRanges[
      queryHits(target.intersect.feature), ] %>%
      mcols %>%
      as.data.frame %>%
      add_column("target_gene" = overlapped.gene_names) %>% # add gene_name to matched ranges
      dplyr::select(target_gene,A1:A12) %>% # selecct read count columns only
      filter(target_gene %in% genelist) %>% # keep rows with genes in the genelist
      group_by(target_gene) %>%
      summarise_all(sum) # sum read counts per gene, as multiple peaks can fall into one gene


    return(overlapped.target.readcount) # return a dataframe
}



#------------------------------------------------------------------------------
# READ IN INPUTS
#------------------------------------------------------------------------------

# consensus peaks read counts
raw_counts = fread(inputfiles$raw_counts, header = TRUE,
  col.names = c("chrom", "start", "end", paste("A", 1:12, sep = ""))
)

# each peak must > 15 reads across all samples
keep_rows = which(rowSums(raw_counts[, -c("chrom", "start", "end")]) > 15)
raw_counts = raw_counts[keep_rows]

# raw peak counts
raw_peak_counts <- makeGRangesFromDataFrame(raw_counts, keep.extra.columns = T)

# normalize counts sqrt of reads per million
peaks = raw_counts[, .(chrom, start, end)]
norm_counts = raw_counts[, A1:A12] %>%
  apply(2, function(x) sqrt(x * 1e6 / sum(x))) %>%
  as.data.table
norm_peak_counts = cbind(peaks, norm_counts) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

# JASPAR 2020 TF motifs
jaspar_2020 = readJASPARMatrix(inputfiles$jaspar_2020, "PFM")

# CRISPR screen HITS and TF list (in screen library)
hitlist = fread(inputfiles$hitlist)
tflist = fread(inputfiles$tflist)

# make motif lookup table
motif_lookup <- name(jaspar_2020) %>%
                data.frame(stringsAsFactors = F) %>%
                rownames_to_column("motif") %>%
                mutate(gene_name = str_extract(`.`, "[a-zA-Z:0-9\\-]+")) %>%
                select(motif, gene_name)
# add two columns to indicate whether the motif is a TF or is a Hit
motif_lookup = mutate(motif_lookup,
  is_hit = mystrfunc(gene_name, hitlist$gene_name),
  is_tf = mystrfunc(gene_name, tflist$gene_name)
)
# motif names of (almost) all TFs
tf.motif.list <- motif_lookup %>% filter(is_tf == "Yes") %>% pull(motif)

# get annotated promoter region bed
# note the latest protein coding genes are updated with gene names
# and gene_id to be consistent with hgnc

# check first line, specifically column2 of first line
# if it has header, it would be character, otherwise all digits
firstline = readLines(inputfiles$gene_region, n=1) %>%
  str_split("\t", simplify = T) %>%
  .[1,2]

if (str_detect(firstline, "^\\d+$")) {
  gene_region <- fread(inputfiles$gene_region, header = F,
    col.names = c("seqname","start","end","gene_id","gene_name","strand")) %>%
    .[, .(seqname, start, end, strand, gene_id, gene_name)]
} else {
  gene_region <- fread(inputfiles$gene_region, header = T,
    col.names = c("seqname","start","end","gene_id","gene_name","strand")) %>%
    .[, .(seqname, start, end, strand, gene_id, gene_name)]
}

# convert to GRanges object
gene_region <- makeGRangesFromDataFrame(gene_region, keep.extra.columns = T)



#------------------------------------------------------------------------------
# Look for motif matches in peaks
#------------------------------------------------------------------------------

# get motif matches with peaks, (takes about 97 seconds to run)

if (norm_flag) {
  cat("Compute motif footprints on targets using normalized ATAC counts...\n")
  Jaspar_ix <- matchMotifs(jaspar_2020, norm_peak_counts,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    out = "matches", p.cutoff = 5e-5)
} else {
  cat("Compute motif footprints on targets using raw ATAC counts...\n")
  Jaspar_ix <- matchMotifs(jaspar_2020, raw_peak_counts,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    out = "matches", p.cutoff = 5e-5)

}



#------------------------------------------------------------------------------
# WORK WITH TFs - from REGULATOR (motifs) to (find) TARGETS
#------------------------------------------------------------------------------

# takes a couple minutes
cat("from REGULATORS look for their TARGETS...\n")
tictoc::tic()

tf_motifs = filter(motif_lookup, is_tf == "Yes") %>% pull(motif) %>% sort
names(tf_motifs) = tf_motifs


# Given a list of motifs, find their binding target genes, ie regulator -> targets
reg2tar <- imap_dfr(tf_motifs,
  ~ getTargetsFromMotif(.x, Jaspar_ix, gene_region, 100, tflist$gene_name) %>%
  add_column("regulator_motif" = .y, .before = "target_gene")) %>%
  unique

# Get gene name from motif name and summarize read 
# counts per interaction (regulator - target combinations)
reg2tar <- left_join(reg2tar, motif_lookup[, 1:2],
                     by = c("regulator_motif" = "motif")) %>%
  dplyr::rename("regulator_gene" = "gene_name") %>%
  select(regulator_gene, target_gene, A1:A12) %>%
  group_by(regulator_gene, target_gene) %>%
  summarise_all(sum)

# Combine replicate 1 and replicate 2 by taking the mean
# Result in just ES and 5 time points
rep1 <- paste0(rep("A", 6), seq(1, 11, 2))
rep2 <- paste0(rep("A", 6), seq(2, 12, 2))
reg2tar <- map2_df(reg2tar[, rep1], reg2tar[, rep2], 
                   ~ map2_dbl(.x, .y, ~ mean(c(.x, .y))) ) %>%
  add_column("regulator_gene" = reg2tar$regulator_gene,
             "target_gene" = reg2tar$target_gene,
             .before = "A1")
reg2tar <- unique(reg2tar)
names(reg2tar) <- c("regulator_gene", "target_gene", "ES", "H1", "H4", "H16", "D1", "D4")

tictoc::toc()



#------------------------------------------------------------------------------
# WRITE OUTPUT
#------------------------------------------------------------------------------

cat(paste0("Write computed interactions to excel file: ", outputfiles$interactions, "\n"))
WriteXLS::WriteXLS(
  x = c("reg2tar"),
  ExcelFileName = outputfiles$interactions,
  SheetNames = c("Regulator_to_Target")
)




#------------------------------------------------------------------------------
# WORK WITH TFs - from TARGETS to (find) their REGULATORS (DEPRECATED)
#------------------------------------------------------------------------------
# THIS IS NOT NEEDED AS I have confirmed using regulator to find targets
# finds the exact same interaction pairs and reads as from targets to 
# regulators, while saving a lot of time. 

# takes half an hour
# cat("from TARGETS look for their REGULATORS...\n")
# tictoc::tic()

# # given a list of gene, find the regulators of each gene
# tar2reg <- mclapply(tflist$gene_name,
#   function(x) {
#     getRegulatorsOfGene(x,Jaspar_ix, gene_region, 100, tf_motifs)
#     }, mc.cores = NCore)

# # remove null elements
# tar2reg <- map_lgl(tar2reg, ~ ! is_null(.x)) %>% which %>% tar2reg[.]
# # combine into one dataframe
# tar2reg <- do.call(rbind, tar2reg)

# # combine replicate 1 and replicate 2 by taking the mean
# # Result in just ES and 5 time points
# rep1 <- paste0(rep("A", 6), seq(1, 11, 2))
# rep2 <- paste0(rep("A", 6), seq(2, 12, 2))

# tar2reg <- map2_df(tar2reg[, rep1], tar2reg[, rep2],
#   ~ map2_dbl(.x, .y, ~ mean(c(.x, .y)))) %>%
#   add_column("regulator_motif" = tar2reg$regulator_motif,
#              "target_gene" = tar2reg$target_gene, .before = "A1")
# names(tar2reg) <- c("regulator_motif", "target_gene",
#   "ES", "H1", "H4", "H16", "D1", "D4")

# # replace motif names with gene names
# tar2reg <- left_join(tar2reg, motif_lookup[, 1:2],
#                      by = c("regulator_motif" = "motif")) %>%
#   select(gene_name, target_gene:D4) %>%
#   group_by(gene_name, target_gene) %>%
#   summarise_all(sum) %>%
#   dplyr::rename("regulator_gene" = "gene_name") %>%
#   unique %>%
#   ungroup

# tictoc::toc()



# getRegulatorsOfGene <- function(gene_name, motif.matches, genomic.features, min.overlap = 10, motif.list) {
    
#     #' @param gene_name: string, gene name
#     #' @param motif.matches: matchMotif object, e.g. motifMatch() result
#     #' of peaks and jaspar 2020 motifs
#     #' @param countGRanges: GenomicRanges object, e.g. a normalized read matrix,
#     #'        converted into GRanges object
#     #' @param genomic.features: GenomicRanges object, e.g. a bed formated
#     #' promoter region coordinates, converted into GRanges object
#     #' @param min.overlap: integer, minimum overlap in base pairs
#     #' @param motif.list: vector, a list of motif names. Only return results if
#     #'        found motifs are part of this list

#     #' @returns regulators.readcount: list, each element's name is a
#     #' regulator's motif name, the values are summarised readcounts of
#     #' matching/overlapping peaks for this regulator
    
#     match.matrix <- motifMatches(motif.matches) # motif match logic matrix

#     # peak coordiantes, including read counts in mcols
#     match.matrix.GRanges <- rowRanges(motif.matches)

#     # target gene's annotation coordinates
#     gene.GRanges <- genomic.features[genomic.features$gene_name == gene_name]
    
#     # intersect peak coordinates with target gene's annotation coordinates
#     # to get peaks that belong to target gene
#     match.intersect.gene <- GenomicRanges::findOverlaps(
#       query = match.matrix.GRanges,
#       subject = gene.GRanges,
#       minoverlap = min.overlap,
#       ignore.strand = T
#       )
    
#     # row index of peaks that fall in target gene coordinates
#     overlapped.match.index <- queryHits(match.intersect.gene)
    
#     # Once extracted peaks that fall into target gene region,
#     # get the motifs that bind to these peaks.
#     # Note slight difference in extracting motif names when there
#     # are only 1 peak versus more than 1 peaks
#     if (length(overlapped.match.index) == 1) { # if only one peak in target gene
#         gene.regulators <- match.matrix[overlapped.match.index, ] %>%
#           .[.] %>% names
#     } else { # if multiple peaks in target gene
#         gene.regulators <- match.matrix[overlapped.match.index, ] %>%
#           colSums %>% .[.>0] %>% names
#     }
    
#     # Only keep motifs that are part of a given list
#     gene.regulators <- gene.regulators[gene.regulators %in% motif.list] 
    
#     # Get the read counts of the binding sites for each motif (regulator)
#     regulators.readcount <- list()
#     for (regulator in gene.regulators) {
#         readcount <- match.matrix.GRanges[overlapped.match.index, ] %>% # readcount of each peak in target gene region
#           as.data.frame %>% select(A1:A12) %>% # extract read counts columns
#           `*`(match.matrix[overlapped.match.index, regulator]) %>% # multiply 0 or 1 based on if the peak has a match to this regulator or not
#           colSums() # sum up all the reads from peaks that match to the regulator binding site
        
#         regulators.readcount[[regulator]] <- readcount
#     }
    
#     # Return a dataframe, each row gives the name of the motif that 
#     # binds to the target gene, along with observed read count assciated
#     # with the matching peaks. Some genes may have no matches
#     if (length(regulators.readcount) > 0) {
#         regulators.readcount <- do.call(rbind, regulators.readcount) %>%
#           as.data.frame %>%
#           add_column("target_gene" = gene_name, .before = "A1") %>%
#           rownames_to_column("regulator_motif")
        
#         return(regulators.readcount)
#     }
# }


