#' @title         Infer interactions
#' @author        Chao Dai
#' @description   Infer interactions of regulator-target TFs based on scanning
#'                of DNA binding motifs on ATAC-seq peak region



#--------------------------------------------------------------------------
# Set up ------------------------------------------------------------------
#--------------------------------------------------------------------------

library(data.table)
library(tidyverse)
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(Matrix)
library(BSgenome.Hsapiens.UCSC.hg38)
library(naturalsort)
library(parallel)
ncores = min(10, detectCores())

set.seed(2019)

#--------------------------------------------------------------------------
# Snakemake objects -------------------------------------------------------
#--------------------------------------------------------------------------

SnakeMode = T

if (SnakeMode) {
  print("### Running in snakemake script mode")
  peak.file = snakemake@input[["peaks"]]
  counts.file = snakemake@input[["atac_counts"]]
  jaspar.file = snakemake@input[["jaspar2020"]]
  tflist.file = snakemake@input[["tflist"]]
  hitlist.file = snakemake@input[["hitlist"]]
  genome.features = snakemake@input[["features"]]
  dge.file = snakemake@input[["deseq2"]]
  rna_sample_anno = snakemake@input[["rna_sample_names"]]
  #out.file1 = snakemake@output[[1]] # nested
  out.file = snakemake@output[[1]] # expanded
} else {
  print("### running in interactive mode")
  peak.file = "/gpfs/commons/groups/sanjana_lab/cdai/TFscreen/atac/diffbind/diffbind_consensu_min2overlap_NO_RECENTER.bed"
  counts.file = "/gpfs/commons/groups/sanjana_lab/cdai/TFscreen/atac/diffbind/diffbind_consensu_min2overlap_NO_RECENTER.txt"
  jaspar.file = "/gpfs/commons/groups/sanjana_lab/cdai/TFscreen/atac/JASPAR2020_combined_matrices_20191030.txt"
  tflist.file = "/gpfs/commons/groups/sanjana_lab/cdai/TFscreen/TFlist_20191230.csv"
  hitlist.file = "/gpfs/commons/groups/sanjana_lab/cdai/TFscreen/Hitlist_20191230.csv"
  genome.features = "resources/annotations/hs38/gencode_v31_intergeneic_TSSup100kb-TSSup2kb.bed"
  dge.file = "/gpfs/commons/groups/sanjana_lab/cdai/TFscreen/RNA-seq_timepoint_deseq_result_20200102.xlsx"
  rna_sample_anno = "/gpfs/commons/groups/sanjana_lab/cdai/TFscreen/RNASeqSampleNames.csv"
}

#--------------------------------------------------------------------------
# Support functions -------------------------------------------------------
#--------------------------------------------------------------------------

# custom function to check if a motifname contains a hit gene name or TF gene name
# the reason is motifname sometimes contains "(var.2)" or other characters
mystrfunc = function(motifname, target_match) {
  #' @param motifname: scalar string, e.g. TFAP4(var.2)
  #' @param target_match: string vector, e.g. hitlist, or tflist

  x = str_match_all(motifname, "[A-Z0-9\\-]{2,13}") %>% map(~.x[, 1]) %>% .[[1]]

  in_list = x %in% target_match %>% purrr::reduce(., `|`)

  return(in_list)
}


getTargetsFromMotif = function(motif_name, matchMotifs_obj, genomic.features, min.overlap = 10, genelist) {
  #' @desscription            Given a motif name (JASPAR), find its binding targets with read counts, observed in ATAC-seq
  #' @param motif_name:       string, this is the jaspar motif name, not gene name
  #' @param motif.matches:    matchMotif object, e.g. motifMatch() result of peaks and jaspar 2020 motifs
  #' @param genomic.features: GenomicRanges object, e.g. a bed formated promoter region coordinates, converted
  #'                          into GRanges object
  #' @param min.overlap:      integer, minimum overlap in base pairs
  #' @param genelist:         vector, a list of gene_names. Note these gene names are some times not the same
  #'                          as motif names, hence the need of motif_lookup
  #' @return                  dataframe, containing read counts of genes that are targets of a given regulator.
  #'                          Read counts sum of all peaks that match to the regulator's motif.


  #' Get motif matches matrix [i x j] (rows are peak coordinates, columns are each motif name,
  #' values are logical values indicating match)

  # sparse logic matrix where rows are peak coordinates (GRange) and columns are Jaspar motif names
  # 1 = motif j has a match on peak i.
  match.matrix <- motifMatches(matchMotifs_obj)

  # get Grange coordinates of peaks that has a match to motif j
  GRanges.index = which(match.matrix[, motif_name]) # which row (peak) has a match to motif j
  matched.target.GRanges <- rowRanges(matchMotifs_obj)[GRanges.index] # further convert row index to GRange coordinates

  # Find overlaps between:
  # 1. peaks that match binding sites of a given motif, and
  # 2. coordinates of gene annotation
  target.intersect.feature <- GenomicRanges::findOverlaps(
      query = matched.target.GRanges,
      subject = genomic.features,
      minoverlap = min.overlap,
      ignore.strand = T
    )

  # Get the genes (from annoation) that overlap with peaks (that have binding sites of given motif)
  overlapped.gene_names <- genomic.features[subjectHits(target.intersect.feature), ] %>%
    as.data.frame %>%
    pull(gene_name)

  # Summarize the read counts of matching peaks into read counts per gene for a given motif
  overlapped.target.readcount <- matched.target.GRanges[queryHits(target.intersect.feature), ] %>%
    mcols %>% as.data.frame %>%
    add_column("target_gene" = overlapped.gene_names) %>% # add gene_name to matched ranges
    dplyr::select(target_gene,A1:A12) %>% # selecct read count columns only
    filter(target_gene %in% genelist) %>% # keep rows with genes in the genelist
    group_by(target_gene) %>%
    summarise_all(sum) # sum read counts per gene, as multiple peaks can fall into one gene

  return(overlapped.target.readcount) # return a data.frame
}


my_unnest = function(listcolumn) {
  ncols = length(listcolumn[[1]])
  if (ncols == 4) {
    unlist(listcolumn) %>%
      matrix(ncol = ncols, byrow = T, dimnames = list(NULL, c("ES", "H16", "D1", "D4"))) %>%
      as.data.table
  } else if (ncols == 3) {
    unlist(listcolumn) %>%
      matrix(ncol = ncols, byrow = T, dimnames = list(NULL, c("H16", "D1", "D4"))) %>%
      as.data.table
  }
}


getRNA_counts = function(gene, count_dt, col_name) {
  #' @param gene this could be dimers as in JASPAR, but no ver numbers
  #' @param count_dt dt of rna-seq
  #' @param col_name column name in count_dt that store count data per gene

  genes = str_split(gene, "::", simplify = T) %>% as.vector
  mycount = count_dt[gene_name %in% genes, col_name, with = F] %>%
    as.list %>% flatten
  if (length(mycount) > 0) {
    mycount = my_unnest(mycount) %>%
      colMeans
  } else {
    mycount = c(0,0,0,0)
  }
  return(mycount)
}

getRNA_l2fc_fdr = function(gene, dt) {
  #' @param gene this could be dimers as in JASPAR, but no ver numbers
  #' @param dt dt that contains l2fc and fdr columns

  genes = str_split(gene, "::", simplify = T) %>% as.vector
  dt = dt[gene_name %in% genes]
  if (nrow(dt) > 0) {
    fdr.mins = my_unnest(dt$fdr) %>% as.matrix %>% rowMin
    row.id = which(fdr.mins == min(fdr.mins))

    l = list(l2fc = dt$l2fc[[row.id]],
             fdr = dt$fdr[[row.id]])
  } else {
    l = list(l2fc = c(0,0,0),
             fdr = c(1,1,1))
  }
  return(l)
}

getHitStatus = function(gene, hits) {
  #' @param gene  this could be dimers but no ver numbers
  #' @param hits  hitlist, a vector
  #' @return      a vector of T/F

  genes = str_split(gene, "::", simplify = T) %>% as.vector
  if_else(genes %in% hits, T, F) %>% purrr::reduce(., `|`)
}



#--------------------------------------------------------------------------
# Load in data ------------------------------------------------------------
#--------------------------------------------------------------------------
print("### Load in data")

# genome features
genome.features = fread(genome.features, header = T) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T, ignore.strand = F,
                           seqnames.field = "#chrom")

# consensus peaks
peaks = chromVAR::getPeaks(peak.file, sort_peaks = F)

# Consensus peaks read counts
raw_counts <- fread(counts.file, header = T)

# each peak must> 15 reads across all samples
keep_rows <- which(rowSums(select(raw_counts, A1:A12)) > 15 )
raw_counts <- raw_counts[keep_rows, ]

# peak coordinates
peaks <- raw_counts[, .(Chr, Start, End)]

# normalized counts: squareroot reads per million
norm_counts <- select(raw_counts, A1:A12) %>%
  as.matrix %>%
  apply(., 2, function(x) sqrt(x * 1000000 / sum(x)) ) %>%
  as.data.table
norm_peak_counts <- cbind(peaks, norm_counts) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = T)
raw_peak_counts = makeGRangesFromDataFrame(raw_counts, keep.extra.columns = T)

# JASPAR Annotations, downloaded JASPAR2020
jaspar_2020 <- TFBSTools::readJASPARMatrix(jaspar.file, matrixClass = "PFM")


# read in hit list and TF list, note the gene names are matched and transformed to be consistent with approved symbols according to HGNC genenames.org
# hit list from crispr screen
hitlist <- read.csv(hitlist.file, stringsAsFactors=F) %>% pull(hgnc_symbol)
# genomewide TF list
tflist <- read.csv(tflist.file, stringsAsFactors=F) %>% pull(hgnc_symbol)


# construct base motif_lookup table
motif_lookup <- TFBSTools::name(jaspar_2020) %>%
  data.table(motif = .)
motif_lookup[, gene_name := str_extract(motif, "[a-zA-Z:0-9\\-]+")]

# add two columns to indicate whether the motif is a TF or is a Hit
motif_lookup[, `:=`(
  is_hit = map_chr(gene_name, ~ if_else(mystrfunc(.x, hitlist), "Yes", "No")),
  is_tf = map_chr(gene_name, ~ if_else(mystrfunc(.x, tflist), "Yes", "No"))
)]



# motif names of (almost) all TFs
tf.motif.list = motif_lookup[is_tf == "Yes", motif]

# get motif matches with peaks, (takes about 97 seconds to run)
Jaspar_ix <- matchMotifs(jaspar_2020, norm_peak_counts,
                         genome = BSgenome.Hsapiens.UCSC.hg38,
                         out = "matches", p.cutoff = 5e-5)

# same but with raw readcounts
Jaspar_ix_raw <- matchMotifs(jaspar_2020, raw_peak_counts,
                         genome = BSgenome.Hsapiens.UCSC.hg38,
                         out = "matches", p.cutoff = 5e-5)


#--------------------------------------------------------------------------
# Given motifs of regulators, find their targets --------------------------
#--------------------------------------------------------------------------
print("### Find motif targets")
# get list of regulator TF motif names
tf_motifs = motif_lookup[is_tf == "Yes", motif] %>% sort
names(tf_motifs) = tf_motifs

# Search targets of each motif --------------------------------------------
motif_to_targets = mclapply(
  tf_motifs,
  function(x) {getTargetsFromMotif(x, Jaspar_ix, genome.features, 100, tflist) %>%
    add_column("regulator_motif" = x, .before = "target_gene") %>%
    as.data.table},
    mc.cores = ncores
)
motif_to_targets = rbindlist(motif_to_targets)

# search targets of each motif - raw
motif_to_targets_raw = mclapply(
  tf_motifs,
  function(x) {getTargetsFromMotif(x, Jaspar_ix_raw, genome.features, 100, tflist) %>%
    add_column("regulator_motif" = x, .before = "target_gene") %>%
    as.data.table},
    mc.cores = ncores
)

motif_to_targets_raw = rbindlist(motif_to_targets_raw)

# Combine replicate 1 and replicate 2 by taking the mean
# Result in just ES and 5 time points
rep1 <- paste0(rep("A", 6), seq(1, 11, 2))
rep2 <- paste0(rep("A", 6), seq(2, 12, 2))



#  all putative interactions between regulator and target TFs -------------
print("### Construt interaction matrix with motifmatchr and results")

# takes about 4 minutes
interactions = motif_to_targets[,
  map2(.SD[, rep1, with = F],
       .SD[, rep2, with = F],
       ~ map2_dbl(.x, .y, function(a,b) mean(c(a,b)))),
  by = c("regulator_motif", "target_gene")
]

# rename columns to timepoints
interactions = interactions[, .(
  regulator_gene = regulator_motif,
  target_gene,
  ES = A1,
  H1 = A3,
  H4 = A5,
  H16 = A7,
  D1 = A9,
  D4 = A11
)]

# fix regulator_gene names some times has version number
interactions = interactions[
  motif_lookup[, 1:2], on = c(regulator_gene = "motif"), nomatch = NULL
  ][, .(regulator_gene = gene_name, target_gene, ES, H1, H4, H16, D1, D4)][
    , .(ES = mean(ES), H1 = mean(H1), H4 = mean(H4),
        H16 = mean(H16), D1 = mean(D1), D4 = mean(D4)),
    by = .(regulator_gene, target_gene)
  ]



# all putative interactions - raw counts ----------------------------------
print("### integrate motifmatchr results with raw ATAC-seq counts")

interactions_raw = motif_to_targets_raw[,
                                map2(.SD[, rep1, with = F],
                                     .SD[, rep2, with = F],
                                     ~ map2_dbl(.x, .y, function(a,b) mean(c(a,b)))),
                                by = c("regulator_motif", "target_gene")
]

interactions_raw = interactions_raw[, .(
  regulator_gene = regulator_motif,
  target_gene,
  ES = A1,
  H1 = A3,
  H4 = A5,
  H16 = A7,
  D1 = A9,
  D4 = A11
)]

# fix regulator_gene names some times has version number
interactions_raw = interactions_raw[
  motif_lookup[, 1:2], on = c(regulator_gene = "motif"), nomatch = NULL
  ][, .(regulator_gene = gene_name, target_gene, ES, H1, H4, H16, D1, D4)][
    , .(ES = mean(ES), H1 = mean(H1), H4 = mean(H4),
        H16 = mean(H16), D1 = mean(D1), D4 = mean(D4)),
    by = .(regulator_gene, target_gene)
  ]



# nest normalized motif reads on target gene ------------------------------

interactions = interactions[, -c("H1", "H4")] %>%
  .[, .(atac = pmap(list(ES, H16, D1, D4), function(a,b,c,d) c(ES = a,H16 = b, D1 = c, D4 = d))),
    by = c("regulator_gene", "target_gene")]

interactions_raw = interactions_raw[, -c("H1", "H4")] %>%
  .[, .(atac = pmap(list(ES, H16, D1, D4),
                    function(a,b,c,d) {c(ES = round(a),
                                        H16 = round(b),
                                        D1 = round(c),
                                        D4 = round(d))
                    })), by = c("regulator_gene", "target_gene")]


# bring in raw atac counts
interactions = interactions[, .(regulator_gene, target_gene, atac.norm = atac)
  ][interactions_raw[, .(regulator_gene, target_gene, atac.cnt = atac)],
    on = c("regulator_gene", "target_gene"),
    nomatch = NULL]



#--------------------------------------------------------------------------
# Integrate target gene expression, compute correlation -------------------
#--------------------------------------------------------------------------
print("### integrate gene expression data")

# DeSeq2 results for 3 time points
sheets = c("H15_vs_ES", "D1_vs_ES", "D4_vs_ES")
dge = map(sheets, ~ readxl::read_excel(dge.file, .x) %>%
            select(gene_id, gene_name, log2FoldChange, pvalue, padj)
          )
names(dge) <- factor(c("H16", "D1", "D4"), levels = c("H16", "D1", "D4"))
dge <- map(dge, ~ drop_na(.x, gene_name)) # remove rows with no gene name

# raw counts
dge.rawCounts <- readxl::read_excel(dge.file, "rawCounts")
gene_id_names = select(dge.rawCounts, gene_id, gene_name)

rna_sample_anno <- fread(rna_sample_anno)
rna_sample <- c("S01_B1", "S01_B2", "S02_B1", "S02_B2", "S15_B1", "S15_B2", "S16_B1", "S16_B2",
                "S17_B1", "S17_B2", "S18_B1", "S18_B2", "S07_B1", "S07_B2", "S08_B1", "S08_B2")


# Normalize read counts: first calc reads per million then take sqrt
dge.normCounts <- select(dge.rawCounts, - gene_id, - gene_name) %>%
  apply(2, function(x) sqrt(x * 1e6 / sum(x)) ) %>% as.data.frame

# replicate names
rna.rep1 <- c("S01_B1", "S15_B1", "S17_B1", "S07_B1")
rna.rep2 <- c("S01_B2", "S15_B2", "S17_B2", "S07_B2")
rna.rep3 <- c("S02_B1", "S16_B1", "S18_B1", "S08_B1")
rna.rep4 <- c("S02_B2", "S16_B2", "S18_B2", "S08_B2")




# normalized rna-seq counts -----------------------------------------------


# combine 4 replicates into 1, using average
dge.normCounts = pmap(
    list(dge.normCounts[, rna.rep1], dge.normCounts[, rna.rep2],
         dge.normCounts[, rna.rep3], dge.normCounts[, rna.rep4]),
    ~ pmap_dbl(list(..1, ..2, ..3, ..4), function(a,b,c,d) mean(c(a,b,c,d)))) %>%
  as.data.table %>%
  .[, c(list(gene_id = gene_id_names$gene_id,
             gene_name = gene_id_names$gene_name),
        .SD)]

dge.normCounts = dge.normCounts[
  , .(gene_id, gene_name,
      ES = S01_B1,
      H16 = S15_B1,
      D1 = S17_B1,
      D4 = S07_B1
  )]


# some gene_ids correspond to multiple or NA gene names,
# remove gene_names = NA, and combine (sum) rows with the same gene_name
dge.normCounts = dge.normCounts[, -c("gene_id"),][, lapply(.SD, sum), by = gene_name][
  !is.na(gene_name)]

# Create nested data.table with each time point expression values in a list
dge.normCounts = dge.normCounts[
  , .(rna.norm = pmap(list(ES, H16, D1, D4),
                      function(a,b,c,d) c(ES = a,H16 = b, D1 = c, D4 = d))
      ), by = gene_name]



# raw rna-seq counts ------------------------------------------------------

dge.rawCounts = pmap(
  list(dge.rawCounts[, rna.rep1], dge.rawCounts[, rna.rep2],
       dge.rawCounts[, rna.rep3], dge.rawCounts[, rna.rep4]),
  ~ pmap_dbl(list(..1, ..2, ..3, ..4),
             function(a,b,c,d) mean(c(a,b,c,d)))) %>%
  as.data.table %>%
  .[, c(list(gene_id = gene_id_names$gene_id,
             gene_name = gene_id_names$gene_name),
        .SD)]

dge.rawCounts = dge.rawCounts[
  , .(gene_id, gene_name,
      ES = round(S01_B1),
      H16 = round(S15_B1),
      D1 = round(S17_B1),
      D4 = round(S07_B1)
  )]

dge.rawCounts = dge.rawCounts[, -c("gene_id"),][, lapply(.SD, sum), by = gene_name][
  !is.na(gene_name)]

dge.rawCounts = dge.rawCounts[
  , .(rna.cnt = pmap(list(ES, H16, D1, D4),
                      function(a,b,c,d) c(ES = a,H16 = b, D1 = c, D4 = d))
  ), by = gene_name]



# integrate target gene expression (normalized) ---------------------------

interactions = interactions[dge.normCounts[, .(gene_name, rna.tar.norm = rna.norm)],
                            on = c(target_gene = "gene_name"), nomatch = NULL]


# integrate target gene expression (raw) ----------------------------------

interactions = interactions[dge.rawCounts[, .(gene_name, rna.tar.cnt = rna.cnt)],
             on = c(target_gene = "gene_name"), nomatch = NULL]


# integrate regulator gene expression (raw) -------------------------------
# some regulators do not have a match in gene expression raw counts
interactions = interactions[,
  .(regulator_gene, target_gene, atac.norm, atac.cnt, rna.tar.norm, rna.tar.cnt,
    rna.reg.cnt = mclapply(
      regulator_gene,
      function(x) getRNA_counts(x, dge.rawCounts, 'rna.cnt') %>% round,
      mc.cores = ncores
      ))]


# integrate regulator gene expression (normalized) ------------------------
interactions = interactions[,
                            .(regulator_gene, target_gene, atac.norm, atac.cnt,
                              rna.tar.norm, rna.tar.cnt, rna.reg.cnt,
                              rna.reg.norm = mclapply(
                                regulator_gene,
                                function(x) getRNA_counts(x, dge.normCounts, 'rna.norm'),
                                mc.cores = ncores
                              ))]



# Integrate DESeq2  log2FC, FDR at timepoints against ES ------------------
print("### Integrate DESeq2 results: log2fc, FDR")

# first remove any genes with NA gene ID, or 1 gene_name to many gene_ids
# in case of multi gene_id, select the one with the smallest pvalue
dge = mclapply(dge,
  function(x) {
    filter(x, !is.na(gene_name)) %>%
    select(-gene_id) %>%
    group_by(gene_name) %>%
    mutate(rk = rank(pvalue, ties.method = "first")) %>%
    filter(rk <= 1) %>%
    select(-rk)},
  mc.cores = ncores
)

# extract l2fc and fdr (padj)
l2fc = data.table(gene_name = dge$H16$gene_name,
                  H16 = dge$H16$log2FoldChange,
                  D1 = dge$D1$log2FoldChange,
                  D4 = dge$D4$log2FoldChange) %>%
  .[,
    .(l2fc = pmap(list(H16, D1, D4),
                  function(a,b,c) c(H16 = a, D1 = b,D4 = c))
      ), by = gene_name]

impute_fdr = function(fdr) {
  if_else(is.na(fdr), 1, fdr)
}

fdr = data.table(gene_name = dge$H16$gene_name,
                 H16 = impute_fdr(dge$H16$padj),
                 D1 = impute_fdr(dge$D1$padj),
                 D4 = impute_fdr(dge$D4$padj)) %>%
  .[,
    .(fdr = pmap(list(H16, D1, D4),
                 function(a,b,c) c(H16 = a, D1 = b, D4 = c))
      ), by = gene_name]

fc_fdr <- cbind(l2fc, fdr[,2])

# get l2fc & fdr for regulators
regulators = interactions$regulator_gene %>% unique
names(regulators) = regulators
regulators.fc_fdr = mclapply(regulators,
                             function(x) getRNA_l2fc_fdr(x, fc_fdr),
                             mc.cores = ncores)
# convert to dt
regulators.fc_fdr = data.table(regulator_gene = regulators,
           l2fc = map(regulators.fc_fdr, ~.x$l2fc),
           fdr = map(regulators.fc_fdr, ~.x$fdr))

# integrate regulator l2fc and fdr
interactions = interactions[
  regulators.fc_fdr[, .(regulator_gene, rna.reg.l2fc = l2fc, rna.reg.fdr = fdr)],
  on = c("regulator_gene"), nomatch = NULL]

# integrate target l2fc and fdr
interactions = interactions[
  fc_fdr[, .(gene_name, rna.tar.l2fc = l2fc, rna.tar.fdr = fdr)],
  on = c(target_gene = "gene_name"), nomatch = NULL
]


# Run correlation test on ATAC-seq read counts and RNA-seq read co --------
print("### Run correlation test")

# about a minute
cor_test = map2(interactions$atac.norm, interactions$rna.tar.norm,
                ~ cor.test(.x, .y, method = "p"))

# Extract pearson correlation and p-values and store along with read counts for each interactions
interactions[, `:=`(
  corr = map_dbl(cor_test, ~.x$estimate),
  cor_p = map_dbl(cor_test, ~.x$p.value)
)]

# Add in isHit info
print("### Label CRISPR screen hit status")
interactions[, `:=`(
  reg.isHit = mclapply(regulator_gene,
                       function(x) getHitStatus(x, hitlist),
                       mc.cores = ncores),
  tar.isHit = (target_gene %in% hitlist)
  )]

# set order
setorder(interactions, regulator_gene, target_gene)

print("### final interaction matrix")

# column explanations -----------------------------------------------------

column_explanation <- tribble(~name, ~notes,
                              "regulator_gene", "regulator gene name",
                              "target_gene","target gene name",
                              "atac.norm", "ATAC-seq observed read counts of regulator-target interaction, normalized sqrt(reads per million)",
                              "atac.cnt",  "ATAC-seq observed raw read counts of regulator-target interaction",
                              "rna.tar.norm", "gene expression read counts of target gene (normalized to sqrt(reads per million))",
                              "rna.tar.cnt", "gene expression raw read counts of target gene",
                              "rna.reg.cnt", "gene expression raw read counts of regulator gene",
                              "rna.reg.norm", "gene expression read counts of regulator gene(normalized to sqrt(reads per million)",
                              "rna.reg.l2fc", "log2 fold change of regulator gene expression, all against ES",
                              "rna.reg.fdr", "fdr or adjusted p value of regulator gene being differentially expressed, all against ES",
                              "rna.tar.l2fc", "log2 fold change of target gene expression, all against ES",
                              "rna.tar.fdr", "fdr or adjusted p value of regulator gene being differentially expressed, all against ES",
                              "corr", "Pearson correlation between atac.norm and rna.tar.norm",
                              "cor_p", "p value of correlation",
                              "reg.isHit; tar.isHit", "TRUE/FALSE indicates whether regulator or target is a hit in CRISPR screen",
                              "ES, H16, D1, D4", "conditions: ES, 16 Hour, 1 Day, and 4 Day")


print("### column explanations:    ")
print("----------------------------")
print(column_explanation)
print("----------------------------")

# remove large objects
rm(interactions_raw, jaspar_2020, Jaspar_ix, Jaspar_ix_raw,
  l2fc, motif_to_targets, motif_to_targets_raw, norm_counts, 
  norm_peak_counts, peaks, raw_counts, raw_peak_counts, 
  cor_test, dge, dge.normCounts, dge.rawCounts, fc_fdr, 
  fdr, genome.features)

gc()

# output file -------------------------------------------------------------
print("### write out files")

listcols = c("atac.norm", "atac.cnt", "rna.tar.norm",
             "rna.tar.cnt", "rna.reg.cnt", "rna.reg.norm",
             "rna.reg.l2fc", "rna.reg.fdr", "rna.tar.l2fc", "rna.tar.fdr")
# unnest to write out expanded file
interactions_out = interactions[,
  c(list(regulator_gene = regulator_gene, target_gene = target_gene),
    lapply(.SD, function(x) my_unnest(x)),
    list(corr = corr, cor_p = cor_p, reg.isHit = reg.isHit, tar.isHit = tar.isHit)),
  .SDcols = listcols]

# fwrite(interactions, out.file1, sep = "\t") # nested
fwrite(interactions_out, out.file, sep = "\t") # expanded

print("### ------ ALL DONE!")




