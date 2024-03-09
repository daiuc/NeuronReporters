### Scan reference genome hg38 fasta, with motif PWM and produce TFBS coordinates

library(tidyverse)
library(motifmatchr)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)
register(MulticoreParam(4))


# 1. read in JASPAR 2020 motif PWMs and protein coding gene coordinates

JAS2020 = readJASPARMatrix("/gpfs/commons/groups/sanjana_lab/cdai/TFscreen/atac/JASPAR2020_combined_matrices_20191030.txt", 
                                matrixClass = "PFM")

GENES = data.table::fread("/gpfs/commons/groups/sanjana_lab/cdai/TFscreen/atac/Protein_coding_genes_Up_2k_20191230.bed",
                          col.names = c("chr","start","end","gene_id","gene_nm","strand"))

# convert to genomic ranges object
GENES = GenomicRanges::makeGRangesFromDataFrame(GENES,keep.extra.columns = T)


# 2. for each motif, look for all the TFBS in gene body

# get motif matches with peaks, (takes about 97 seconds to run)

#TFBS = matchMotifs(JAS2020[[1]], GENES[1:10], genome = BSgenome.Hsapiens.UCSC.hg38, out = "positions")


# convert granges to bed
GRangestoBed = function(GRangeObj) {
    return(data.frame(chr = seqnames(GRangeObj),
                      start = start(GRangeObj),
                      end = end(GRangeObj),
                      name = rep(".", length(GRangeObj)),
                      scores=rep(".", length(GRangeObj)),
                      strand = strand(GRangeObj)
                      ))
}


#x = matchMotifs(JAS2020[1:2], GENES[1:10], genome = BSgenome.Hsapiens.UCSC.hg38, out = "positions")

#Jaspar_ix <- matchMotifs(jaspar_2020, norm_peak_counts, 
#                         genome = BSgenome.Hsapiens.UCSC.hg38,
#                         out = "matches", p.cutoff = 5e-5)


# 3. for each gene, look for all the motifs that has a TFBS within its gene body region


GENES.sample = GENES[1:10]

TFBS_PER_GENE = lapply(seq_along(GENES), 
                       function(i) matchMotifs(JAS2020, GENES[i], 
                                               genome=BSgenome.Hsapiens.UCSC.hg38, 
                                               out = "positions"))

names(TFBS_PER_GENE) = mcols(GENES.sample)$gene_id

# remove results that have 0 genomic ranges
TFBS_PER_GENE = map(TFBS_PER_GENE, ~ .x[unlist(lapply(.x, function(gr) length(gr) > 0))])
TFBS_PER_GENE = TFBS_PER_GENE[unlist(lapply(TFBS_PER_GENE, function(l) length(l) > 0))]
lapply(TFBS_PER_GENE, function(ls) lapply(ls, function(l) GRangestoBed(l)))

TFBS_ALL_GENE = list()
for (nm in names(TFBS_PER_GENE)) {
    ls = TFBS_PER_GENE[[nm]]
    ls = map2_df(as.list(ls), names(ls), ~ GRangestoBed(.x) %>% add_column("gene_id" = nm, "motif" = .y) )
    TFBS_ALL_GENE = c(TFBS_ALL_GENE, list(ls))
    
}
TFBS_ALL_GENE = do.call(rbind, TFBS_ALL_GENE)

GENE_LOOKUP = as.data.frame(mcols(GENES))

TFBS_ALL_GENE = left_join(TFBS_ALL_GENE, GENE_LOOKUP, by = "gene_id")

data.table::fwrite(TFBS_ALL_GENE, "Results/Analyses/JASPAR2020_TFBS_on_hg38_geneBody_up2k.txt")
