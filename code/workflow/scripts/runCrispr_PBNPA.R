library(tidyverse)
library(data.table)
library(PBNPA)
library(furrr)


rc = fread("resources/crispr/counts/wrangled_raw_counts.tsv")

colnames(rc) = str_replace_all(colnames(rc),
                               pattern = "(sgRNA)|(Gene)",
                               replacement = c("sgRNA" = "sgRNA_name", "Gene" = "Gene_name"))

# get sgRNA index
sgRNA.index = rownames(rc) %>% as.numeric

# get Gene index
unique.gene.list = unique(rc$Gene_name)
Gene.index = map_int(rc$Gene_name, ~ which(unique.gene.list == .x))

rc[, `:=`(
  sgRNA = sgRNA.index,
  Gene = Gene.index
)]

Tdneg.vs.iN = rc[, .(sgRNA, Gene, initial.count = Day14_iN, final.count = Day14_Tdneg)] %>% list
Tdpos.vs.iN = rc[, .(sgRNA, Gene, initial.count = Day14_iN, final.count = Day14_Tdpos)] %>% list
Tdneg.vs.Tdpos = rc[, .(sgRNA, Gene, initial.count = Day14_Tdpos, final.count = Day14_Tdneg)] %>% list

results = list(Tdneg.vs.iN = Tdneg.vs.iN,
               Tdpos.vs.iN = Tdpos.vs.iN,
               Tdneg.vs.Tdpos = Tdneg.vs.Tdpos)

plan(multisession, workers = min(4, availableCores()))

results = future_map(results, ~ PBNPA(.x, sim.no = 100, fdr = .1))

intersect(unique.gene.list[results$Tdneg.vs.iN$neg.gene],
unique.gene.list[results$Tdneg.vs.Tdpos$neg.gene]) %>% intersect(hitlist$gene_name)

hitlist = fread("resources/crispr/Hitlist_20191230.csv")
