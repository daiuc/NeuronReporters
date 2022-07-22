# combine RSEM counts from each sample into a single count table file
# this script produces a combined count table at
# "Results/RNAseq/RSEM/RSEM_combined.genes.results"


library(tidyverse)
library(tximport)

#Log = file(snakemake@log[[1]], open="wt")
#sink(log, type = c("output", "message"))

output.file = snakemake@output[[1]]

RUN_MODE = 1
if (RUN_MODE == 1) {
    # snakemake mode
    input.files = flatten_chr(snakemake@input)
} else {
    # Rscript mode
    input.files = list.files("Results/RNAseq/RSEM", "*genes.results", full.names = T) %>% sort
}

print(input.files)

getColumnMeta = function(input.files) {
    df = str_extract_all(input.files, "[a-z0-9]+_[a-z0-9]{2}_[0-9]{1}_[A-Z]{1}") %>% flatten_chr %>% str_split("_", simplify = 
    T)  %>% as.data.frame
    colnames(df) = c("treatment", "timepoint", "rep", "batch")
    return(df)
}

count.table = tximport(input.files, type = "rsem", txIn = FALSE, txOut = FALSE)
colMeta = getColumnMeta(input.files)

cts.cols = str_extract_all(input.files, "[a-z0-9]+_[a-z0-9]{2}_[0-9]{1}_[A-Z]{1}") %>% flatten_chr
cts = count.table$counts %>% 
    as.data.frame() %>% 
    mutate(across(.cols = everything(), round))
colnames(cts) = cts.cols
cts = rownames_to_column(cts, "gene_id")


# first remove gene_ids that has "_PAR_Y" in the end
# then remove any gene_id digits after "."
cts = filter(cts, !str_detect(gene_id, "PAR_Y")) %>% 
        mutate_at("gene_id", ~ str_split(.x, "\\.", simplify = T) %>% .[,1])

write_tsv(cts, output.file)
