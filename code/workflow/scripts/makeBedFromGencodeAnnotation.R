library(tidyverse)
library(data.table)
library(yaml)

Snake_Mode = TRUE

if (Snake_Mode) {
    print("### Running in snakemake pipeline")
    gencode = snakemake@input[[1]]
    tflist = snakemake@input[[2]]
    protein_coding_bed = snakemake@output[[1]]
    tf_bed = snakemake@output[[2]]

} else {
    print("### Running in interactive mode")
    config_data = yaml.load_file("config/config.yaml")
    gencode = config_data$GENCODE_HS38
    tflist = config_data$TF_LIST
}

print("## Read in gencode V31 annotations.")
gencode = fread(gencode)
tflist = fread(tflist)


# remove scaffolds
CHROMS = c(paste("chr", 1:22, sep = ""), "chrX", "chrY")
print("## Extract protein coding genes")
# protein coding genes only at gene level
ptgenes = gencode[
    gene_type %in% c("protein_coding") & 
    feature %in% c("gene") &
    seqname %in% CHROMS,
    .(seqname, start, end, gene_name, score = 999, strand)
    ] %>% unique
print(paste0("## Protein coding genes: ", nrow(ptgenes), " rows."))

# transcription factors only
print("## Extract transcription factors.")
tfgenes = ptgenes[
    gene_name %in% union(tflist$gene_name, tflist$hgnc_symbol),
]
print(paste0("## TF: ", nrow(tfgenes), " rows."))

# fix column names for BED format
colnames(ptgenes) = str_replace(colnames(ptgenes), "seqname", "#chrom")
colnames(tfgenes) = str_replace(colnames(tfgenes), "seqname", "#chrom")


print("## Write BED files.")
# write protein coding
write_tsv(ptgenes, protein_coding_bed)

# write transcription factors
write_tsv(tfgenes, tf_bed)
print("## Done!")



