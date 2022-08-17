suppressMessages(library(ATACseqQC))

SNAKEMAKE_MODE = T

if (SNAKEMAKE_MODE) {
  print("### Running in snakemake script.")
  bam = snakemake@input[[1]]
  outfile = snakemake@output[[1]]

} else {
  setwd("/gpfs/commons/groups/sanjana_lab/cdai/NeuronReporters/code")
  bam = "results/ATACseq/MarkDups/t0_1.bam"
  outfile = "results/ATACseq/plotFragmentSizes/t0_1.bam"
}

#estimateLibComplexity(readsDupFreq(bamFile = bam))

print(paste0("### Input file: ", bam))
pdf(outfile)

fragSize <- fragSizeDist(bamFiles = bam,
                         bamFiles.labels = sub(".bam", "", basename(bam)))

dev.off()
