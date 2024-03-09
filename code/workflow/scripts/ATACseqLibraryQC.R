suppressMessages(library(ATACseqQC))

SNAKEMAKE_MODE = T

if (SNAKEMAKE_MODE) {
  print("### Running in snakemake script.")
  bam = snakemake@input[[1]]
  bamlabel = snakemake@params[[1]]
  outfile = snakemake@output[[1]]

} else {
  setwd("/gpfs/commons/groups/sanjana_lab/cdai/NeuronReporters/code")
  bam = "results/ATACseq/MarkDups/t0_1.bam"
  outfile = "results/ATACseq/plotFragmentSizes/t0_1.bam"
}

#estimateLibComplexity(readsDupFreq(bamFile = bam))

cat(paste0("### Input file: ", bam, " , label: ", bamlabel))
pdf(outfile)

fragSize <- fragSizeDist(bamFiles = bam,
                         bamFiles.labels = bamlabel
                        )

dev.off()
