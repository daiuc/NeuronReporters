library(tidyverse)
library(DiffBind)


RUN_MODE = 1

if (RUN_MODE == 1) {
    diffbind.samplesheet = snakemake@params[['SAMPLESHEET']]
    out.bed = snakemake@output[['bed']]
    out.counts = snakemake@output[['counts']]
    MIN_OVERLAP = snakemake@params[['MIN_OVERLAP']]

    print(paste0("Using sample sheet: ", diffbind.samplesheet))
    print(paste0("output: ", out.bed, " , ", out.counts))
    print(paste0("DiffBind using min overlap of ", MIN_OVERLAP, " samples"))
}

# ATAC metadata
diffbind.samplesheet = data.table::fread(diffbind.samplesheet)

# create a DiffBind object, including bams and macs2 peak files
peaks = dba(sampleSheet = diffbind.samplesheet, peakCaller = "macs")

# get counts, but the point here is to construct consensus peakset
consensus_peaks = dba.count(peaks, minOverlap = MIN_OVERLAP)

consensus_peaks.matrix = consensus_peaks$peaks[[1]] %>% 
    filter(str_detect(Chr, "chr[0-9XY]+"))

# export consensus peaks coordinates into a bed file
select(consensus_peaks.matrix, Chr, Start, End) %>%
    data.table::fwrite(file = out.bed, quote = F, sep = "\t", row.names = F, col.names = F)

# construct count matrix for the peaks
consensus_peaks.counts = map(consensus_peaks$peaks, ~ pull(.x, Reads))
names(consensus_peaks.counts) = consensus_peaks$samples$SampleID
consensus_peaks.counts = do.call(bind_cols, consensus_peaks.counts)

# export consensu peaks read counts to a tab file
data.table::fwrite(x = consensus_peaks.counts, file = out.counts,
                   quote = F, sep = "\t", row.names = F, col.names = T)



### --------------------------------------------------------------------
# using minOverlap = x, results in y consensus peaks
# 2 : 107800
# 3 : 70783
# 4 : 54530
# 5 : 44036
# 6 : 36834
# 7 : 31531
# 8 : 27295
# 9 : 23600
# 10 : 20287
# 11 : 16775
# 12 : 13413