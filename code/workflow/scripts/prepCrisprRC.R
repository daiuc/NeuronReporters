
SNAKE_MODE = T

if (SNAKE_MODE) {
    print("### snakemake mode")
    input = snakemake@input[[1]]
    output_raw = snakemake@output[[1]]
    output_norm = snakemake@output[[2]]

} else {
    print("### interactive mode!! are you sure?!!")
    input = 'resources/crispr/counts/RC.csv'
}


library(tidyverse)
library(RNOmni)

rc = read_csv(input)
rc = select(rc, ID:Gene, S6:S8)
rc = rename(rc, sgRNA = ID,
            Day14_iN = S6,
            Day14_Tdpos = S7,
            Day14_Tdneg = S8
            )


rc.norm = mutate(rc, across(Day14_iN:Day14_Tdneg, .fns = RankNorm))


write_tsv(rc, output_raw)
write_tsv(rc.norm, output_norm)
