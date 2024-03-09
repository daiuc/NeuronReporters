library(googlesheets4)
library(tidyverse)

# read in from my googlesheet
atac = read_sheet("https://docs.google.com/spreadsheets/d/1G47ii0JBJMhTyZzyu2c9eTvVjlfMUAHwbx-JOc5z9rA/edit#gid=1233847355", 
                  sheet = 1)
rna = read_sheet("https://docs.google.com/spreadsheets/d/1G47ii0JBJMhTyZzyu2c9eTvVjlfMUAHwbx-JOc5z9rA/edit#gid=1439670628", 
                 sheet = 2)


arrange(atac, timepoint, rep, read) %>% 
    write_tsv("~/cdai/Neurons/config/atac.tsv")

arrange(rna, treatment, timepoint, rep, batch) %>% 
    write_tsv("~/cdai/Neurons/config/rna.tsv")

