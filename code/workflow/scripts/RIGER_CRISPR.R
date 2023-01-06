# -------------------------------------------------------------------------
#             SET UP
# -------------------------------------------------------------------------

library(purrr)
library(dplyr)
library(data.table)
library(stringr)
library(parallel)

SnakeMode = FALSE

if (SnakeMode) {
  print("### Running in snakemake script mode")
  ncores = min(snakemake@threads[[1]], detectCores())
} else {
  print("### Running in interactive/debugging mode")
  ncores = min(10, detectCores())

  library(yaml)
  config = yaml.load_file("config/config.yaml")
  l2ratio.file = config$CRISPR_LOG2_RATIOS[[1]]
  cnts.file = config$CRISPR_COUNTS
  # cnts.raw.file = 'resources/crispr/counts/wrangled_raw_counts.tsv'
  # cnts.norm.file = 'resources/crispr/counts/wrangled_ranknorm_counts.tsv'
}

Add1s = function(x) {
  if_else(x == 0, as.integer(1), x)
}

ratios = fread(l2ratio.file)
cnts = fread(cnts.file)[, .(ID, Gene, S6, S7, S8)]

cnts[, `:=`(
  S6 = Add1s(S6),
  S7 = Add1s(S7),
  S8 = Add1s(S8)
)]


cnts.norm = cnts[, c(list(ID = ID, Gene=Gene),
                map(list(S6=S6, S7=S7, S8=S8), ~ .x*1e6/sum(.x))
                )]

cnts.norm[, `:=`(
  r76 = log2(S7/S6),
  r86 = log2(S8/S6),
  r87 = log2(S8/S7)
)]


cnts.norm[]

cnts.norm = melt(cnts.norm[, -c("S6", "S7", "S8")],
     id.vars = c("ID", "Gene"),
     variable.name = "Condition",
     value.name = "l2fc") %>%
  .[, .(ID, l2fc,
        rnk = rank(-l2fc, ties.method = "first")),
    by = c("Gene", "Condition")]
#%>%
#  .[rnk < 3]

cnts.norm[, .(
  ws = map2_dbl(l2fc, rnk, ~)
), by = c("Gene", "Condition")]
