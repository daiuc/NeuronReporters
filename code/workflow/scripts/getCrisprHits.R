#' @title      Find CRISPR screen hits
#' @author     Chao Dai
#' @details    --from manuscript method section (Congyi Lu)...
#'             analyzed changes in sgRNA distribution in tdTomato- cells vs. either
#'             unsorted inputs or tdTomato+ cells... For each comparison, first
#'             assessed how many sgRNAs for each TF were enriched above the top
#'             10th percentile of non-targeting sgRNAs (empirical FDR < 0.1). With
#'             these enriched sgRNAs, TFs were sorted based on the number of gRNAs
#'             targeting each TF present. Enriched TFs with 2+ enriched sgRNAs
#'             were selected for comparison. Neuron-essential human TFs were
#'             defined as those TFs shared between the two groups.


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
  # cnts.raw.file = 'resources/crispr/counts/wrangled_raw_counts.tsv'
  # cnts.norm.file = 'resources/crispr/counts/wrangled_ranknorm_counts.tsv'
}


# -------------------------------------------------------------------------
#             Support Functions
# -------------------------------------------------------------------------

findHits = function(dt, min_hit_sgRNA = 2) {
  #' @param dt data.table of the ratios along with sgRNA IDs, and Gene name
  #' @param min_hit_sgRNA minimum number of guides passing ratio cutoff
  #' @return a vector of gene names that's a hit

  # use NTC to set up thresholds
  # r76 hits < NonTargeting's 10th percentile ratio
  # r86, r87, hits > Nontargeting's 90th percentile
  r86_cutoff = dt[Gene %in% c("NTC"), r86] %>% quantile(.9)
  r87_cutoff = dt[Gene %in% c("NTC"), r87] %>% quantile(.9)

  # set flag = T if the sgRNA pass cutoff
  dt = dt[, .(
    ID, Gene, r76, r86, r87,
    r86flag = if_else(r86 > r86_cutoff & ! Gene %in% c("NTC"), T, F),
    r87flag = if_else(r87 > r87_cutoff & ! Gene %in% c("NTC"), T, F)
  )]

  # for each gene, count number of sgRNA passing cutoff
  hits.dt = dt[! Gene %in% "NTC"][
    , .(nHits_86 = sum(r86flag),
        nHits_87 = sum(r87flag)),
    by = Gene][,
      .(Gene, nHits_86, nHits_87,
        nHits_mean = map2_dbl(nHits_86, nHits_87, ~mean(c(.x,.y))))
    ]

  hit.genes = intersect(
    hits.dt[nHits_87 >= min_hit_sgRNA, Gene],
    hits.dt[nHits_86 >= min_hit_sgRNA, Gene])

  #return(hit.genes)
  return(list(hits.dt = hits.dt, hit.genes = hit.genes))
}

randomizeRowValues = function(dt) {
  #' @param dt original ratio data.table
  #' @return row lables(ID and Gene) randomized, but r76, r86, r87 bootstrapped

  value_idx = sample(1:nrow(dt), nrow(dt), replace = T)
  all_values = c(dt$r76, dt$r86, dt$r87)
  dt = dt[, .(
    ID = ID,
    Gene = Gene,
    r76 = sample(all_values, length(ID), replace = T),
    r86 = sample(all_values, length(ID), replace = T),
    r87 = sample(all_values, length(ID), replace = T)
  )]


  return(dt)
}

# -------------------------------------------------------------------------
#             Load Data
# -------------------------------------------------------------------------


# log2 ratios of crispr counts, done by C. Lu
ratios = fread(l2ratio.file)

ratios = ratios[
  , .(ID,
    Gene = if_else(str_detect(Gene, "NonTargeting"), "NTC", Gene),
    r76, r86, r87)]


# -------------------------------------------------------------------------
#             Find hit TFs
# -------------------------------------------------------------------------

# find hits
my.hits = findHits(ratios)
names(my.hits$hit.genes) = my.hits$hit.genes


# -------------------------------------------------------------------------
#             Compute p-value using permutation
#
#     Compute the proportion of permutations of which we find the same
#     or similar number of mean_hit_guides comparing to observation from
#     true dataset.
#
# -------------------------------------------------------------------------


N_perm = 5000
perm.l = mclapply(1:N_perm,
                  function(x) {
                    randomizeRowValues(ratios) %>%
                      findHits %>%
                      .$hits.dt
                    }, mc.cores = ncores)


x = mclapply(
  my.hits$hit.genes,
  function(hitgene) {
    map_lgl(perm.l,
            ~ .x[Gene == hitgene, nHits_mean] >= my.hits$hits.dt[Gene == hitgene, nHits_mean]
            ) %>% sum
    }, mc.cores = ncores)


p = unlist(x) / N_perm









