library(tidyverse)


hyper_enrich2 <- function(test.set, enrichment.set, pop.set) {
    # population: all black & white balls; all genes
    # m: all white balls; all hit genes
    # n: all black balls; all non-hit genes
    # q: number of white balls drawn; number of hit genes found in a particular disease set
    # k: number of balls drawn; number of genes in a particular disease set

    q <- intersect(test.set, enrichment.set) %>% length() # hit genes that are also in disease set genes
    m <- length(unique(enrichment.set)) # all hit genes
    n <- length(unique(pop.set)) - m # all non-hit genes
    k <- length(unique(test.set)) # all genes in disease set

    return(phyper(q = q - 1, m = m, n = n, k = k, lower.tail = F))
}


ggplot(mpg) +
    geom_histogram(aes(cty))

View(mpg)
?str_detect
