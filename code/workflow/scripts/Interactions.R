# custom function to check if a motifname contains a hit gene name or TF gene name
# the reason is motifname sometimes contains "(var.2)" or other characters
mystrfunc <- function(motifname, target_match) {
    # motifname: scalar string value, name of the motif: eg. TFAP4(var.2)
    # target_match: vector, either Tf list of Hit list, or any gene list
    # returns T/F
    checkname <- function(x) {
        if (str_detect(motifname, paste0("^", x, "$"))) {
            check <- TRUE
        } else if (str_detect(motifname, paste0(":", x, "$"))) {
            check <- TRUE
        } else if (str_detect(motifname, paste0("^", x, ":"))) {
            check <- TRUE
        } else {
            check <- FALSE
        }
    }
    test <- map_lgl(target_match, ~ checkname(.x)) %>% any
    #test <- map_lgl(target_match, ~ if (str_detect(motifname, paste0("^", .x, "$"))) {TRUE} else if (str_detect(motifname, paste0(":", .x, "$"))) {TRUE} else if (str_detect(motifname, paste0("^", .x, ":"))) {TRUE} else {FALSE})
    #test <- map_lgl(target_match, ~ str_detect(string = motifname, pattern = paste0("^", .x, "$")) | str_detect(string = motifname, pattern = paste0(":", .x, "$")) | str_detect(string = motifname, pattern = paste0("^", .x, ":")))
    return(test)
}

getTargetsFromMotif <- function(motif_name, motif.matches, genomic.features, min.overlap = 10, genelist) {
    # Given a motif name (JASPAR), find its binding targets with read counts, observed in ATAC-seq
    #--------------- INPUT: -----------------
    # motif_name: string, this is the jaspar motif name, not gene name
    # motif.matches: matchMotif object, e.g. motifMatch() result of peaks and jaspar 2020 motifs
    # genomic.features: GenomicRanges object, e.g. a bed formated promoter region coordinates, converted into GRanges object
    # min.overlap: integer, minimum overlap in base pairs
    # genelist: vector, a list of gene_names. Note these gene names are some times not the same as motif names, hence the need of motif_lookup
    #--------------- OUTPUT: -----------------
    # overlapped.target.readcount: dataframe, containing read counts of genes that are targets of a given regulator. 
    # Read counts sum of all peaks that match to the regulator's motif.
    # -----------------------------------------
    
    # Get motif matches matrix (rows are peak coordianes, columns are each motif name, values are logical values indicating match)
    match.matrix <- motifMatches(motif.matches)
    
    # Get peaks that have binding sites matching given motif, each GRange also include read count columns
    GRanges.index <- which(match.matrix[ , motif_name]) 
    matched.target.GRanges <- rowRanges(motif.matches)[GRanges.index]
    
    # Find overlaps between: 
    # 1. peaks that match binding sites of a given motif, and 
    # 2. coordinates of gene annotation
    target.intersect.feature <- GenomicRanges::findOverlaps(query = matched.target.GRanges, subject = genomic.features,
                                                           minoverlap = min.overlap, ignore.strand = T)
    
    # Get the genes (from annoation) that overlap with peaks (that have binding sites of given motif)
    overlapped.gene_names <- genomic.features[subjectHits(target.intersect.feature), ] %>% as.data.frame %>% pull(gene_name) 
    
    # Summarize the read counts of matching peaks into read counts per gene for a given motif
    overlapped.target.readcount <- matched.target.GRanges[queryHits(target.intersect.feature), ] %>%
                                mcols %>% as.data.frame %>%
                                add_column("target_gene" = overlapped.gene_names) %>% # add gene_name to matched ranges
                                dplyr::select(target_gene,A1:A12) %>% # selecct read count columns only
                                filter(target_gene %in% genelist) %>% # keep rows with genes in the genelist
                                group_by(target_gene) %>% summarise_all(sum) # sum read counts per gene, as multiple peaks can fall into one gene
                                
    return(overlapped.target.readcount) # return a dataframe
}

getRegulatorsOfGene <- function(gene_name, motif.matches, genomic.features, min.overlap = 10, motif.list) {
    ###### INPUTS: #####
    # gene_name: string, gene name
    # motif.matches: matchMotif object, e.g. motifMatch() result of peaks and jaspar 2020 motifs
    # countGRanges: GenomicRanges object, e.g. a normalized read matrix, converted into GRanges object
    # genomic.features: GenomicRanges object, e.g. a bed formated promoter region coordinates, converted into GRanges object
    # min.overlap: integer, minimum overlap in base pairs
    # motif.list: vector, a list of motif names. Only return results if found motifs are part of this list
    
    ##### OUTPUTS: #####
    # regulators.readcount: list, each element's name is a regulator's motif name, the values are
    # summarised readcounts of matching/overlapping peaks for this regulator
    
    match.matrix <- motifMatches(motif.matches) # motif match logic matrix
    match.matrix.GRanges <- rowRanges(motif.matches) # peak coordiantes, including read counts in mcols
    
    # target gene's annotation coordinates
    gene.GRanges <- genomic.features[genomic.features$gene_name == gene_name]
    
    # intersect peak coordinates with target gene's annotation coordinates
    # to get peaks that belong to target gene
    match.intersect.gene <- GenomicRanges::findOverlaps(query = match.matrix.GRanges, subject = gene.GRanges,
                                                       minoverlap = min.overlap, ignore.strand = T)
    # row index of peaks that fall in target gene coordinates
    overlapped.match.index <- queryHits(match.intersect.gene) 
    
    # Once extracted peaks that fall into target gene region, get the motifs that bind to these peaks.
    # Note slight difference in extracting motif names when there are only 1 peak versus more than 1 peaks
    if (length(overlapped.match.index) == 1) { # if only one peak in target gene
        gene.regulators <- match.matrix[overlapped.match.index, ] %>% .[.] %>% names 
    } else { # if multiple peaks in target gene
        gene.regulators <- match.matrix[overlapped.match.index, ] %>% colSums %>% .[.>0] %>% names
    }
    
    # Only keep motifs that are part of a given list
    gene.regulators <- gene.regulators[gene.regulators %in% motif.list] 
    
    # Get the read counts of the binding sites for each motif (regulator)
    regulators.readcount <- list()
    for (regulator in gene.regulators) {
        readcount <- match.matrix.GRanges[overlapped.match.index, ] %>% # readcount of each peak in target gene region
                        as.data.frame %>% select(A1:A12) %>% # extract read counts columns
                        `*`(match.matrix[overlapped.match.index, regulator]) %>% # multiply 0 or 1 based on if the peak has a match to this regulator or not
                        colSums() # sum up all the reads from peaks that match to the regulator binding site
        regulators.readcount[[regulator]] <- readcount
    }
    
    # Return a dataframe, each row gives the name of the motif that binds to the target gene, 
    # along with observed read count assciated with the matching peaks. Some genes may have no matches
    if (length(regulators.readcount) > 0) {
        regulators.readcount <- do.call(rbind, regulators.readcount) %>% 
                                as.data.frame %>% 
                                add_column("target_gene" = gene_name, .before = "A1") %>%
                                rownames_to_column("regulator_motif")
        return(regulators.readcount)
    }
}

