library(tidyverse)


# functions used by multiple scripts --------------------------------------

read_and_filter_blast_results <- function(tsv_file){
  blast <- read_tsv(tsv_file, col_types = "ccccdddddddddddddddddccccccccc")
  
  # edit kingdom to seven categories 
  blast <- blast %>%
    mutate(kingdom = ifelse(kingdom == "unclassified Bacteria kingdom", "Bacteria", kingdom),
           kingdom = ifelse(kingdom %in% c('Bamfordvirae', 'unclassified Viruses kingdom', 'Heunggongvirae', 'Orthornavirae'), "Virus", kingdom),
           kingdom = ifelse(kingdom == "unclassified Archaea kingdom", "Archaea", kingdom),
           kingdom = ifelse(kingdom == "unclassified Eukaryota kingdom", "Other Eukaryota", kingdom)) %>%
    # filter out unclassified. Matches are to synthetic constructs that don't have taxonomies, or to clusters that don't have a superkingdom LCA (taxid 0)
    filter(!kingdom %in% c("unclassified unclassified entries kingdom",
                           "unclassified other entries kingdom",
                           "unclassified root kingdom",
                           "unclassified cellular organisms kingdom"))
  
  # filter the BLAST matches
  blast <- blast %>%
    # aerolysin example: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3424411/
    # based on aerolysin example: if corrected_bitscore is less than 100, make sure the length of the match is >70% of the original protein
    # filter out matches with a low bitscore that don't also have a high coverage
    mutate(keep = ifelse(corrected_bitscore >= 100, "keep",
                         ifelse(qcovhsp >=0.7, "keep", "filter"))) %>%
    filter(keep == "keep") %>%
    select(-keep)
  
  return(blast)
}

# blastp kingdom-level HGT prediction functions ---------------------------

ahs_index <- function(sum_bitscore_donor, sum_bitscore_acceptor){
  diff <- sum_bitscore_donor - sum_bitscore_acceptor
  return(diff)
}

lca <- function(lineage_df){
  # given a data frame with lineage columns kingdom, phylum, class, order,
  # genus, and species, calculate the taxonomy level that the lowest common 
  # ancestor of all observations.
  lineage_df <- lineage_df %>%
    # group by query protein id
    group_by(qseqid) %>%
    # count how many of each taxonomic lineage level was observed
    mutate(n_distinct_kingdom = n_distinct(kingdom),
           n_distinct_phylum = n_distinct(phylum),
           n_distinct_class = n_distinct(class),
           n_distinct_order = n_distinct(order),
           n_distinct_genus = n_distinct(genus),
           n_distinct_species = n_distinct(species)) %>%
    ungroup() %>%
    # keep only the query id and tallying columns
    select(qseqid, starts_with("n_distinct")) %>%
    # filter to distinct
    distinct() %>%
    # figure out which level the LCA occurs at
    mutate(acceptor_lca_level = ifelse(n_distinct_species == 1, "species",
                                       ifelse(n_distinct_genus == 1, "genus",
                                              ifelse(n_distinct_order == 1, "order",
                                                     ifelse(n_distinct_class == 1, "class",
                                                            ifelse(n_distinct_phylum == 1, "phylum", "kingdom")))))) %>%
    select(qseqid, acceptor_lca_level)
  
  return(lineage_df)
}

normalized_bitscore_01 <- function(bitscore, max_bitscore){
  # calculate a normalized bitscore on a 0-1 scale
  # bitscore: the bitscore to normalize
  # max_bitscore: the maximum bitscore observed for that query
  normalized_bitscore <- bitscore / max_bitscore
  return(normalized_bitscore)
}

