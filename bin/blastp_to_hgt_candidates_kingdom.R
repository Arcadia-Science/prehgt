#!/usr/bin/env Rscript
library(tidyverse)
source("bin/functions.R")

# command line args -------------------------------------------------------

# read from command line arguments and set global variables
args <- commandArgs(trailingOnly = TRUE)
blast_tsv <- args[1]
blast_hgt_out <- args[2]
gene_lst_out <- args[3]

# definitions -------------------------------------------------------------

# * acceptor group: The kingdom of the BLASTP query protein
# * donor group: For this script, the donor group is defined as six of seven NCBI kingdoms.
#   (bacteria, fungi, archaea, viruses, metazoans, plants, and other eukaryote species). 
#   The held-out kingdom is the kingdom of the acceptor group

# functions ---------------------------------------------------------------

# define a function to calculate the AI value
alien_index <- function(e_donor, e_acceptor) {
  # * alien_index: "An AI > 0 indicates a better hit to candidate donor than 
  #   recipient taxa and a possible HGT. Higher AI represent higher gap of E-values 
  #   between candidate donor and recipient and a more likely HGT" (10.3390/genes8100248).
  #   Because it is based on e_value, the alien index is affected by the size of 
  #   the BLAST database.
  # * e_donor: best (lowest) e_value among all BLAST matches for a given donor group
  # * e_acceptor: best (lowest) e_value among all BLAST matches for the acceptor group
  ai <- log(e_acceptor + 1e-200) - log(e_donor + 1e-200)
  return(ai)
}

# define a function to calculate the percent difference between two bitscores
hgt_index <- function(bitscore_donor, bitscore_acceptor) {
  # * hgt_index: Should not be affected by size of the BLAST database.
  # * bitscore_donor: best (highest) bitscore_value among all BLAST matches for a given donor group
  # * bitscore_acceptor: best (highest) bitscore_value among all BLAST matches for the acceptor group
  # the 200 is a simplification of /2 * 100
  diff <- (bitscore_donor - bitscore_acceptor) * 200 / (bitscore_donor + bitscore_acceptor)
  return(diff)
}

# define a function to calculate how distributed a protein is among different kingdoms
donor_distribution_index <- function(n, num_matches_per_group){
  # * donor distribution_index: a measure of the distribution of the gene across the 
  #   tree of life. Excluding matches within the acceptor group, the top 50 BLAST
  #   matches should belong to species in the actual donor group. A distribution_index > 0.8
  #   "indicates that most of the top 50 hits belong to a given donor group" (10.1016/j.molp.2022.02.001).
  #   Adapted from the tissue specificity index developed in 10.1093/bioinformatics/bti042.
  # * n: the number of donor groups
  # * num_matches_per_group: number of species per donor group. 
  #   Formatted as a vector of the number of matches for all of the donor group.
  #   If the length of the vector is < n, the function will back-add zeroes, 
  #   indicating no observations for some donor groups. 
  # * max_num_matches: maximum number of matches per group among all groups
  
  # if the length of the input vector given is less than n, backfill it with zeros (e.g. no observations for some donor groups)
  if(length(num_matches_per_group) < n){
    backadd <- rep(x = 0, times = n - length(num_matches_per_group))
    num_matches_per_group <- c(num_matches_per_group, backadd)
  }
  max_num_matches <- max(num_matches_per_group)
  donor_distribution_index <- sum(1 - (num_matches_per_group/max_num_matches)) / (n - 1)
  return(donor_distribution_index)
}

gini <- function(x) {
  # calculate gini coefficient
  x_sorted <- sort(x)
  n <- length(x_sorted)
  return (1 - 2 * (sum((1:n) * x_sorted) / sum(x_sorted) - (n + 1) / 2) / n)
}

# add more (better?) indices that report on how distributed a protein is across groups
group_specificity_indices <- function(df, kingdoms) {
  # * df: a data frame with qseqid and kingdom columns
  # * kingdoms: a vector of kingdom values, e.g. c("Fungi", "Viridiplantae", "Bacteria", "Other Eukaryota", "Metazoa", "Archaea", "Virus")
  # * entropy: a measure of the uncertainty or randomness of a set of probabilities. 
  #   Gives a higher number for a qseqid that is evenly distributed across all kingdoms.
  #   Gives a lower number for a qseqid that is mostly associated with a single kingdom.
  # * entropy_normalized: the maximum value occurs when all kingdoms are equally represented. 
  #   The entropy is -log2(1/K), where K is the number of kingdoms. 
  #   Given 7 kingdoms, the max entropy would be log2(7).
  #   Entropy is normalized to 0-1 range by dividing by the maximum potential value.
  #   A normalized entropy close to 0 indicates that the qseqid is mostly associated with one kingdom.
  #   A normalized entropy close to 1 indicates that the qseqid is uniformly distributed across all kingdoms.
  # * gini: Gini Coefficient, a measure of inequality among values of a frequency distribution. 
  #   A Gini coefficient close to 0 indicates that the qseqid is uniformly distributed across all kingdoms.
  #   A Gini coefficient close to 1 indicates that the qseqid is specific to one kingdom.
  
  # calculate the counts for each qseqid and kingdom
  df <- df %>%
    count(qseqid, kingdom) %>%
    ungroup()
  
  # include missing kingdoms, even if the count for that group is 0
  df <- expand_grid(qseqid = unique(df$qseqid), kingdom = kingdoms) %>%
    left_join(df, by = c("qseqid", "kingdom")) %>%
    replace_na(list(n = 0))
  
  # calculate specificity (fraction of matches to a specific group)
  df <- df %>%
    group_by(qseqid) %>%
    mutate(specificity = n / sum(n)) %>%
    ungroup()
  
  # use specificity to calculate entropy and counts to calculate gini coefficient
  all <- df %>%
    group_by(qseqid) %>%
    mutate(entropy = -sum(specificity * log2(specificity), na.rm = TRUE),
           entropy_normalized = entropy / log2(n_distinct(kingdom))) %>%
    mutate(gini = gini(n)) %>%
    ungroup() %>%
    select(qseqid, entropy, entropy_normalized, gini) %>%
    distinct()
  
  return(all)
}

# read BLAST results ---------------------------------------------------

# read in BLAST results
#blast_tsv <- "~/github/2023-rehgt/outputs/blast_diamond/Psilocybe_vs_clustered_nr_lineages.tsv"
blast <- read_and_filter_blast_results(blast_tsv)

# set acceptor and donor groups -------------------------------------------

# figure out the donor and acceptor groups from the BLAST results
groups <- c('Viridiplantae', 'Bacteria', 'Other Eukaryota', 'Metazoa', 'Archaea', 'Virus', 'Fungi')
# tally the total number of blast hits per kingdom in the full set of results
kingdom_tally <- blast %>% 
  group_by(kingdom) %>%
  tally()

# set the value of the acceptor group to the group with the highest number of BLAST matches
acceptor_group <- kingdom_tally[which.max(kingdom_tally$n), 1]$kingdom
# set the value of the donor groups to everything except the acceptor group
donor_groups <- groups[!groups %in% acceptor_group]  

# parse BLAST results -----------------------------------------------------

blast <- blast %>%
  # filter out exact matches as these should not be used for calculation of indices. 
  # bc of clustering, not every input protein will have a perfect match in the database
  # This is a hacky way of doing this -- it looks for whether the nr match is in the query header. 
  # The CDS from genbank should have these.
  # I like it better than filtering on pident/evalue bc those might be legitimate matches
  filter(!str_detect(string = qseqid, pattern = sseqid)) %>%
  # filter out matches to groups outside of the defined donor/acceptor groups
  # and those with missing values
  filter(kingdom %in% c(donor_groups, acceptor_group)) 

# calculate how specific matches are to a given kingdom using gini and entropy
group_specificity <- blast %>%
  group_specificity_indices(kingdoms = groups)
  
# calculate the max_bitscore and the min_evalue for each donor group and the acceptor group per query
best_match_per_group_values <- blast %>%
  group_by(qseqid, kingdom) %>%
  summarise(max_bitscore = max(corrected_bitscore),
            min_evalue = min(evalue))

# retrieve the best match and its pident and lineage for each group
best_match_per_group_identities <- blast %>%
  group_by(qseqid, kingdom) %>%
  slice_max(corrected_bitscore) %>%
  slice_min(evalue) %>%
  slice_max(pident) %>%
  slice_head(n = 1) %>% # arbitrarily select the first if there are multiple
  mutate(best_match_lineage = paste(superkingdom, kingdom, phylum, class,
                                    order, family, genus, species, sep = ";")) %>%
  select(qseqid, kingdom, best_match = sseqid, best_match_lineage, best_match_pident = pident)

# calculate the number of matches observed per group per query
num_matches_per_group <- blast %>%
  group_by(qseqid, kingdom) %>% 
  tally() %>%
  select(qseqid, kingdom, num_matches_per_group = n)

combined <- left_join(best_match_per_group_values, best_match_per_group_identities, by = c("qseqid", "kingdom")) %>%
  left_join(num_matches_per_group, by = c("qseqid", "kingdom"))

# calculate HGT indices that rely on a single value and predict candidate CDSs ------------------------

# filter out BLAST results where all matches were within the acceptor group
only_acceptor_group_matches <- combined %>%
  group_by(qseqid) %>%
  tally() %>%
  filter(n == 1)

candidates <- combined %>%
  filter(!qseqid %in% only_acceptor_group_matches$qseqid)

# separate candidates into donor and acceptor groups and reformat, and then re-join together to calculate indices
candidates_acceptor <- candidates %>%
  filter(kingdom == acceptor_group) %>%
  select(qseqid, acceptor_lineage_at_hgt_taxonomy_level = kingdom, 
         acceptor_max_bitscore = max_bitscore,
         acceptor_min_evalue = min_evalue, 
         acceptor_num_matches_at_lineage = num_matches_per_group,
         acceptor_max_pident = best_match_pident,
         acceptor_best_nonself_match_id = best_match)

candidates_donor <- candidates %>%
  filter(kingdom %in% donor_groups) %>%
  select(qseqid, donor_lineage_at_hgt_taxonomy_level = kingdom, 
         donor_max_bitscore = max_bitscore,
         donor_min_evalue = min_evalue, 
         donor_num_matches_at_lineage = num_matches_per_group,
         donor_best_match_id = best_match, 
         donor_best_match_full_lineage = best_match_lineage, 
         donor_best_match_pident = best_match_pident)

# from candidate donor groups, calculate the donor distribution index
donor_dist_index <- candidates_donor %>%
  group_by(qseqid) %>%
  mutate(donor_distribution_index = donor_distribution_index(n = length(donor_groups), 
                                                             num_matches_per_group = donor_num_matches_at_lineage))

# rejoin information together to calculate alien index and horizontal gene transfer index
candidates <- left_join(candidates_acceptor, candidates_donor, by = "qseqid", multiple = "all")

# calculate alien index and hgt index
candidates <- candidates %>%
  mutate(alien_index = alien_index(e_donor = donor_min_evalue, e_acceptor = acceptor_min_evalue),
         hgt_index = hgt_index(bitscore_donor = donor_max_bitscore, bitscore_acceptor = acceptor_max_bitscore))

# from the acceptor group, calculate how wide spread the protein is among other genomes in the acceptor group
acceptor_lca <- blast %>% 
  filter(qseqid %in% candidates_acceptor$qseqid) %>%
  filter(kingdom %in% acceptor_group) %>%
  lca()

# calculate aggregate indicies --------------------------------------------

# normalize bitscores and then sum up bitscores within each group
sum_bitscore <- blast %>%
  # remove queries that are only fungi
  filter(!qseqid %in% only_acceptor_group_matches$qseqid) %>%
  # normalize bitscores
  group_by(qseqid) %>%
  mutate(max_bitscore = max(corrected_bitscore),
         normalized_bitscore_01 = normalized_bitscore_01(corrected_bitscore, max_bitscore)) %>%
  ungroup() %>%
  # sum over normalized bitscores within groups (ingroup, outgroups)
  group_by(qseqid, kingdom) %>%
  summarize(sum_bitscore_per_group_01 = sum(normalized_bitscore_01))

# separate out acceptor group sum bitscore values
acceptor_sum_bitscore <- sum_bitscore %>%
  filter(kingdom %in% acceptor_group) %>%
  select(qseqid,
         acceptor_lineage_at_hgt_taxonomy_level = kingdom,
         acceptor_sum_bitscore_per_group_01 = sum_bitscore_per_group_01)

# separate out donor group sum bitscore values
donor_sum_bitscore <- sum_bitscore %>%
  filter(!kingdom %in% acceptor_group) %>%
  select(qseqid, 
         donor_lineage_at_hgt_taxonomy_level = kingdom,
         donor_sum_bitscore_per_group_01 = sum_bitscore_per_group_01)  

# calculate AHS using 0-1 normalized bitscores (departure from 10.1371/journal.pcbi.1010686)
ahs <- left_join(acceptor_sum_bitscore, donor_sum_bitscore, by = "qseqid") %>%
  mutate(ahs_01_index = ahs_index(sum_bitscore_donor = donor_sum_bitscore_per_group_01, sum_bitscore_acceptor = acceptor_sum_bitscore_per_group_01))

# combine results ---------------------------------------------------------

# join everything together
candidates <- candidates %>%
  left_join(acceptor_lca, by = "qseqid") %>%
  left_join(donor_dist_index, by = c("qseqid", "donor_lineage_at_hgt_taxonomy_level", "donor_max_bitscore", 
                                     "donor_min_evalue", "donor_num_matches_at_lineage",
                                     "donor_best_match_id", "donor_best_match_full_lineage",
                                     "donor_best_match_pident")) %>%
  left_join(group_specificity, by = "qseqid") %>%
  left_join(ahs, by = c("qseqid", "acceptor_lineage_at_hgt_taxonomy_level", "donor_lineage_at_hgt_taxonomy_level")) %>%
  mutate(hgt_taxonomy_level = "kingdom") # explicitly label hgt taxonomy level

# predict candidate HGT events based on results ---------------------------

# reorder outputs
candidates <- candidates %>%
  select(qseqid, hgt_taxonomy_level, acceptor_lineage_at_hgt_taxonomy_level, 
         acceptor_num_matches_at_lineage, acceptor_lca_level, acceptor_best_nonself_match_id,
         acceptor_max_pident, acceptor_max_bitscore, acceptor_min_evalue,
         donor_lineage_at_hgt_taxonomy_level, donor_num_matches_at_lineage,
         donor_best_match_full_lineage, donor_best_match_id, donor_best_match_pident,
         donor_max_bitscore, donor_min_evalue,   
         alien_index, hgt_index, donor_distribution_index, entropy, entropy_normalized, 
         gini, acceptor_sum_bitscore_per_group_01, donor_sum_bitscore_per_group_01, ahs_01_index)

# filter to genes that have the potential to be HGT events. 
# * alien index > 0 "indicates a better hit to candidate donor than recipient taxa and a possible HGT" (10.3390/genes8100248).
# * ahs > 0 "a positive AHS score suggests a potential HGT candidate" (10.1371/journal.pcbi.1010686);
#   While the paper used a more sophisticated bitscore normalization method, we found that 0-1 normalization and looking at positive scores identified what look like real HGT candidates
candidates <- candidates %>%  
  filter(alien_index > 0 | ahs_01_index > 0 | hgt_index > 0) 

write_tsv(candidates, blast_hgt_out)
write_tsv(candidates[ , 1], gene_lst_out, col_names = FALSE)
