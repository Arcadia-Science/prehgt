#!/usr/bin/env Rscript
library(tidyverse)
source("bin/functions.R")

# command line args -------------------------------------------------------

# read from command line arguments and set global variables
args <- commandArgs(trailingOnly = TRUE)
blast_tsv <- args[1]
blast_hgt_out <- args[2]
gene_lst_out <- args[3]

# read BLAST results ---------------------------------------------------

# read in BLAST results
#blast <- read_tsv("~/github/prehgt/outputs/blast_diamond/Agrocybe_vs_clustered_nr_lineages.tsv", col_types = "ccccdddddddddddddddddccccccccc")
# blast_tsv <- "~/github/prehgt/outputs/blast_diamond/Psilocybe_vs_clustered_nr_lineages.tsv"
blast <- read_and_filter_blast_results(blast_tsv)

filter_cols <- c("superkingdom", "kingdom", "phylum", "class", "order", "family")

# at the sub-kingdom level, unclassified hits will inappropriately look different, when that's not something we can actually calculate from the provided information.
# this removes strain, species & genus hits since we don't calculate HGT down to that level, and then filters any hits that were unclassified at higher tax levels
df <- blast %>%
  # find the maximum corrected bitscore for each group; this will be a proxy for self hit, which we won't always have bc of the clustered db
  group_by(qseqid) %>% 
  mutate(max_corrected_bitscore = max(corrected_bitscore)) %>%
  ungroup() %>% 
  # now filter out hits that don't have full taxonomic information
  # # remove species and strain columns
  # select(-strain, -species, -genus) %>%
  # filter out any rows that say unclassified, as these will throw off our calculations
  filter(!if_any(any_of(filter_cols), ~ grepl("unclassified", .)))

# functions ---------------------------------------------------------------

set_qlineage <- function(df, taxonomy_level){
  # rename the taxonomy_level column (e.g. genus) to taxlevel
  # this allows the column to be consistently referred to
  df <- df %>%
    rename(taxlevel = all_of(taxonomy_level))
  
  # Infer q_lineage from lineages with most hits
  qtaxlevel <- df %>%
    count(taxlevel, sort = TRUE) %>%
    slice(1) %>%
    pull(taxlevel)
  
  qlineage <- df %>%
    filter(taxlevel == qtaxlevel) %>%
    slice(1) %>%
    select(superkingdom:taxlevel) %>%
    unlist() %>%
    paste(collapse = "; ")
  return(qlineage)
}

# Define function to calculate similarity ratio of bitscores
calc_similarity_ratio <- function(df) {
  df <- df %>% 
    mutate(similarity_ratio = corrected_bitscore / max_corrected_bitscore) %>%
    return(df)
}

calc_taxonomic_distance_lineages <- function(qlineage, slineage) {
  qsplit <- strsplit(qlineage, split = ";")[[1]]
  ssplit <- strsplit(slineage, split = ";")[[1]]
  
  # determine the first taxonomic lineage that qsplit and ssplit differ
  for (i in 1:length(qsplit)) {
    if(qsplit[i] != ssplit[i]) {
      steps <- 6 - i + 1
      # return only the first level from which the taxonomic lineages differ
      steps <- max(steps)
      return(steps/6) # there are 7 possible lineage steps, normalize to 0-1 range
    }
  }
  # if there are no differences between qsplit and ssplit, return 0
  return(0)
}

calc_taxonomic_distance <- function(df, qlineage) {
  df <- df %>%
    replace_na(list(superkingdom = "unknown",
                    kingdom = "unknown",
                    phylum = "unknown",
                    class = "unknown",
                    order = "unknown",
                    family = "unknown")) %>%
    rowwise() %>%
    mutate(slineage = paste(c(superkingdom, kingdom, phylum, class, order, family, genus), collapse = "; "),
           taxonomic_distance = calc_taxonomic_distance_lineages(qlineage, slineage))
  
  return(df)
}

calc_transfer_index <- function(df) {
  num_matches_per_qseqid <- df %>%
    group_by(qseqid) %>%
    tally() %>%
    select(qseqid, num_hits = n)
  
  intermediate <- df %>% 
    group_by(qseqid) %>%
    arrange(evalue) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    mutate(intermediate = (similarity_ratio*taxonomic_distance)/rank) %>%
    group_by(qseqid) %>%
    summarize(sum_intermediate = sum(intermediate)) %>%
    distinct() %>%
    ungroup()
  
  transfer_index <- intermediate %>%
    left_join(num_matches_per_qseqid, by = "qseqid") %>%
    mutate(transfer_index = sum_intermediate/sqrt(num_hits)) %>%
    select(qseqid, transfer_index)
  
  return(transfer_index)
}

calculate_pvalues_adjusted <- function(transfer_indices, p_adjust_method = "fdr") {
  # calculate mean and standard deviation
  mean_value <- mean(transfer_indices)
  sd_value <- sd(transfer_indices)
  
  # calculate p-values based on the normal distribution
  p_values <- pnorm(transfer_indices, mean = mean_value, sd = sd_value, lower.tail = FALSE)
  
  # apply p-value adjustment with user-specified method
  adjusted_p_values <- p.adjust(p_values, method = p_adjust_method)
  
  return(list(p_values = p_values, adjusted_p_values = adjusted_p_values))
}

# Process data
hgtfinder <- function(df, similarity_ratio_threshold, qlineage){
  df <- df %>%
    calc_similarity_ratio() %>%
    filter(similarity_ratio >= similarity_ratio_threshold) %>%
    calc_taxonomic_distance(qlineage) %>%
    calc_transfer_index() %>%
    mutate(p_value = calculate_pvalues_adjusted(transfer_index)$p_values,
           adjusted_p_value = calculate_pvalues_adjusted(transfer_index)$adjusted_p_values)
  return(df)
}

run_hgtfinder <- function(df, filter_pvalues = TRUE) {
  results_all <- data.frame()
  for(threshold in c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)){
    for(taxlevel in c("kingdom", "phylum", "class", "order", "family")){
      qlineage <- set_qlineage(df = df, taxonomy_level = taxlevel)
      results <- hgtfinder(df, similarity_ratio_threshold = threshold, qlineage = qlineage) %>%
        mutate(similarity_ratio_threshold = threshold,
               taxonomy_level = taxlevel, 
               acceptor_lineage_at_hgt_taxonomy_level = gsub("; ", ";", qlineage))
      results_all <- bind_rows(results_all, results)
    }
  }
  
  # filter pvalues
  if(filter_pvalues == TRUE){
    results_all <- results_all %>%
      filter(adjusted_p_value <= 0.01)
  }
  
  return(results_all)
}

report_lineage_at_hgt_taxonomy_level <- function(lineage, taxonomy_level){
  lineage <- gsub("; ", ";", lineage)
  lineage <- str_split(lineage, pattern = ";")[[1]]
  if(taxonomy_level == "kingdom"){
    lineage_at_hgt_taxonomy_level <- paste(lineage[1], lineage[2], sep = ";")
  } else if(taxonomy_level == "phylum"){
    lineage_at_hgt_taxonomy_level <- paste(lineage[1], lineage[2], lineage[3], sep = ";")
  } else if(taxonomy_level == "class"){
    lineage_at_hgt_taxonomy_level <- paste(lineage[1], lineage[2], lineage[3], lineage[4], sep = ";")
  } else if(taxonomy_level == "order"){
    lineage_at_hgt_taxonomy_level <- paste(lineage[1], lineage[2], lineage[3], lineage[4], lineage[5], sep = ";")
  } else if(taxonomy_level == "family"){
    lineage_at_hgt_taxonomy_level <- paste(lineage[1], lineage[2], lineage[3], lineage[4], lineage[5], lineage[6], sep = ";")
  } 
  return(lineage_at_hgt_taxonomy_level)
}

# run hgtfinder algorithms ------------------------------------------------

results <- run_hgtfinder(df)

# format results and write TSV --------------------------------------------

# select the highest taxonomic level that something is seen at
# e.g., it's HGT at kingdom, the highest level is kingdom
results_best_tax <- results %>%
  group_by(qseqid) %>%
  # this arranges things so that for each qseqid, kingdom will be at top, then phylum, then class, order, family.
  arrange(factor(taxonomy_level, levels = c("kingdom", "phylum", "class", "order", "family"))) %>%
  # take only the top hit
  slice_head(n = 1)

# identify the best non-self hit for each HGT candidate
qlineage <- set_qlineage(df = df, taxonomy_level = "family") # set the query lineage for the given set of BLAST results
tax_dist_df <- data.frame(taxonomic_distance = c(0, 1/6, 2/6, 3/6, 4/6, 5/6, 6/6),
                          taxonomy_level = c("self", "family", "order", "class", "phylum", "kingdom", "kingdom")) # make a data frame of results
top_hits_all_tax_levels <- df %>%
  filter(qseqid %in% results$qseqid) %>%            # only look at the queries that are in the results
  calc_taxonomic_distance(qlineage = qlineage) %>%  # calculates the taxonomic distance for each BLAST match, using a family level qlineage
  left_join(tax_dist_df, by = "taxonomic_distance") # join with the taxonomy label information

top_hits_donor <- top_hits_all_tax_levels %>%
  group_by(qseqid, taxonomy_level) %>%         # select the highest bitscore match for each group
  arrange(desc(corrected_bitscore)) %>%
  slice_head(n = 1) %>%
  rowwise() %>%
  mutate(donor_best_match_full_lineage = paste(superkingdom, kingdom, phylum, class, order, family, genus, species, sep = ";")) %>%
  ungroup() %>%
  select(qseqid, donor_max_bitscore = corrected_bitscore, donor_min_evalue = evalue,
         donor_best_match_pident = pident, donor_best_match_full_lineage, taxonomy_level,
         donor_best_match_id = sseqid)

# join results to best hit at given taxonomic level & report in a tsv file
results_best_tax_formatted1 <- left_join(results_best_tax, top_hits_donor, by = c("qseqid", "taxonomy_level")) %>%
  rowwise() %>%
  mutate(donor_lineage_at_hgt_taxonomy_level = report_lineage_at_hgt_taxonomy_level(lineage = donor_best_match_full_lineage,
                                                                                    taxonomy_level = taxonomy_level))

# get acceptor information ------------------------------------------------

qlineage <- gsub("; ", ";", qlineage)
blast_tmp <-  blast %>%
  # filter self hits so acceptor max scores are the same as in the kingdom script
  filter(!str_detect(string = qseqid, pattern = sseqid)) %>%
  # filter to HGT candidates
  filter(qseqid %in% results_best_tax$qseqid) %>%
  # join with HGT taxonomy level
  left_join(results_best_tax) %>%
  rowwise() %>%
  # calculate lineage information
  mutate(slineage = paste(c(superkingdom, kingdom, phylum, class, order, family, genus), collapse = "; "), 
         lineage_at_hgt_taxonomy_level = report_lineage_at_hgt_taxonomy_level(lineage = slineage, taxonomy_level = taxonomy_level)) 

# filter to hits that match at the correct taxonomic lineage (this will include unclassified hits that were removed before HGT calculation)
acceptor_tmp <- blast_tmp %>%
  filter(grepl(pattern = lineage_at_hgt_taxonomy_level, x = qlineage)) 

# pull top hit information within the acceptor group
acceptor_top_stats <- acceptor_tmp %>%
  group_by(qseqid) %>%
  arrange(desc(corrected_bitscore)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(qseqid, acceptor_max_bitscore = corrected_bitscore, acceptor_min_evalue = evalue,
         acceptor_max_pident = pident, acceptor_lineage_at_hgt_taxonomy_level = lineage_at_hgt_taxonomy_level,
         acceptor_best_nonself_match_id = sseqid)

# calculate the number of hits that have the acceptor lineage at the hgt taxonomy level
acceptor_num_hits <- acceptor_tmp %>%
  group_by(qseqid, lineage_at_hgt_taxonomy_level) %>%
  tally() %>%
  select(qseqid, acceptor_lineage_at_hgt_taxonomy_level = lineage_at_hgt_taxonomy_level, 
         acceptor_num_matches_at_lineage = n)

# from the acceptor group, calculate how wide spread the protein is among other genomes in the acceptor group
acceptor_lca <- acceptor_tmp %>% 
  lca()

# combine the information together
acceptor_information <- left_join(acceptor_top_stats, acceptor_num_hits) %>%
  left_join(acceptor_lca)

# get donor number of hits ------------------------------------------------

donor_num_hits <- blast_tmp %>%
  left_join(results_best_tax_formatted) %>%
  filter(lineage_at_hgt_taxonomy_level == donor_lineage_at_hgt_taxonomy_level) %>%
  group_by(qseqid, donor_lineage_at_hgt_taxonomy_level) %>%
  tally() %>%
  select(qseqid, donor_lineage_at_hgt_taxonomy_level, donor_num_matches_at_lineage = n)

# get total number of hits ------------------------------------------------

num_hits_per_qseqid <- blast_tmp %>%
  group_by(qseqid) %>%
  tally() %>%
  select(qseqid, total_num_matches = n)

# combine and write -------------------------------------------------------

results_best_tax_formatted <- results_best_tax_formatted1 %>%
  left_join(acceptor_information) %>%
  left_join(num_hits_per_qseqid) %>%
  left_join(donor_num_hits) %>%
  select(qseqid, hgt_taxonomy_level = taxonomy_level, acceptor_lineage_at_hgt_taxonomy_level, 
         acceptor_lca_level, acceptor_best_nonself_match_id,
         acceptor_max_pident, acceptor_max_bitscore, acceptor_min_evalue,
         acceptor_num_matches_at_lineage, donor_num_matches_at_lineage, total_num_matches,
         donor_lineage_at_hgt_taxonomy_level, donor_num_matches_at_lineage,
         donor_best_match_full_lineage, donor_best_match_id, donor_best_match_pident,
         donor_max_bitscore, donor_min_evalue, transfer_index, transfer_index_p_value = p_value, 
         transfer_index_adjusted_p_value = adjusted_p_value)
         #acceptor_sum_bitscore_per_group_01, donor_sum_bitscore_per_group_01, ahs_01_index)

write_tsv(results_best_tax_formatted, blast_hgt_out)
write_tsv(results_best_tax_formatted %>% select(qseqid) %>% pull, gene_lst_out)
