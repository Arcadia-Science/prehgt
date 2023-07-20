#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(janitor)

# command line args -------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
# inputs:
compositional_tsv     <- args[1]
blast_kingdom_tsv     <- args[2]
blast_subkingdom_tsv  <- args[3]
genomes_csv           <- args[4]
pangenome_cluster_tsv <- args[5]
gff_tsv               <- args[6]
kofamscan_tsv         <- args[7]
hmmscan_tblout        <- args[8]
# outputs:
all_results_tsv       <- args[9]
method_tally_tsv      <- args[10]

# paths to test locally:
# compositional_tsv <- "~/github/prehgt/out_test/compositional/Bigelowiella_clusters.tsv"
# blast_kingdom_tsv <- "~/github/prehgt/out_test/blastp/Bigelowiella_blastp_kingdom_scores.tsv"
# blast_subkingdom_tsv <- "~/github/prehgt/out_test/blastp/Bielowiella_blastp_subkingdom_scores.tsv"
# genomes_csv <- "~/github/prehgt/out_test/download/Bigelowiella_genomes.csv"
# pangenome_cluster_tsv <- "~/github/prehgt/out_test/build/Bigelowiella_cluster.tsv"
# gff_tsv <- "~/github/prehgt/out_test/combine/Bigelowiella_gff_info.tsv"
# kofamscan_tsv <- "~/github/prehgt/out_test/kofamscan/Bigelowiella_kofamscan.tsv"
# hmmscan_tblout <- "~/github/prehgt/out_test/hmmscan/Bigelowiella.tblout"

# functions ---------------------------------------------------------------

read_tblout <- function(file, type){
  # modified from https://github.com/arendsee/rhmmer/blob/master/R/parse.R  
  col_types <- readr::cols(
    domain_name         = readr::col_character(),
    domain_accession    = readr::col_character(),
    query_name          = readr::col_character(),
    query_accession     = readr::col_character(),
    sequence_evalue     = readr::col_double(),
    sequence_score      = readr::col_double(),
    sequence_bias       = readr::col_double(),
    best_domain_evalue  = readr::col_double(),
    best_domain_score   = readr::col_double(),
    best_domain_bis     = readr::col_double(),
    domain_number_exp   = readr::col_double(),
    domain_number_reg   = readr::col_integer(),
    domain_number_clu   = readr::col_integer(),
    domain_number_ov    = readr::col_integer(),
    domain_number_env   = readr::col_integer(),
    domain_number_dom   = readr::col_integer(),
    domain_number_rep   = readr::col_integer(),
    domain_number_inc   = readr::col_character()
  )
  
  N <- length(col_types$cols)
  
  # the line delimiter should always be just "\n", even on Windows
  lines <- readr::read_lines(file, lazy=FALSE, progress=FALSE)
  
  table <- sub(pattern = sprintf("(%s).*", paste0(rep('\\S+', N), collapse=" +")),
               replacement = '\\1',
               x=lines,
               perl = TRUE) %>%
    gsub(pattern="  *", replacement="\t") %>%
    paste0(collapse="\n") %>%
    readr::read_tsv(col_names=names(col_types$cols), comment='#', na='-', 
                    col_types = col_types, lazy=FALSE, progress=FALSE)
  
  descriptions <- lines[!grepl("^#", lines, perl=TRUE)] %>%
    sub(pattern = sprintf("%s *(.*)", paste0(rep('\\S+', N), collapse=" +")),
        replacement = '\\1', 
        perl = TRUE)
  
  table$description <- descriptions[!grepl(" *#", descriptions, perl=TRUE)]
  
  return(table)
}

read_blast_hgt_candidates_kingdom <- function(file){
  # col_types retrieved using spec_tsv() on test file
  col_types <- cols(qseqid = col_character(),
                    hgt_taxonomy_level = col_character(),
                    acceptor_lineage_at_hgt_taxonomy_level = col_character(),
                    acceptor_lca_level = col_character(),
                    acceptor_best_nonself_match_id = col_character(),
                    acceptor_max_pident = col_double(),
                    acceptor_max_bitscore = col_double(),
                    acceptor_min_evalue = col_double(),
                    acceptor_num_matches_at_lineage = col_double(),
                    donor_num_matches_at_lineage = col_double(),
                    total_num_matches = col_double(),
                    donor_lineage_at_hgt_taxonomy_level = col_character(),
                    donor_best_match_full_lineage = col_character(),
                    donor_best_match_id = col_character(),
                    donor_best_match_pident = col_double(),
                    donor_max_bitscore = col_double(),
                    donor_min_evalue = col_double(),
                    alien_index = col_double(),
                    hgt_index = col_double(),
                    donor_distribution_index = col_double(),
                    entropy = col_double(),
                    entropy_normalized = col_double(),
                    gini = col_double(),
                    acceptor_sum_bitscore_per_group_01 = col_double(),
                    donor_sum_bitscore_per_group_01 = col_double(),
                    ahs_01_index = col_double())
  df <- read_tsv(file, col_types = col_types)
  return(df)
}

read_blast_hgt_candidates_subkingdom <- function(file){
  col_types <- cols(qseqid = col_character(),
                    hgt_taxonomy_level = col_character(),
                    acceptor_lineage_at_hgt_taxonomy_level = col_character(),
                    acceptor_lca_level = col_character(),
                    acceptor_best_nonself_match_id = col_character(),
                    acceptor_max_pident = col_double(),
                    acceptor_max_bitscore = col_double(),
                    acceptor_min_evalue = col_double(),
                    acceptor_num_matches_at_lineage = col_double(),
                    total_num_matches = col_double(),
                    donor_num_matches_at_lineage = col_double(),
                    donor_lineage_at_hgt_taxonomy_level = col_character(),
                    donor_best_match_full_lineage = col_character(),
                    donor_best_match_id = col_character(),
                    donor_best_match_pident = col_double(),
                    donor_max_bitscore = col_double(),
                    donor_min_evalue = col_double(),
                    transfer_index = col_double(),
                    transfer_index_p_value = col_double(),
                    transfer_index_adjusted_p_value = col_double())
  df <- read_tsv(file, col_types = col_types)
  return(df)
}

# read in outputs from compositional scans --------------------------------

compositional <- compositional_tsv %>%
  set_names() %>%
  map_dfr(read_tsv, col_types = "ccd", .id = "genus") %>%
  mutate(genus = gsub("_clusters.tsv", "", basename(genus))) %>%
  rename(RAAU_cluster = cluster) %>%
  distinct()

# read in BLAST results ---------------------------------------------------

# check the number of rows for each input file
# remove the files that only have 1 row, which means no BLAST results
return_blast_files_with_results <- function(files){
  blast_files_with_results <- character()
  i <- 1
  for(file in files){
    if(length(count_fields(file, tokenizer_tsv())) > 1){
      blast_files_with_results[i] <- file
      i <- i + 1
    }
  }
  return(blast_files_with_results)
}

#blast_kingdom_tsv <- "~/github/prehgt/out_test/blastp/Bigelowiella_blastp_scores.tsv"
blast_kingdom <- return_blast_files_with_results(blast_kingdom_tsv)
if(length(blast_kingdom) > 0){
  blast_kingdom <- blast_kingdom %>%
    set_names() %>%
    map_dfr(read_blast_hgt_candidates_kingdom, .id = "genus") %>%
    rename_with( ~ paste0("blast_", .x)) %>%
    rename(hgt_candidate = blast_qseqid, genus = blast_genus) %>%
    mutate(genus = gsub("_blastp_kingdom_scores.tsv", "", basename(genus))) %>%
    distinct() %>%
    mutate(blast_algorithm_type = "kingdom", .after = hgt_candidate)
}

blast_subkingdom <- return_blast_files_with_results(blast_subkingdom_tsv)
if(length(blast_subkingdom) > 0){
  blast_subkingdom <- blast_subkingdom %>%
    set_names() %>%
    map_dfr(read_blast_hgt_candidates_subkingdom, .id = "genus") %>%
    rename_with( ~ paste0("blast_", .x)) %>%
    rename(hgt_candidate = blast_qseqid, genus = blast_genus) %>%
    mutate(genus = gsub("_blastp_subkingdom_scores.tsv", "", basename(genus))) %>%
    distinct() %>%
    mutate(blast_algorithm_type = "sub-kingdom", .after = hgt_candidate)
  
}

if(!is.null(nrow(blast_kingdom)) & !is.null(nrow(blast_subkingdom))){
  blast <- bind_rows(blast_kingdom, blast_subkingdom)
} else if(!is.null(nrow(blast_kingdom))){
  blast <- blast_kingdom
} else {
  blast <- blast_subkingdom
}
# read in and parse pangenome information ---------------------------------

pangenome_size <- genomes_csv %>%
  map_dfr(read_csv, col_types = "cc") %>%
  group_by(genus) %>%
  tally() %>%
  select(genus, pangenome_size = n) %>%
  ungroup() %>%
  distinct()

# report number of genes in pangenome cluster
pangenome_cluster_sizes <- pangenome_cluster_tsv %>%
  set_names() %>%
  map_dfr(read_tsv, col_types = "cc", col_names = c("rep", "member"), .id = "genus") %>%
  mutate(genus = gsub("_cluster.tsv", "", basename(genus))) %>%
  mutate(rep = paste0(rep, "_1")) %>%
  filter(rep %in% c(compositional$hgt_candidate, blast$hgt_candidate)) %>%
  group_by(genus, rep) %>%
  tally() %>%
  select(genus, hgt_candidate = rep, pangenome_num_genes_in_cluster = n) %>%
  distinct()

# join with gff information -----------------------------------------------

gff <- gff_tsv %>%
  map_dfr(read_tsv)

hgt_candidates <- data.frame(hgt_candidate = c(blast$hgt_candidate, compositional$hgt_candidate))

gff <- gff %>%
  rowwise() %>%
  mutate(hgt_candidate = list(hgt_candidates$hgt_candidate[grepl(gff_protein_id, hgt_candidates$hgt_candidate)])) %>%
  ungroup() %>%
  unnest(hgt_candidate) %>%
  distinct()
  
# read in kofamscan annotations ---------------------------------------------

kofamscan <- kofamscan_tsv %>%
  set_names() %>%
  map_dfr(read_tsv, col_types = "cccdddc", comment = "#", .id = "genus", col_names = c("tmp", "hgt_candidate", "ko", "threshold", "score", "evalue", "ko_definition")) %>%
  select(-tmp) %>% # rm col with only asterisk
  rename_with( ~ paste0("kofamscan_", .x)) %>%
  rename(hgt_candidate = kofamscan_hgt_candidate, genus = kofamscan_genus) %>%
  mutate(genus = gsub("_kofamscan.tsv", "", basename(genus))) %>%
  # only select the best match for each hgt candidate
  group_by(genus, hgt_candidate) %>%
  arrange(desc(kofamscan_score)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  distinct()

# read in hmm annotations -------------------------------------------------

# read in annotation metadata files
hmm_descriptions <- read_tsv("https://raw.githubusercontent.com/antismash/antismash/master/antismash/detection/hmm_detection/data/hmmdetails.txt") %>%
  clean_names() %>%
  select(name = condensation, hmmscan_description = condensation_domain)
vog_descriptions <- read_tsv("http://fileshare.csb.univie.ac.at/vog/latest/vog.annotations.tsv.gz") %>%
  clean_names() %>%
  select(name = number_group_name, hmmscan_description = consensus_functional_description)
hmm_descriptions <- bind_rows(hmm_descriptions, vog_descriptions)

hmmscan <- hmmscan_tblout %>%  
  set_names() %>%
  map_dfr(read_tblout, .id = "genus") %>%
  mutate(genus = gsub(".tblout", "", basename(genus))) %>%
  group_by(genus, query_name) %>%
  # limit report to the best annotation for each gene
  slice_min(sequence_evalue) %>%
  slice_min(best_domain_evalue) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(genus, hgt_candidate = query_name, hmmscan_domain_name = domain_name,
         hmmscan_sequence_evalue = sequence_evalue, hmmscan_sequence_score = sequence_score, 
         hmmscan_sequence_bias = sequence_bias, hmmscan_best_domain_evalue = best_domain_evalue, 
         hmmscan_best_domain_score = best_domain_score) %>%
  left_join(hmm_descriptions,  by = c("hmmscan_domain_name" = "name")) %>%
  relocate(hmmscan_description, .after = hmmscan_domain_name) %>%
  distinct()

# combine all information -------------------------------------------------

all_candidates <- full_join(compositional, blast, by = c("hgt_candidate", "genus"))

# add in the rest of the information
all_candidates <- left_join(all_candidates, kofamscan, by = c("hgt_candidate", "genus")) %>%
  left_join(hmmscan, by = c("hgt_candidate", "genus")) %>%
  left_join(pangenome_cluster_sizes, by = c("hgt_candidate", "genus"))  %>%
  left_join(pangenome_size, by = "genus") %>%
  left_join(gff, by = "hgt_candidate") %>%
  # label detection method
  mutate(method = ifelse(hgt_candidate %in% blast$hgt_candidate, "blast", NA),
         method = ifelse(hgt_candidate %in% compositional$hgt_candidate, "raau", method),
         method = ifelse(hgt_candidate %in% blast$hgt_candidate & hgt_candidate %in% compositional$hgt_candidate, "both", method),
         .after = hgt_candidate) %>%
  distinct()

# label HGT candidates based on all information ---------------------------

# contamination labeling rules
# * if it’s only in one genome, mark as likely contamination if:
#    * percent identity is > 90% OR (rationale: not enough info to support that it’s real, and a long contig could be a chimeric assembly)
#    * percent identity > 70% and the contig length is less than 20kbp (it’s both short and high-ish identity, so don’t trust it)
# * if it’s in more than one genome, mark as likely contamination if:
#    * percent identity is >90% AND contig length is less than 20kbp (assumes that the same contaminant snuck into multiple genomes of the same genus)

label_contamination <- function(all_candidates_df){
  all_candidates_df <- all_candidates_df %>%
    mutate(blast_contamination = ifelse(pangenome_size == 1 & blast_donor_best_match_pident >= 90, "0 likely contamination",
                                        ifelse(pangenome_size == 1 & blast_donor_best_match_pident >= 70 & gff_seqid_length <= 20000, "0 likely contamination",
                                               ifelse(pangenome_size > 1 & pangenome_num_genes_in_cluster > 1 & blast_donor_best_match_pident >= 90 & gff_seqid_length <= 20000, "0 likely contamination", "likely not contamination"))))
  return(all_candidates_df$blast_contamination)
}

# label BLAST. Note that this logic almost exclusively relies on alien index. 
# until we have run this many times and cross checked our results with tree-based approaches, I think this is good enough for now.
# will require very good documentation to make this decision clear, and to educate around HGT index etc.
all_candidates <- all_candidates %>%
  mutate(blast_contamination = label_contamination(.)) %>%
  mutate(blast_HGT_score = ifelse(blast_alien_index >= 45, "3 highly likely HGT", NA),
         blast_HGT_score = ifelse(blast_alien_index < 45 & blast_alien_index > 15, "2 likely HGT", blast_HGT_score),
         blast_HGT_score = ifelse(blast_alien_index < 15, "1 possible HGT", blast_HGT_score),
         # relabel potential contaminants
         blast_HGT_score = ifelse(blast_contamination == "0 likely contamination", "0 likely contamination", blast_HGT_score)) %>%
  select(-blast_contamination) %>%
  relocate(blast_HGT_score, .after = blast_algorithm_type)

# write outputs -----------------------------------------------------------

write_tsv(all_candidates, all_results_tsv)

method_tally <- all_candidates %>%
  group_by(genus, method) %>%
  tally()

write_tsv(method_tally, method_tally_tsv)
