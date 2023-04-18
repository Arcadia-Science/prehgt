library(readr)
library(dplyr)
library(purrr)
library(janitor)

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

# read in outputs from compositional scans --------------------------------

compositional <- unlist(snakemake@input[['compositional']]) %>%
#compositional <- Sys.glob("outputs/compositional_scans_hgt_candidates/*_clusters.tsv") %>%
  set_names() %>%
  map_dfr(read_tsv, .id = "genus") %>%
  mutate(genus = gsub("_clusters.tsv", "", basename(genus))) %>%
  rename(RAAU_cluster = cluster)

# read in BLAST results ---------------------------------------------------

# check the number of rows for each input file
# remove the files that only have 1 row, which means no BLAST results
blast_files <- unlist(snakemake@input[['blast']])
#blast_files <- Sys.glob("outputs/blast_hgt_candidates/*_blast_scores.tsv")
blast_files_with_results <- character()
i <- 1
for(file in blast_files){
  if(length(count_fields(file, tokenizer_tsv())) > 1){
    blast_files_with_results[i] <- file
    i <- i + 1
  }
}

blast <- blast_files_with_results %>%
  set_names() %>%
  map_dfr(read_tsv, .id = "genus") %>%
  rename_with( ~ paste0("blast_", .x)) %>%
  rename(hgt_candidate = blast_qseqid, genus = blast_genus) %>%
  mutate(genus = gsub("_blast_scores.tsv", "", basename(genus)))


# read in eggnog annotations ----------------------------------------------

eggnog <- unlist(snakemake@input[['eggnog']]) %>%
#eggnog <- Sys.glob("outputs/hgt_candidates_annotation/eggnog/*.emapper.annotations") %>%
  set_names() %>%
  map_dfr(read_tsv, skip = 4, comment = "##", .id = "genus") %>%
  clean_names() %>%
  rename_with( ~ paste0("eggnog_", .x)) %>%
  rename(hgt_candidate = eggnog_number_query, genus = eggnog_genus) %>%
  mutate(genus = gsub(".emapper.annotations", "", basename(genus)))

# read in hmm annotations -------------------------------------------------

# read in annotation metadata files
hmm_descriptions <- read_tsv("https://raw.githubusercontent.com/antismash/antismash/master/antismash/detection/hmm_detection/data/hmmdetails.txt") %>%
  clean_names() %>%
  select(name = condensation, hmmscan_description = condensation_domain)
vog_descriptions <- read_tsv("http://fileshare.csb.univie.ac.at/vog/latest/vog.annotations.tsv.gz") %>%
  clean_names() %>%
  select(name = number_group_name, hmmscan_description = consensus_functional_description)
hmm_descriptions <- bind_rows(hmm_descriptions, vog_descriptions)

hmmscan <- unlist(snakemake@input[['hmmscan']]) %>%
  #hmmscan <- Sys.glob("outputs/hgt_candidates_annotation/hmmscan/*.tblout") %>%
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
  relocate(hmmscan_description, .after = hmmscan_domain_name)


# read in and parse pangenome information ---------------------------------

pangenome_size <- unlist(snakemake@input[['acc_to_genus']]) %>%
  #Sys.glob("outputs/accessions_to_genus/*.csv") %>%
  map_dfr(read_csv) %>%
  group_by(genus) %>%
  tally() %>%
  select(genus, pangenome_size = n) %>%
  ungroup()

# report number of genes in pangenome cluster
pangenome_cluster_sizes <- unlist(snakemake@input[['pangenome_cluster']]) %>%
  # Sys.glob("outputs/genus_pangenome_clustered/*_cds_cluster.tsv") %>%
  set_names() %>%
  map_dfr(read_tsv, col_names = c("rep", "member"), .id = "genus") %>%
  mutate(genus = gsub("_cds_cluster.tsv", "", basename(genus))) %>%
  mutate(rep = paste0(rep, "_1")) %>%
  filter(rep %in% c(compositional$hgt_candidate, blast$hgt_candidate)) %>%
  group_by(genus, rep) %>%
  tally() %>%
  select(genus, hgt_candidate = rep, pangenome_num_genes_in_cluster = n)


# join with gff information -----------------------------------------------

gff <- unlist(snakemake@input[['gff']]) %>%
  map_dfr(read_tsv)

hgt_candidates <- data.frame(hgt_candidate = c(blast$hgt_candidate, compositional$hgt_candidate))

gff_to_join <- hgt_candidates %>%
  # split on underscore in rep sequence name and identify the fraction of the split string that matches to the proteins IDs in the gff file 
  mutate(gff_protein_id = sapply(strsplit(hgt_candidates$hgt_candidate, "_"), 
                                 function(x) first(intersect(gff$gff_protein_id, x)))) %>%
  distinct()

gff <- gff %>%
  right_join(gff_to_join, by = c("gff_protein_id")) %>%
  distinct()

# combine and write outputs -----------------------------------------------
all_candidates <- full_join(compositional, blast, by = c("hgt_candidate", "genus"))

all_candidates <- left_join(all_candidates, eggnog, by = c("hgt_candidate", "genus")) %>%
  left_join(hmmscan, by = c("hgt_candidate", "genus")) %>%
  left_join(pangenome_cluster_sizes, by = c("hgt_candidate", "genus"))  %>%
  left_join(pangenome_size, by = "genus") %>%
  left_join(gff, by = "hgt_candidate") %>%
  # label detection method
  mutate(method = ifelse(hgt_candidate %in% blast$hgt_candidate, "blast", NA),
         method = ifelse(hgt_candidate %in% compositional$hgt_candidate, "raau", method),
         method = ifelse(hgt_candidate %in% blast$hgt_candidate & hgt_candidate %in% compositional$hgt_candidate, "both", method),
         .after = hgt_candidate)

write_tsv(all_candidates, snakemake@output[['all_results']])

method_tally <- all_candidates %>%
  group_by(genus, method) %>%
  tally()

write_tsv(method_tally, snakemake@output[['method_tally']])
