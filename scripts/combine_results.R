library(readr)
library(dplyr)
library(purrr)
library(janitor)

compositional <- unlist(snakemake@input[['compositional']]) %>%
#compositional <- Sys.glob("outputs/compositional_scans_hgt_candidates/*_clusters.tsv") %>%
  set_names() %>%
  map_dfr(read_tsv, .id = "genus") %>%
  mutate(genus = gsub("_clusters.tsv", "", basename(genus))) %>%
  rename(RAAU_cluster = cluster)

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

all_candidates <- full_join(compositional, blast, by = c("hgt_candidate", "genus"))

eggnog <- unlist(snakemake@input[['eggnog']]) %>%
#eggnog <- Sys.glob("outputs/hgt_candidates_annotation/eggnog/*.emapper.annotations") %>%
  set_names() %>%
  map_dfr(read_tsv, skip = 4, comment = "##", .id = "genus") %>%
  clean_names() %>%
  rename_with( ~ paste0("eggnog_", .x)) %>%
  rename(hgt_candidate = eggnog_number_query, genus = eggnog_genus) %>%
  mutate(genus = gsub(".emapper.annotations", "", basename(genus)))

all_candidates <- left_join(all_candidates, eggnog, by = c("hgt_candidate", "genus")) %>%
  mutate(method = ifelse(hgt_candidate %in% blast$hgt_candidate, "blast", NA),
         method = ifelse(hgt_candidate %in% compositional$hgt_candidate, "raau", method),
         method = ifelse(hgt_candidate %in% blast$hgt_candidate & hgt_candidate %in% compositional$hgt_candidate, "both", method),
         .after = hgt_candidate)
write_tsv(all_candidates, snakemake@output[['all_results']])

method_tally <- all_candidates %>%
  group_by(genus, method) %>%
  tally()

write_tsv(method_tally, snakemake@output[['method_tally']])

