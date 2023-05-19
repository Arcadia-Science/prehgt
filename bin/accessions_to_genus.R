library(tidyverse)

metadata <- read_tsv(snakemake@input[['metadata']]) %>%
  filter(source == "genome")

metadata <- metadata %>% 
  filter(genus == snakemake@wildcards[['genus']])

write_csv(metadata , snakemake@output[['genus']])
