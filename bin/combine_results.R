#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(purrr)

# command line args -------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
# output
out_tsv     <- args[1]
in_tsvs     <- args[2]:args[length(args)]


# read, combine, and write results files ----------------------------------

all_results <- in_tsvs %>%
  map_dfr(read_tsv)

write_tsv(all_results, out_tsv)
