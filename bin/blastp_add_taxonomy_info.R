#!/usr/bin/env Rscript
library(readr)
library(tidyr)
library(dplyr)
library(dbplyr)
library(DBI)
library(RSQLite)

# read from command line arguments and set global variables
args <- commandArgs(trailingOnly=TRUE)
lca_db_path <- args[1]
blast_tsv <- args[2]
blast_out <- args[3]

# connect to the database
lca_db <- DBI::dbConnect(RSQLite::SQLite(), lca_db_path)
#lca_db <- DBI::dbConnect(RSQLite::SQLite(), "inputs/nr_cluster_taxid_formatted_final.sqlite")

# define the table to query
lca <- tbl(lca_db, "nr_cluster_taxid_table")

# read in the BLAST results
blast <- read_tsv(blast_tsv,
                  col_names = c("qseqid", "qtitle", "sseqid", "stitle", "pident", "approx_pident", "length", 
                                "mismatch", "gapopen", "qstart", "qend", "qlen", "qcovhsp", "sstart", "send", 
                                "slen", "scovhsp", "evalue", "bitscore", "score", "corrected_bitscore"))

# make a vector for the sequence ids to retrieve lineages from the sql db for
tmp_query <- unique(blast$sseqid)

# run the SQL query
blast_lca <- lca %>%
  select(rep, lca_lineage_named) %>%
  filter(rep %in% tmp_query) %>%
  collect() # required to tell the sql query to run and return a dataframe

# join the lineage information to the BLAST results and separate the lineage into separate columns
blast <- blast %>%
  left_join(blast_lca, by = c("sseqid" = "rep")) %>%
  separate(lca_lineage_named, into = c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species", "strain"), sep = ";")

write_tsv(blast, blast_out)
