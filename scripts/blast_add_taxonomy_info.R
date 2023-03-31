library(readr)
library(tidyr)
library(dplyr)
library(dbplyr)
library(RSQLite)

# connect to the database
lca_db <- DBI::dbConnect(RSQLite::SQLite(), snakemake@input[['sqldb']])
#lca_db <- DBI::dbConnect(RSQLite::SQLite(), "inputs/nr_cluster_taxid_formatted_final.sqlite")

# define the table to query
lca <- tbl(lca_db, snakemake@params[['sqldb_tbl']])
#lca <- tbl(lca_db, "nr_cluster_taxid_table")

# read in the BLAST results
#blast <- read_tsv("sandbox/try_blast/Schizosaccharomyces_aa_rep_seq_diamond_blastp_vs_nr_rep_seq.txt",
blast <- read_tsv(snakemake@input[['tsv']],
                   col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                                 "qstart", "qend", "sstart", "send", "evalue", "bitscore"))

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

write_tsv(blast, snakemake@output[['tsv']])
