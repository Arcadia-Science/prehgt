#!/usr/bin/env Rscript
library(tidyverse)
library(fastcluster)

# command line args -------------------------------------------------------

# read from command line arguments and set global variables
args <- commandArgs(trailingOnly = TRUE)
pepstats_txt <- args[1]
pepstats_hgt_out <- args[2]
gene_lst_out <- args[3]

# functions ---------------------------------------------------------------

# write a function to parse the output of emboss pepstats to get relative amino acid frequencies
parse_pepstats_to_amino_acid_frequencies <- function(content) {
  protein_ids <- c()
  amino_acid_data <- list()
  
  current_protein_id <- ""
  current_amino_acid_info <- list()
  
  for (line in content) {
    if (grepl("^PEPSTATS", line)) {
      if (current_protein_id != "") {
        current_amino_acid_info[["ID"]] <- current_protein_id
        amino_acid_data[[current_protein_id]] <- current_amino_acid_info
        current_amino_acid_info <- list()
      }
      current_protein_id <- str_extract(line, "(?<=of\\s)[^\\s]+")
      protein_ids <- append(protein_ids, current_protein_id)
    }
    
    if (grepl("^[A-Z] =", line)) {
      parts <- str_split(line, "\\s+") %>% unlist()
      amino_acid <- parts[1]
      mole_percentage <- as.numeric(parts[5])
      current_amino_acid_info[[amino_acid]] <- mole_percentage
    }
  }
  
  if (current_protein_id != "") {
    current_amino_acid_info[["ID"]] <- current_protein_id
    amino_acid_data[[current_protein_id]] <- current_amino_acid_info
  }
  
  amino_acid_df <- do.call(rbind, lapply(amino_acid_data, as.data.frame))
  return(amino_acid_df)
}

# read in data and parse --------------------------------------------------

file_content <- readLines(pepstats_txt)
raau <- parse_pepstats_to_amino_acid_frequencies(file_content)
print("RAAU parsing done.")

# create a distance matrix of genes which will estimate all-by-all amino acid usage per gene
d <- raau %>%
  select(-ID) %>% # rm bc info is in rownames
  as.matrix() %>% # convert to matrix
  dist() # compute and return a distance matrix
print("Distance matrix calculated.")

# cluster genes based on the amino acid frequency table
hr <- fastcluster::hclust(d)
print("Hierarchical clustering done.")

# define clusters by cutting the dendogram
clusters <- cutree(hr, h=max(hr$height/1.5))

# turn cluster annotations into a data frame
clusters <- clusters %>%
  as.data.frame() %>%
  select(cluster=".") %>% # rename the column to cluster
  rownames_to_column("cds") # move gene name row name to column

# filter to clusters that contain less than 150 genes
# note 150 is a totally empiric guess from limited data gazing
small_clusters <- clusters %>%
  group_by(cluster) %>%
  tally() %>%
  filter(n < 150)

# filter to genes in small clusters
clusters <- clusters %>%
  filter(cluster %in% small_clusters$cluster) %>%
  rename(hgt_candidate = cds)

# write out cluster membership file
write_tsv(clusters, pepstats_hgt_out)
# write out a gene list to use extract CDS sequences of HGT candidates from FASTA
write_tsv(clusters[1], gene_lst_out, col_names = F)

