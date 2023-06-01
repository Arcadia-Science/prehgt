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

pepstats_txt <- "~/github/2023-rehgt/outputs/compositional_scans_pepstats/Ophiocordyceps_pepstats.txt"
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

# calculate cluster threshold
min_cluster_size <- round(0.001 * nrow(raau)) # clusters must contain 0.1% of total proteins or fewer

small_clusters <- clusters %>%
  group_by(cluster) %>%
  tally() %>%
  filter(n <= min_cluster_size)

# filter to genes in small clusters
clusters <- clusters %>%
  filter(cluster %in% small_clusters$cluster) %>%
  rename(hgt_candidate = cds)

# write out cluster membership file
write_tsv(clusters, pepstats_hgt_out)
# write out a gene list to use extract CDS sequences of HGT candidates from FASTA
write_tsv(clusters[1], gene_lst_out, col_names = FALSE)


# try new outlier detection -----------------------------------------------

tmp <- raau %>% 
  select(-ID) %>% 
  select_if(negate(function(col) is.numeric(col) && sum(col) == 0)) %>%
  as.matrix()

distances <- mahalanobis(tmp, center = colMeans(tmp), cov = cov(tmp))
threshold <- qchisq(.99999999, df = ncol(tmp))
outlier_indices <- which(distances > threshold)
length(distances)
length(outlier_indices)

# PCA ---------------------------------------------------------------------

# Perform PCA
pca_res <- prcomp(tmp, center = TRUE, scale. = TRUE)

# Get the scores of the first two principal components
scores <- pca_res$x[,1:2]

# Plot the scores
# plot(scores, xlab = "PC1", ylab = "PC2", main = "PCA plot of RAAU")

# Identify potential HGT genes based on standard deviation threshold
outlier_threshold <- 2 # You might need to adjust this threshold
outliers <- row.names(tmp)[(scores[,1] > outlier_threshold*sd(scores[,1])) | 
                             (scores[,2] > outlier_threshold*sd(scores[,2]))]

length(outliers)

# compare three methods ---------------------------------------------------
# compare against blast as well
blast <- read_tsv("~/github/2023-rehgt/outputs/hgt_candidates_final/results_fungi.tsv") %>%
  filter(method %in% c("blast", "both")) %>%
  filter(genus == "Ophiocordyceps")

library(sourmashconsumr)
upset_df <- sourmashconsumr::from_list_to_upset_df(list(pca = outliers,
                                                        mahal = names(outlier_indices),
                                                        clusters = clusters$hgt_candidate,
                                                        blast = blast$hgt_candidate))
UpSetR::upset(upset_df)


shared <- outliers[outliers %in% names(outlier_indices)]
shared <- shared[shared %in% clusters$hgt_candidate]
length(shared)
# get protein accessions
tmp <- data.frame(shared = shared) %>%
  separate(shared, into = c("chr", "cds", "accession", "end", "start"), sep = "_")
cat(tmp$accession, sep = "\n")

