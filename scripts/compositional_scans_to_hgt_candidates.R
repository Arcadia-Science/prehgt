library(readr)
library(dplyr)
library(janitor)
library(tibble)
library(tidyr)

# read in the relative amino acid usage data frame produced by codonw and fix the column names and gene name
codonw_raau <- read_delim(snakemake@input[['raau']], delim= " ") %>%
  clean_names()  %>% # standardize column name capitalization and remove spaces
  dplyr::mutate(gene_name = gsub("_$", "", gene_name)) # remove trailing underscore introduced by codonw

# remove frequencies of unknown amino acids and create a distance matrix of genes which will estimate all-by-all amino acid usuage per gene
d <- codonw_raau %>%
  select(-unk) %>% # remove unknown codon frequencies
  column_to_rownames("gene_name") %>% # move the gene name to rows
  as.matrix() %>% # convert to matrix
  dist() # compute and return a distance matrix

# cluster genes based on the amino acid frequency table
hr <- hclust(d)

# define clusters by cutting the dendogram
clusters <- cutree(hr, h=max(hr$height/1.5))

# turn cluster annotations into a data frame
clusters <- clusters %>%
  as.data.frame() %>%
  select(cluster=".") %>% # rename the column to cluster
  rownames_to_column("cds") # move gene name row name to column

# filter to clusters that contain less than 150 genes
small_clusters <- clusters %>%
  group_by(cluster) %>%
  tally() %>%
  filter(n < 150)

# filter to genes in small clusters
clusters <- clusters %>%
  filter(cluster %in% small_clusters$cluster)

# write out a gene list to use extract CDS sequences of HGT candidates from FASTA
write_tsv(clusters[1], snakemake@output[['gene_lst']])
# write out cluster membership file
write_tsv(clusters, snakemake@output[['tsv']])

# commented out plotting code ---------------------------------------------

# # get a color palette equal to the number of clusters
# clusterCols <- rainbow(length(unique(mycl)))
# 
# # create vector of colors for side bar
# myClusterSideBar <- clusterCols[mycl]
# 
# # choose a color palette for the heat map
# myheatcol <- rev(gplots::redgreen(75))
# 
# # draw the heat map
# gplots::heatmap.2(d %>% as.matrix, main="Hierarchical Cluster", Rowv=as.dendrogram(hr), 
#                   Colv=NA, dendrogram="row", scale="row", col=myheatcol, 
#                   density.info="none", trace="none", RowSideColors= myClusterSideBar)

