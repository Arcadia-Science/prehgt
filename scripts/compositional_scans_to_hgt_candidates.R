library(tidyverse)

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

file_content <- readLines(snakemake@input[['raau']])
#file_content <- readLines("sandbox/emboss/tmp.pepstats")
raau <- parse_amino_acid_frequencies(file_content)

# create a distance matrix of genes which will estimate all-by-all amino acid usage per gene
d <- raau %>%
  select(-ID) %>% # rm bc info is in rownames
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
# note 150 is a totally empiric guess from limited data gazing
small_clusters <- clusters %>%
  group_by(cluster) %>%
  tally() %>%
  filter(n < 150)

# filter to genes in small clusters
clusters <- clusters %>%
  filter(cluster %in% small_clusters$cluster) %>%
  rename(hgt_candidate = cds)

# write out a gene list to use extract CDS sequences of HGT candidates from FASTA
write_tsv(clusters[2], snakemake@output[['gene_lst']], col_names = F)
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

