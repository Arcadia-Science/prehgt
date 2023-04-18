library(tidyverse)

# function ----------------------------------------------------------------

parse_gff_attribute_field <- function(gff) {
  # GFF attribute fields encode a lot of information
  # each field is named and separated from its value by an equals ("=") sign
  # fields are separated from each other by a semi colon (";").
  # this function separates each attribute into a column in a data frame.
  # note that the column must be named "attribute"
  
  # split the attributes column by semicolon
  gff_split <- gff %>% 
    mutate(split_attr = str_split(string = .data$attribute, pattern = ";")) %>%
    unnest(split_attr)
  
  # Separate the key-value pairs
  gff_key_value <- gff_split %>%
    separate(split_attr, into = c("key", "value"), sep = "=", extra = "merge") %>%
    spread(key, value)  # Spread the key-value pairs to wide format
  
  return(gff_key_value)
}

# read and parse ----------------------------------------------------------

# accessions <- c("GCA_027627235.1", "GCA_013053245.1", "GCA_015484485.1")
# gffs <- paste0("inputs/genbank/", accessions, "_genomic.gff.gz")
gff <- unlist(snakemake@input[['gff']]) %>%
  map_dfr(read_tsv, comment = "##", skip = 5,
          col_names = c("seqid", "source", "feature", "start", "end", 
                        "score", "strand", "frame", "attribute"),
          col_types = "cccddcccc")

# make contig/chr length a column
seq_length <- gff %>%
  filter(feature == "region") %>%
  select(seqid, seqid_length = end)

# add column back to gff df
gff <- gff %>%
  filter(feature == "CDS") %>%
  left_join(seq_length, by = "seqid")

# parse the gff attribute field into columns in a data frame
gff <- parse_gff_attribute_field(gff)

# rename columns to label source of info as "gff"
gff <- gff %>%
  rename_with( ~ paste0("gff_", .x))

# count how many frames a protein has
gff_frame_tally <- gff %>%
  group_by(gff_protein_id) %>%
  tally() %>%
  select(gff_protein_id, gff_frame_tally = n)

# select only the first frame
gff <- gff %>%
  left_join(gff_frame_tally, by = "gff_protein_id") %>%
  group_by(gff_protein_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  distinct()

# write out info
write_tsv(gff, snakemake@output[['gff']])
