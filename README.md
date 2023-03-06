# rehgt: Recent Eukaryotic Horizontal Gene Transfer

The goal of this repository is to develop a lightweight, automated, and scalable approach for pre-screening Eukaryotic genomes or transcriptomes for recent horizontal gene transfer (HGT).
Phylogenetic analysis involving species and gene tree construction are considered the gold standard methods for HGT detection, but these methods are computationally expensive and often don't scale to thousands of genomes.
Given the barrier to scaling, users typically need to know what genomes to include prior to analysis.
The goal of this pipeline is to identify candidate donor and acceptor genomes HGT events that can then be further investigated with tree-based approaches.

We chose to focus on "recent" HGT events, reasoning that recent HGT events likely reflect recent adaptation to an organisms environment and may highlight adaptive genes with relevant technological properties, areas of genomes that are susceptible to genetic engineering, or genes that could be beneficially transfered to other closely related organisms.
We assume that recent events will contain non-ameliorated DNA sequence and won't be shared by all members of a clade (although we are unsure what level of a taxonomic lineage constitutes a clade at the moment).

## Conceptual overview of the pipeline

This pipeline screens Eukaryotic genomes or transcriptomes for horizontal gene transfer.
At the moment, only the genome screening arm has been developed, so it is described below.
Following this description, we provide deviations we expect when we add in transcriptome processing.

The pipeline begins from a metadata table that records genome accessions and lineages.
The coding domain sequences (`_CDS_from_genomic.fna.gz`) and gene annotation files (`_gff.gz`) are downloaded and then combined (0.9 percent identity threshold) into a genus-level pseudo-pangenome (the clustering is very fast, but doesn't deal with ortholog/paralog problems -- the goal of this approach is not to identify a true pangenome, but to reduce the number of genes that are screened).
Each gene in the pangenome is then screened using composition metrics (codon usage, amino acid usage, etc.).
The results of the composition screen are then parsed to identify HGT candidates (TBD) and to remove likely contaminants (TBD).
The remaining HGT candidates are then BLASTED (DIAMOND) against the NCBI clustered non-redundant database and the results are parsed to identify genes with discrepant lineage matches (exact approach TBD).

Potential changes or additions for transcriptome data:
1. The download module will need to change, as `ncbi-genome-download` does not interface with [TSA](https://www.ncbi.nlm.nih.gov/genbank/tsa/).
2. We may need to introduce gene calling and annotation modules. [Dammit](https://github.com/dib-lab/dammit/tree/v2_staging/dammit/workflows) will be a good resource to use to decide what tools to use.
3. Different contamination screening, probably looking for Eukaryotic-specific RNA processing signals (5' cap, 3' cap, kozac sequence, etc.).
4. Integration between genomes and transcriptomes.
