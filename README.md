# prehgt: generating a *pre*liminary list of horizontal gene transfer (*HGT*) candidates using compositional and phylogenetic implicit approaches 

preHGT a pipeline for lightweight, automated, and scalable approach for pre-screening genomes for horizontal gene transfer (HGT).

TODO: add pipeline overview DAG/figure.


Phylogenetic analysis involving species and gene tree construction are considered the gold standard methods for HGT detection, but these methods are computationally expensive and often don't scale to thousands of genomes.
Given the barrier to scaling, users typically need to know what genomes to include prior to analysis.
The goal of this pipeline is to identify candidate donor and acceptor genomes HGT events that can then be further investigated with tree-based approaches.

We're not picky about capitalization; choose whatever feels right to you (prehgt, preHGT) and feel free to capitalizate at the beginning of a sentence if it brings you joy (PreHGT).

For more scientific details about the preHGT, see [this pub](TODO: add doi link).

## Conceptual overview of the pipeline

This pipeline screens genomes for horizontal gene transfer.

The pipeline begins from a metadata table that records genome accessions and lineages.
The coding domain sequences (`_CDS_from_genomic.fna.gz`) and gene annotation files (`_gff.gz`) are downloaded and then combined (0.9 percent identity threshold) into a genus-level pseudo-pangenome (the clustering is very fast, but doesn't deal with ortholog/paralog problems -- the goal of this approach is not to identify a true pangenome, but to reduce the number of genes that are screened).
Each gene in the pangenome is then screened using composition metrics (codon usage, amino acid usage, etc.).
The results of the composition screen are then parsed to identify HGT candidates (TBD) and to remove likely contaminants (TBD).
The remaining HGT candidates are then BLASTED (DIAMOND) against the NCBI clustered non-redundant database and the results are parsed to identify genes with discrepant lineage matches (exact approach TBD).

## Quick Start

The pipeline is implemented as either a Snakemake or a Nextflow workflow, with environment management with conda.

### Nextflow

### Snakemake

## Additional documentation

### Outputs

The main output of the pipeline is a TSV file summarizing the 

Below we provide a description of each column output column.

* **hgt_candidate**: HGT candidate gene sequence name. The name is derived from the original FASTA CDS from genomic file downloaded from NCBI and contains the chromosome name, gene number, and GenBank protein ID if it exists. This is also the BLAST query sequence ID.
* **hgt_taxonomy_level**: the taxonomic level at which the HGT event was detected (kingdom, phylum, class, order, family)
* **blast_acceptor_lineage_at_hgt_taxonomy_level**: the taxonomic lineage of the acceptor genome up to the HGT level
* **blast_acceptor_num_matches_at_lineage**: number of BLAST hits at the acceptor_lineage_at_hgt_taxonomy_level
* **blast_acceptor_lca_level**: within the acceptor group, what level of taxonomy does the lowest common ancestor occur among all matches? If it’s at the phylum level, the HGT event is probably older than if it’s at the genus level. Or, if the HGT is only observed in two phyla, perhaps the HGT happened twice.
* **blast_acceptor_max_pident**: excluding self matches, the maximum percent identity of matches within the acceptor group.
* **blast_acceptor_max_bitscore**: excluding self matches, the maximum corrected bitscore of matches within the acceptor group.
* **blast_acceptor_min_evalue**: excluding self matches, the minimum evalue of matches within the acceptor group
* **blast_donor_lineage_at_hgt_taxonomy_level**: the lineage of the predicted donor group at the hgt_taxonomy_level
* **blast_donor_num_matches_at_lineage**: number of BLAST hits at the donor_lineage_at_hgt_taxonomy_level
* **blast_donor_best_match_full_lineage**: the full taxonomic lineage of the best match at hgt_taxonomy_level
* **blast_donor_best_match_id**: the sequence ide of the best match at hgt_taxonomy_level
* **blast_donor_best_match_pident**: the percent identity of the best match at hgt_taxonomy_level
* **blast_donor_max_bitscore**: the maximum corrected bitscore of the best match at hgt_taxonomy_level
* **blast_donor_min_evalue**: the minimum evalue of the best match at hgt_taxonomy_level
* **blast_alien_index**: a score of HGT probability based on e-value. > 0 possible HGT, >15 likely HGT, >45 highly likely HGT. (kingdom-level only)
* **blast_hgt_index**: a score of HGT probability based on bit score. Basically the same as alien index (kingdom-level only)
* **blast_donor_distribution_index**: within the donor groups, how frequently does the HGT candidate occur in all the donor groups? Is it only seen in bacteria, or is it in bacterial, viruses, and archaea? >0.8 means more specific to a donor group (kingdom-level only)
* **blast_entropy**: a measure of the uncertainty or randomness of a set of probabilities (kingdom-level only). Gives a higher number for a qseqid that is evenly distributed across all kingdoms. Gives a lower number for a qseqid that is mostly associated with a single kingdom.
* **blast_entropy_normalized**: the maximum value occurs when all kingdoms are equally represented. The entropy is -log2(1/K), where K is the number of kingdoms. Given 7 kingdoms, the max entropy would be log2(7). Entropy is normalized to 0-1 range by dividing by the maximum potential value.
* **blast_gini**: Gini Coefficient, a measure of inequality among values of a frequency distribution (kingdom-level only). A Gini coefficient close to 0 indicates that the qseqid is uniformly distributed across all kingdoms. A Gini coefficient close to 1 indicates that the qseqid is specific to one kingdom.
* **blast_acceptor_sum_bitscore_per_group_01**: the sum of all normalized (0-1) corrected bitscores within the acceptor_lineage_at_hgt_taxonomy_level
* **blast_donor_sum_bitscore_per_group_01**: the sum of all normalized (0-1) corrected bitscores within the donor_lineage_at_hgt_taxonomy_level
* **blast_ahs_01_index**: Aggregate hit score index, calculated by subtracting the sum of bitscores in the donor group from the sum of bitscores in the acceptor group when bitscores are 0-1 normalized. Modified from the AvP software.
* **blast_transfer_index**: Transfer index, calculated from bitscore ratios, taxonomic distances, and rank and total number of BLAST hits. Modified from the HGT-finder software.
* **blast_transfer_index_p_value**: p value for Transfer index
* **blast_transfer_index_adjusted_p_value**: FDR (p < 0.01) for transfer index p value
* **blast_total_num_matches**: the total number of BLAST hits returned. The maximum is 100, set earlier in the pipeline as the cutoff.

### Full usage

## Citations
