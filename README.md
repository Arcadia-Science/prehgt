# preHGT: generating a *pre*liminary horizontal gene transfer (_HGT_) candidates using compositional and phylogenetic implicit approaches

[![Snakemake](https://img.shields.io/badge/snakemake--green)](https://snakemake.readthedocs.io/en/stable/)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/Arcadia-Science/prehgt)

preHGT a pipeline quickly pre-screens genomes across the tree of life for horizontal gene transfer (HGT).

The preHGT pipeline wraps existing parametric and implicit phylogenetic methods for HGT detection and reports multiple metrics about input genome contamination.
It quickly produces a "good-enough" candidate list of genes that researchers can further investigate with more stringent HGT detection methods, different data modalities, or wet lab experimentation.

For more scientific details about the preHGT, see [this pub](TODO: add doi link).

We're not picky about capitalization; choose whatever feels right to you (prehgt, preHGT) and feel free to capitalizate at the beginning of a sentence if it brings you joy (PreHGT).

## Overview of the pipeline

The preHGT pipeline combines compositional scans, pangenome inference, and BLAST-based searches.
For the interpretation of the BLAST results, we took advantage of a rich literature of BLAST-based HGT predictor indices and re-implemented many with the goal of providing the most information possible from one BLAST run with a single tool.
To reduce overall run times, the pipeline employs clustering heuristics, including construction of a genus-level pangenome to reduce query size and searches against a clustered database to reduce database size.
To reduce false positives, we included multiple screens for contamination based on similarity to database matches, position of gene in the contiguous sequence, and homolog presence in closely related genomes.

TODO: add pipeline overview DAG/figure.

## Quick Start

The pipeline is implemented as either a Snakemake or a Nextflow workflow, with environment management with conda.

### Installing conda

Software and environment management is performed with conda, regardless of if you run the pipeline with snakemake or nextflow.
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html) and for installing mamba [here](https://mamba.readthedocs.io/en/latest/).

### Obtaining databases

If you choose to run the workflow with Nextflow, you will provide the databases as command line parameters at run time.
If you choose to run the workflow with Snakemake, you will need the databases to be located in specific folders.
The download code below reflects the relative paths of the databases where the snakemake workflow will look for them in the `prehgt` workflow directory (see [Running the pipeline with Snakemake](###running-the-pipeline-with-snakemake) below for instructions on how to obtain the directory locally).

```
curl -JLo
curl -JLo
curl -JLo
```

### Input sample sheet

The input sample sheet is the same for both workflows.
It's a TSV file with the genus or genera you would like to run through the preHGT pipeline.
An example is included [here](https://github.com/Arcadia-Science/test-datasets/blob/main/rehgt/bigelowiella_test.tsv) and reproduced below:

```
genus
Bigelowiella
Schizosaccharomyces
```

### Running the pipeline with Nextflow

If you would like to run the pipeline with Nextflow, begin by installing Nextflow.

```
mamba env create -n prehgtnf nextflow=22.10.6
conda activate prehgtnf
```

Then, you can test the workflow with a small data set using:

```
nextflow run Arcadia-Science/prehgt -profile test,conda --outdir <OUTDIR>
```

To run the full pipeline, including supplying paths to databases, run:

```
nextflow run Arcadia-Science/prehgt -profile conda --outdir <OUTDIR> --input input.tsv --blast_db path_to_blast_db --blast_db_tax path_to_blast_db_taxonomy_file --ko_list path_to_ko_list_file --ko_profiles path_to_ko_profiles_folder --hmm_db path_to_hmm_db
```

### Running the pipeline with Snakemake

If you would like to run the pipeline with Snakemake, begin by installing Snakemake and other runtime dependencies.

```
mamba env create -n prehgt --file environment.yml
conda activate prehgt
```

Next, clone this repository and `cd` into it:

```
git clone https://github.com/Arcadia-Science/prehgt.git
cd prehgt
```

The input files for the snakemake workflow are not parameterized, so you need to make sure the input sample sheet and databases are in the correct location with the correct file names.
We recap these below.
Note, the KofamScan databases are downloaded by the pipeline itself, so no need to worry about those.

```
TO LIST OUT PATHS
```

To start the pipeline, run:

```
snakemake --use-conda -j 8 --rerun-incomplete
```

where:

- `--use-conda` tells snakemake to manage software dependencies with conda.
- `-j` tells snakemake the maximum number of threads to use at once.
- `--rerun-incomplete` tells snakemake to check whether any files were interupted and to re-make them if so.

## Additional documentation

### Outputs

The main output is a TSV file summarizing the results of each part of the pipeline. Below we provide a description of each column output column.

<summary> <b>Overview columns</b> </summary>

- **genus**: the input genus for a given HGT candidate
- **hgt_candidate**: HGT candidate gene sequence name. The name is derived from the original FASTA CDS from genomic file downloaded from NCBI and contains the chromosome name, gene number, and GenBank protein ID if it exists. This is also the BLAST query sequence ID.
- **method**: the method used to infer the HGT candidate. Currently blast, raau (relative amino acid usage), or both
<details>
<summary> <b>Relative amino acid usage columns</b> </summary>

- **RAAU_cluster**: the cluster that an HGT candidate was in. If multiple HGT candidates are reported for the RAAU method, the cluster might indicate if genes have a similar RAAU (which may indicate shared donor or evolutionary history).
</details>
<details>
<summary> <b>BLAST columns</b> </summary>

- **blast_algorithm_type**: the BLAST algorithm type used to infer the HGT candidate. One of eith kingdom or sub-kingdom.
- **blast_HGT_score**: HGT score inferred from Alien Index. Also reports contamination liklihood. Since Alien index is only calculated for kingdom level transfers, the score will be NA for all sub-kingdom algorithm type results.
- **blast_hgt_taxonomy_level**: the taxonomic level at which the HGT event was detected (kingdom, phylum, class, order, family)
- **blast_acceptor_lineage_at_hgt_taxonomy_level**: the taxonomic lineage of the acceptor genome up to the HGT level
- **blast_acceptor_lca_level**: within the acceptor group, what level of taxonomy does the lowest common ancestor occur among all matches? If it’s at the phylum level, the HGT event is probably older than if it’s at the genus level. Or, if the HGT is only observed in two phyla, perhaps the HGT happened twice.
- **blast_acceptor_max_pident**: excluding self matches, the maximum percent identity of matches within the acceptor group.
- **blast_acceptor_max_bitscore**: excluding self matches, the maximum corrected bitscore of matches within the acceptor group.
- **blast_acceptor_min_evalue**: excluding self matches, the minimum evalue of matches within the acceptor group
- **blast_acceptor_num_matches_at_lineage**: number of BLAST hits at the acceptor_lineage_at_hgt_taxonomy_level
- **blast_donor_num_matches_at_lineage**: number of BLAST hits at the donor_lineage_at_hgt_taxonomy_level
- **blast_total_num_matches**: The total number of BLAST hits returned. The maximum is 100, set earlier in the pipeline as the cutoff.
- **blast_donor_lineage_at_hgt_taxonomy_level**: the lineage of the predicted donor group at the hgt_taxonomy_level
- **blast_donor_best_match_full_lineage**: the full taxonomic lineage of the best match at hgt_taxonomy_level
- **blast_donor_best_match_id**: the sequence ide of the best match at hgt_taxonomy_level
- **blast_donor_best_match_pident**: the percent identity of the best match at hgt_taxonomy_level
- **blast_donor_max_bitscore**: the maximum corrected bitscore of the best match at hgt_taxonomy_level
- **blast_donor_min_evalue**: the minimum evalue of the best match at hgt_taxonomy_level
- **blast_alien_index**: a score of HGT probability based on e-value. > 0 possible HGT, >15 likely HGT, >45 highly likely HGT. (kingdom-level only)
- **blast_hgt_index**: a score of HGT probability based on bit score. Basically the same as alien index (kingdom-level only)
- **blast_donor_distribution_index**: within the donor groups, how frequently does the HGT candidate occur in all the donor groups? Is it only seen in bacteria, or is it in bacterial, viruses, and archaea? >0.8 means more specific to a donor group (kingdom-level only)
- **blast_entropy**: a measure of the uncertainty or randomness of a set of probabilities (kingdom-level only). Gives a higher number for a qseqid that is evenly distributed across all kingdoms. Gives a lower number for a qseqid that is mostly associated with a single kingdom.
- **blast_entropy_normalized**: the maximum value occurs when all kingdoms are equally represented. The entropy is -log2(1/K), where K is the number of kingdoms. Given 7 kingdoms, the max entropy would be log2(7). Entropy is normalized to 0-1 range by dividing by the maximum potential value.
- **blast_gini**: Gini Coefficient, a measure of inequality among values of a frequency distribution (kingdom-level only). A Gini coefficient close to 0 indicates that the qseqid is uniformly distributed across all kingdoms. A Gini coefficient close to 1 indicates that the qseqid is specific to one kingdom.
- **blast_acceptor_sum_bitscore_per_group_01**: the sum of all normalized (0-1) corrected bitscores within the acceptor_lineage_at_hgt_taxonomy_level
- **blast_donor_sum_bitscore_per_group_01**: the sum of all normalized (0-1) corrected bitscores within the donor_lineage_at_hgt_taxonomy_level
- **blast_ahs_01_index**: aggregate hit score index, calculated by subtracting the sum of bitscores in the donor group from the sum of bitscores in the acceptor group when bitscores are 0-1 normalized. Modified from the AvP software.
- **blast_transfer_index**: transfer index, calculated from bitscore ratios, taxonomic distances, and rank and total number of BLAST hits. Modified from the HGT-finder software.
- **blast_transfer_index_p_value**: p value for Transfer index
- **blast_transfer_index_adjusted_p_value**: FDR (p < 0.01) for transfer index p value
</details>
<details>
<summary> <b>Annotation columns</b> </summary>

- **kofamscan_ko**: best KEGG ortholog identifier reported by KofamScan annotation
- **kofamscan_threshold**: KofamScan HMM threshold; family-specific adaptive score calculated by KofamScan for each KO family.
- **kofamscan_score**: KofamScan HMM score. We do not filter annotations that do not meet the threshold, as these may still be useful for sleuthing out potential functions of an HGT candidate. |
- **kofamscan_evalue**: KofamScan HMM evalue
- **kofamscan_ko_definition**: Full definition for best KEGG ortholog
- **hmmscan_domain_name**: domain name for hmmscan hit
- **hmmscan_description**: description of domain for hmmscan hit
- **hmmscan_sequence_evalue**: evalue for hmmscan hit
- **hmmscan_sequence_score**: sequence score for hmmscan hit
- **hmmscan_sequence_bias**: sequence bias for hmmscan hit
- **hmmscan_best_domain_evalue**: best domain evalue for hmmscan hit
- **hmmscan_best_domain_score**: best domain score for hmmscan hit
</details>
<details>
<summary> <b>Pangenome columns</b> </summary>

- **pangenome_num_genes_in_cluster**: total number of genes in the HGT candidate's cluster. Each HGT candidate is the representative sequence for its cluster.
- **pangenome_size**: total number of genomes in the genus-level pangenome.
</details>
<details>
<summary> <b>GFF columns</b> </summary>

- **gff_seqid**: sequence ID in the GFF file
- **gff_source**: source in the GFF file
- **gff_feature**: feature annotation (CDS) in the GFF file
- **gff_start**: sequence start in the GFF file
- **gff_end**: sequence end in the GFF file
- **gff_score**: score in the GFF file
- **gff_strand**: strand in the GFF file
- **gff_frame**: frame in the GFF file. If a CDS is reported in multiple frames, only one is retained and the total number of frames observed in the original GFF file is reported in the column gff_frame_tally.
- **gff_attribute**: GFF attribute
- **gff_seqid_length**: Sequence length for a given seqid in a GFF file.
- **gff_Dbxref**: GFF Dbxref
- **gff_gbkey**: GFF GenBank key
- **gff_gene**: GFF gene name
- **gff_ID**: GFF ID
- **gff_locus_tag**: GFF locus tag
- **gff_Name**: GFF name
- **gff_Note**: GFF note
- **gff_Parent**: GFF parent ID
- **gff_product**: GFF protein product
- **gff_protein_id**: GFF protein ID
- **gff_transl_table**: GFF translation table
- **gff_frame_tally**: total number of frames in which a CDS was observed in the original GFF file. One information for the first frame is retained to reduce reporting redundancy.
</details>

### Full usage

TBD

## Citations
