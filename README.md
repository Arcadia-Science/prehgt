# preHGT: generating a *pre*liminary horizontal gene transfer (_HGT_) candidates using compositional and phylogenetic implicit approaches

[![Snakemake](https://img.shields.io/badge/snakemake--green)](https://snakemake.readthedocs.io/en/stable/)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/Arcadia-Science/prehgt)

preHGT is a pipeline that quickly pre-screens genomes across the tree of life for horizontal gene transfer (HGT).

The preHGT pipeline wraps existing parametric and implicit phylogenetic methods for HGT detection and reports multiple metrics about input genome contamination.
It quickly produces a "good-enough" candidate list of genes that researchers can further investigate with more stringent HGT detection methods, different data modalities, or wet lab experimentation.

For more scientific details about the preHGT, see [this pub](TODO: add doi link).

We're not picky about capitalization; choose whatever feels right to you (prehgt, preHGT) and feel free to capitalizate at the beginning of a sentence if it brings you joy (PreHGT).

## Quick Start

The pipeline is implemented as either a Snakemake or a Nextflow workflow, with environment management with conda.

### Installing conda

Software and environment management is performed with conda, regardless of if you run the pipeline with snakemake or nextflow.
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html) and for installing mamba [here](https://mamba.readthedocs.io/en/latest/).

### Obtaining databases

If you choose to run the workflow with Nextflow, you will provide the databases as command line parameters at run time.
If you choose to run the workflow with Snakemake, you will need the databases to be located in specific folders.
The download code below reflects the relative paths of the databases where the snakemake workflow will look for them in the `prehgt` workflow directory (see [Running the pipeline with Snakemake](#running-the-pipeline-with-snakemake) below for instructions on how to obtain the directory locally).

```
mkdir -p inputs/
curl -JLo inputs/nr_rep_seq.fasta.gz https://osf.io/gqxva/download # 59Gb BLASTp database
curl -JLo inputs/nr_cluster_taxid_formatted_final.sqlite https://osf.io/5qj7e/download # 66Gb BLASTp taxonomy database
curl -JLo inputs/hmms/all_hmms.hmm https://osf.io/f92qd/download # 3.2Gb HMM database
```

If you're using Nextflow, you will also need the KofamScan databases (these are automatically downloaded by the Snakemake pipeline)

```
mkdir -p inputs/kofamscandb
curl -JLo inputs/kofamscandb/ko_list.gz ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz && gunzip inputs/kofamscandb/ko_list.gz -C inputs/kofamscandb/
curl -JLo inputs/kofamscandb/profiles.tar.gz ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz && tar xf inputs/kofamscandb/profiles.tar.gz -C inputs/kofamscandb/
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

### Running the pipeline with [Nextflow Tower](https://cloud.tower.nf/)

[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/Arcadia-Science/2023-rehgt)

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
See the [database download instructions](#obtaining-databases) above to make sure you placed your databases in the correct locations.
Note, the KofamScan databases are downloaded by the pipeline itself, so no need to worry about those.

To start the pipeline, run:

```
snakemake --use-conda -j 8 --rerun-incomplete
```

where:

- `--use-conda` tells snakemake to manage software dependencies with conda.
- `-j` tells snakemake the maximum number of threads to use at once.
- `--rerun-incomplete` tells snakemake to check whether any files were interupted and to re-make them if so.

## Additional documentation

### Overview of the pipeline

The preHGT pipeline combines compositional scans, pangenome inference, and BLAST-based searches.
For the interpretation of the BLAST results, we took advantage of a rich literature of BLAST-based HGT predictor indices and re-implemented many with the goal of providing the most information possible from one BLAST run with a single tool.
To reduce overall run times, the pipeline employs clustering heuristics, including construction of a genus-level pangenome to reduce query size and searches against a clustered database to reduce database size.
To reduce false positives, we included multiple screens for contamination based on similarity to database matches, position of gene in the contiguous sequence, and homolog presence in closely related genomes.

Below, we provide an overview of what each step of the pipeline does.

1. **Retrieving gene sequences and annotation files.** The pipeline begins with the user providing a genus or genera of interest in a TSV file. The pipeline then scans GenBank and RefSeq for matching genomes and downloads relevant files. When a genome is available in both GenBank and RefSeq, only the RefSeq version is retained. This step also parses the input files in preparation for future steps.
   - **[ncbi-genome-download](https://github.com/kblin/ncbi-genome-download):** Loops through GenBank and RefSeq to find and download all genomes with gene models (`*_cds_from_genomic.fna.gz`) and genome annotation files (`*_genomic.gff.gz`).
   - **[delete_gca_files.sh](./bin/delete_gca_files.sh)**: If a genome is in both GenBank and RefSeq, this script deletes the GenBank version and only keep the RefSeq version.
   - **file parsing**: All `*_cds_from_genomic.fna.gz` files are combined into a single file, and a CSV file that records the total number of downloaded genomes per genus is generated.
2. **Building a pangenome.** For each genus, the pipeline then combines genes into a pseudopangenome, which reduces the number of genes that are investigated and provides metadata about the gene.
   - **[`mmseqs easy-cluster`](https://github.com/soedinglab/MMseqs2)**: Nucleotide sequences are clustered at 90% length and identity.
   - **[EMBOSS `transeq`](https://emboss.sourceforge.net/apps/cvs/emboss/apps/transeq.html)**: clustered nucleotide sequences are translated into amino acid sequences.
3. **Detecting HGT candidates.** Using the genes in the pangenome, the pipeline uses two categories of approaches to identify HGT candidates.
   - **Compositional scan.** The first approach uses relative amino acid usage to detect proteins with outlying composition.
     - [**EMBOSS `pepstats`**](https://embossgui.sourceforge.net/demo/manual/pepstats.html): Measures relative amino acid usage (RAAU) for each gene.
     - **[compositional_scans_to_hgt_candidates.R](./bin/compositional_scans_to_hgt_candidates.R)**: Produces a list of genes that have outlying RAAU. Starts by parsing the `pepstats` results with the function `parse_pepstats_to_amino_acid_frequencies()`. Then, produces a distance matrix with the base R function `dist()` and hierarchically clusters the distance matrix with fastcluster’s `hclust()`. It detects outliers by cutting the resultant tree with `height/1.5` and retaining any cluster that contains fewer than 0.1% of the pangenome size.
   - **BLASTp scan.** Uses BLASTp to identify proteins with distant homologs.
     - **[DIAMOND `blastp`](https://github.com/bbuchfink/diamond)**: All genes in the pseudopangenome are BLASTed against a [clustered version of NCBI’s nr database (90% length, 90% identity)](https://github.com/Arcadia-Science/2023-nr-clustering). The clustered database makes the BLASTp step faster and ensures results contain a variety of taxonomic lineages in cases where distant homology exists.
     - **[blastp_add_taxonomy_info.R](./bin/blastp_add_taxonomy_info.R)**: Adds lineage information to the BLASTp search using dplyr, dbplyr, and RSQLite.
     - **[blastp_to_hgt_candidates_kingdom.R](./bin/blastp_to_hgt_candidates_kingdom.R)**:
       - [Alien index](https://doi.org/10.1126/science.1156407): A score of HGT probability based on BLAST hit e-value for the top hit in the kingdom-level acceptor versus donor hits. Because it is based on e-value, it can be biased by BLAST database size.
       - [HGT score](https://doi.org/10.1371/journal.pgen.1003035): A score of HGT probability based on BLAST hit corrected bitscore for the top hit in the kingdom-level acceptor versus donor hits.
       - [Donor distribution index](https://doi.org/10.1016/j.molp.2022.02.001): An index that calculates the specificity of a HGT candidate gene within the donor groups. It estimates how frequently does the HGT candidate occurs in all the donor groups.
       - (Normalized) entropy: A measure of the uncertainty or randomness of a set of probabilities. Since there are seven kingdoms investigated, the maximum entropy is `log2(7)`. Entropy is normalized to 0-1 range by dividing by the maximum potential value.
       - Gini coefficient: A measure of inequality among values of a frequency distribution. A Gini coefficient close to 0 indicates that the qseqid is uniformly distributed across all kingdoms.
       - [Aggregate hit score](https://doi.org/10.1371/journal.pcbi.1010686): A score of HGT liklihood calculated by subtracting the sum of normalized corrected bitscores in the donor group from the sum of normalized corrected bitscores in the acceptor group.
     - **[blastp_to_hgt_candidates_subkingdom.R](./bin/blastp_to_hgt_candidates_subkingdom.R)**:
       - [Transfer index](https://doi.org/10.3390/toxins7104035): A score for HGT probability calculated by comparing the corrected bitscores and taxonomic distance of all BLAST hits against the best BLAST hit to the database.
4. **Annotation.** The pipeline then annotates the HGT candidates.
   - **[combine_and_parse_gff_per_genus.R](./bin/combine_and_parse_gff_per_genus.R)**: All downloaded `*_genomic.gff.gz` are combined and parsed to pull out annotation information for each coding domain sequence.
   - **[KofamScan](https://github.com/takaram/kofam_scan)**: KofamScan uses hidden markov models (HMMs) to perform KEGG ortholog annotation.
   - **[hmmscan](http://hmmer.org/)**: We [built](./make_hmm_db.snakefile) a [custom HMM database](https://osf.io/trgpc/) to scan for annotations of interest using HMMER3 `hmmscan`. The HMM database currently contains VOGs from [VOGDB](https://vogdb.org/) and [biosynthetic genes](./inputs/hmms/hmm_urls.csv), and can be extended in the future to meet user annotation interests.
5. **Reporting.** The last step combines all information that the pipeline has produced and outputs the results in a TSV file using the script [combine_results.R](./bin/combine_results.R). The results include the GenBank protein identifier for the HGT candidate, BLAST and relative amino acid usage scores, pangenome information, gene and ortholog annotations, and contextualizing information about the gene such as position in the contiguous sequence. This script also reports whether an HGT candidate is actually contamination by integrating information about the number of times the gene is observed in the pangenome, the similarity of the gene to other hits in the database, and the length of the contiguous sequence that the gene is located on in the genome.

TODO: add pipeline overview DAG/figure.

### Outputs

The main output is a TSV file summarizing the results of each part of the pipeline. Below we provide a description of each column output column.

#### TSV summary file

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

#### Other outputs

The pipeline also produces a number of intermediate files, including the original `*_cds_from_genomic.fna.gz` files, full BLAST results, the HGT candidate gene sequences, and the pangenome gene cluster membership manifest.
The exact path of these files will depend on whether the pipeline is executed in Snakemake or Nextflow.

## Citations
