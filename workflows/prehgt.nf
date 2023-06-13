/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check input path parameters to see if they exist
def checkPathParamList = [params.input, params.blast_db, params.blast_db_tax, params.ko_list, params.ko_profiles, params.hmm_db]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input)        { ch_input = file(params.input) } else { exit 1, 'TSV file specifying genera not provided!' }
if (params.blast_db)     { ch_BLAST_DB    = Channel.fromPath(params.blast_db) } else { exit 1, 'Path to blast database FASTA file not provided!' }
if (params.blast_db_tax) { ch_BLAST_TAX   = Channel.fromPath(params.blast_db_tax) } else { exit 1, 'Path to blast database taxonomy SQLITE file not provided!' }
if (params.ko_list)      { ch_KO_LIST     = Channel.fromPath(params.ko_list) } else { exit 1, 'Path to ko_list file not provided!' }
if (params.ko_profiles)  { ch_KO_PROFILES = Channel.fromPath(params.ko_profiles) } else { exit 1, 'Path to ko profiles archive not provided!' }
if (params.hmm_db)       { ch_HMM_DB      = Channel.fromPath(params.hmm_db) } else { exit 1, 'Path to hmm database not provided!' }

// Parse input file to retrieve genera to run pipeline on
metadata = ch_input
    .readLines()
    .drop(1)
    .collect { it.split("\t") }

ch_GENUS = Channel.from(metadata.collect { it[0] }).unique()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { download_reference_genomes      } from '../modules/download_reference_genomes'
include { build_genus_pangenome           } from '../modules/build_genus_pangenome'
include { translate_pangenome             } from '../modules/translate_pangenome'
include { blastp_against_clustered_nr     } from '../modules/blastp_against_clustered_nr'
include { blastp_add_taxonomy_info        } from '../modules/blastp_add_taxonomy_info'
include { blastp_to_hgt_candidates_subkingdom } from '../modules/blastp_to_hgt_candidates_subkingdom'
include { blastp_to_hgt_candidates_kingdom} from '../modules/blastp_to_hgt_candidates_kingdom'
include { compositional_scans_pepstats    } from '../modules/compositional_scans_pepstats'
include { compositional_scans_to_hgt_candidates } from '../modules/compositional_scans_to_hgt_candidates'
include { combine_hgt_candidates          } from '../modules/combine_hgt_candidates'
include { extract_hgt_candidates          } from '../modules/extract_hgt_candidates'
include { combine_and_parse_gff_per_genus } from '../modules/combine_and_parse_gff_per_genus'
include { kofamscan_hgt_candidates        } from '../modules/kofamscan_hgt_candidates'
include { hmmscan_hgt_candidates          } from '../modules/hmmscan_hgt_candidates'
include { combine_results                 } from '../modules/combine_results'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// MODULE: Installed directly from nf-core/modules
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREHGT {
    /* DOWNLOAD GENOMES & GFF FILES AND PARSE */

    // Using user-provided genera names, this process uses the ncbi-genome-download tool to download the GFF annotation file
    // and coding domain sequence (CDS) file for all member of that genus in GenBank and RefSeq.
    // It then deletes any genomes that are in RefSeq and GenBank, keeping only the RefSeq version.
    // Next, it combines the raw genome files for that genus into a single FASTA file.
    // Lastly, it reports the genomes that were downloaded into a CSV file.
    download_reference_genomes(ch_GENUS)

    // This process combines all GFF annotation files from all accessions that belong to a given genus into one file.
    // After pangenome clustering, this file will be used to retrieve information (genomic coords, etc) for rep sequences.
    combine_and_parse_gff_per_genus(download_reference_genomes.out.gff)

    /* BUILD PANGENOME */

    // This process clusters all CDS sequences from all accessions in a given genus into a pangenome.
    // This reduces redundancy for subsequent searches, and this information can be leveraged to interpret the HGT candidate genes
    // that will eventually be predicted by the pipeline.
    // selecting clustering threshold:
    // 0.9 https://www.science.org/doi/full/10.1126/sciadv.aba0111
    // 0.98, 0.99 https://www.sciencedirect.com/science/article/pii/S0960982220314263
    // 0.9 https://www.pnas.org/doi/abs/10.1073/pnas.2009974118
    build_genus_pangenome(download_reference_genomes.out.cds_all)

    // Translate the representative pangenome sequences from nucleotide to amino acid
    translate_pangenome(build_genus_pangenome.out.rep_seq)

    /* BLAST-BASED HGT CANDIDATE PREDICTION */

    // This process uses the diamond implementation of BLASTP to compare each CDS in the genus-level pangenome to a clustered version of NR.
    // For more information on the database, see this repository: https://github.com/Arcadia-Science/2023-nr-clustering
    // Using the diamond implementation and the clustered database decreases the time it takes to run this step.
    blastp_against_clustered_nr(ch_BLAST_DB, translate_pangenome.out.aa_rep_seq)

    // This process adds taxonomic lineages to each BLAST match.
    // The taxonomy sheet records the lowest common ancestor for a BLAST match, given that all BLAST matches represent a cluster of proteins.
    // Because the taxonomy sheet is so large (>60GB), the script uses an sql query executed via dplyr and dbplyr to decrease search times.
    blastp_add_taxonomy_info(ch_BLAST_TAX, blastp_against_clustered_nr.out.tsv)

    // This script processes BLAST matches and their taxonomic lineages to identify kingdom-level HGT candidates using alien index, horizontal gene transfer index,
    // donor distribution index, and acceptor lowest common ancestor calculations.
    // It scores all candidates and writes the scores and other relevant information to a TSV file and outputs a list of candidate gene IDs.
    blastp_to_hgt_candidates_kingdom(blastp_add_taxonomy_info.out.tsv)

    // This script processes BLAST matches and their taxonomic lineages to identify sub-kingdom HGT candidates using transfer index.
    // It scores all candidates and writes the scores and other relevant information to a TSV file and outputs a list of candidate gene IDs.
    blastp_to_hgt_candidates_subkingdom(blastp_add_taxonomy_info.out.tsv)

    /* COMPOSITION HGT CANDIDATE PREDICTION */

    // This process uses emboss pepstats to calculate relative amino acid usage per pangenome seq
    compositional_scans_pepstats(translate_pangenome.out.aa_rep_seq)

    // This script performs hierarchical clustering on genes based on their relative amino acid usage and identifies clusters with fewer than 150 genes.
    // The output is a list of genes that belong to these small clusters, which are considered HGT candidates.
    compositional_scans_to_hgt_candidates(compositional_scans_pepstats.out.txt)

    /* ANNOTATE HGT CANDIDATES */

    // This process combines the lists of HGT candidate genes identified through two different methods - BLAST and compositional scans.
    // The output is a single list containing the unique genes from both input lists.
    combine_hgt_candidates(compositional_scans_to_hgt_candidates.out.gene_lst, 
                           blastp_to_hgt_candidates_kingdom.out.gene_lst,
                           blastp_to_hgt_candidates_subkingdom.out.gene_lst)

    // This process extracts the amino acid sequences of the HGT candidate genes from the pangenome clusters.
    // The outputs is a FASTA file containing the sequences of the identified HGT candidate genes.
    extract_hgt_candidates(translate_pangenome.out.aa_rep_seq, combine_hgt_candidates.out.gene_lst)

    // This process performs KEGG ortholog annotation via hmm searches using the tool kofamscan
    kofamscan_hgt_candidates(ch_KO_LIST, ch_KO_PROFILES, extract_hgt_candidates.out.fasta)

    // Using the hmm file built in make_hmm_db.snakefile, this process uses hidden markov models to annotate specific protein classes of interest.
    // At the moment it targets viruses and biosynthetic gene clusters, but the HMM file can be expanded in the aforementioned snakefile as desired.
    hmmscan_hgt_candidates(ch_HMM_DB, extract_hgt_candidates.out.fasta)

    /* COMBINE ALL RESULTS */

    // Combine all of the results into a single mega TSV file.
    // The results are joined either on the genus or on the HGT candidate gene name, derived from the pangenome FASTA file.
    combine_results(compositional_scans_to_hgt_candidates.out.tsv,
                    blastp_to_hgt_candidates_kingdom.out.blast_scores,
                    blastp_to_hgt_candidates_subkingdom.out.blast_scores,
                    download_reference_genomes.out.csv,
                    build_genus_pangenome.out.cluster,
                    combine_and_parse_gff_per_genus.out.tsv,
                    kofamscan_hgt_candidates.out.tsv,
                    hmmscan_hgt_candidates.out.tblout)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
