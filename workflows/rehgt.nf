/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check input path parameters to see if they exist
def checkPathParamList = [params.input, params.blast_db, params.blast_db_tax, params.eggnog_db, params.eggnog_dmnd, params.hmm_db]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input)        { ch_input = file(params.input) } else { exit 1, 'TSV file specifying genera not provided!' }
if (params.blast_db)     { ch_BLAST_DB    = Channel.fromPath(params.blast_db) } else { exit 1, 'Path to blast database FASTA file not provided!' }
if (params.blast_db_tax) { ch_BLAST_TAX   = Channel.fromPath(params.blast_db_tax) } else { exit 1, 'Path to blast database taxonomy SQLITE file not provided!' }
if (params.eggnog_db)    { ch_EGGNOG_DB   = Channel.fromPath(params.eggnog_db) } else { exit 1, 'Path to eggnog database db file not provided!' }
if (params.eggnog_dmnd)  { ch_EGGNOG_DMND = Channel.fromPath(params.eggnog_dmnd) } else { exit 1, 'Path to eggnog database dmnd file not provided!' }
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
include { blastp_to_hgt_candidates        } from '../modules/blastp_to_hgt_candidates'
include { compositional_scans_pepstats    } from '../modules/compositional_scans_pepstats'
include { compositional_scans_to_hgt_candidates } from '../modules/compositional_scans_to_hgt_candidates'
include { combine_hgt_candidates          } from '../modules/combine_hgt_candidates'
include { extract_hgt_candidates          } from '../modules/extract_hgt_candidates'
include { combine_and_parse_gff_per_genus } from '../modules/combine_and_parse_gff_per_genus'
include { eggnog_hgt_candidates           } from '../modules/eggnog_hgt_candidates'
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

workflow REHGT {
    // Download genomes & GFF files
    download_reference_genomes(ch_GENUS) 
    // Parse GFF files
    combine_and_parse_gff_per_genus(download_reference_genomes.out.gff)
    // Build pangenomes
    build_genus_pangenome(download_reference_genomes.out.cds_all)
    translate_pangenome(build_genus_pangenome.out.rep_seq)
    // Do BLAST-based HGT candidate prediction
    blastp_against_clustered_nr(ch_BLAST_DB, translate_pangenome.out.aa_rep_seq)
    blastp_add_taxonomy_info(ch_BLAST_TAX, blastp_against_clustered_nr.out.tsv)
    blastp_to_hgt_candidates(blastp_add_taxonomy_info.out.tsv)
    // Do compositional HGT candidate prediction
    compositional_scans_pepstats(translate_pangenome.out.aa_rep_seq)
    compositional_scans_to_hgt_candidates(compositional_scans_pepstats.out.txt)
    // Annotate HGT candidates
    combine_hgt_candidates(compositional_scans_to_hgt_candidates.out.gene_lst, blastp_to_hgt_candidates.out.gene_lst)
    extract_hgt_candidates(translate_pangenome.out.aa_rep_seq, combine_hgt_candidates.out.gene_lst)
    eggnog_hgt_candidates(ch_EGGNOG_DB, ch_EGGNOG_DMND, extract_hgt_candidates.out.fasta)
    hmmscan_hgt_candidates(ch_HMM_DB, extract_hgt_candidates.out.fasta)
    // Combine all the results
    combine_results(compositional_scans_to_hgt_candidates.out.tsv,
                    blastp_to_hgt_candidates.out.blast_scores,
                    download_reference_genomes.out.csv,
                    build_genus_pangenome.out.cluster,
                    combine_and_parse_gff_per_genus.out.tsv,
                    eggnog_hgt_candidates.out.annotations,
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
