process blastp_to_hgt_candidates_kingdom {
    tag "$genus"
    label 'process_low'

    conda "conda-forge::r-tidyverse=2.0.0"
    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/tidy-prehgt:2.0.0':
        '' }"

    input:
    tuple val(genus), path(blast_lineages_tsv)

    output:
    tuple val(genus), path("*_blastp_kingdom_scores.tsv")  , emit: kingdom_blast_scores
    tuple val(genus), path("*_blastp_kingdom_gene_lst.txt"), emit: kingdom_gene_lst

    script:
    def prefix = task.ext.prefix ?: "${genus}"
    """
    blastp_to_hgt_candidates_kingdom.R ${blast_lineages_tsv} ${prefix}_blastp_kingdom_scores.tsv ${prefix}_blastp_kingdom_gene_lst.txt
    """
}
