process blastp_to_hgt_candidates_subkingdom {
    tag "$genus"
    label 'process_low'

    conda "conda-forge::r-tidyverse=2.0.0"
    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/tidy-prehgt:2.0.0':
        '' }"

    input:
    tuple val(genus), path(blast_lineages_tsv)
    val padj_threshold

    output:
    tuple val(genus), path("*_blastp_subkingdom_scores.tsv")  , emit: subkingdom_blast_scores
    tuple val(genus), path("*_blastp_subkingdom_gene_lst.txt"), emit: subkingdom_gene_lst

    script:
    def prefix = task.ext.prefix ?: "${genus}"
    """
    blastp_to_hgt_candidates_subkingdom.R ${blast_lineages_tsv} ${padj_threshold} ${prefix}_blastp_subkingdom_scores.tsv ${prefix}_blastp_subkingdom_gene_lst.txt
    """
}
