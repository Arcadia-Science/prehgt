process combine_results {
    tag "$genus"
    label 'process_low'

    conda "conda-forge::r-tidyverse=2.0.0 conda-forge::r-janitor=2.2.0"
    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/tidy-prehgt:2.0.0':
        '' }"

    input:
    tuple val(genus), path(compositional_tsv)
    tuple val(genus), path(blast_kingdom_tsv)
    tuple val(genus), path(blast_subkingdom_tsv)
    tuple val(genus), path(genomes_csv)
    tuple val(genus), path(pangenome_cluster_tsv)
    tuple val(genus), path(gff_tsv)
    tuple val(genus), path(kofamscan_tsv)
    tuple val(genus), path(hmmscan_tblout)

    output:
    tuple val(genus), path("*_results.tsv")     , emit: all_results
    tuple val(genus), path("*_method_tally.tsv"), emit: method_tally

    script:
    def prefix = task.ext.prefix ?: "${genus}"
    """
    combine_results.R ${compositional_tsv} ${blast_kingdom_tsv} ${blast_subkingdom_tsv} ${genomes_csv} ${pangenome_cluster_tsv} ${gff_tsv} ${kofamscan_tsv} ${hmmscan_tblout} ${prefix}_results.tsv ${prefix}_method_tally.tsv
    """
}
