process blastp_add_taxonomy_info {
    tag "$genus"
    label 'process_low'

    conda "conda-forge::r-tidyverse=2.0.0 conda-forge::r-dbplyr=2.3.2 conda-forge::r-dbi=1.1.3 conda-forge::r-rsqlite=2.3.1"
    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/tidy-prehgt:2.0.0':
        '' }"

    input:
    path(input_sqlite_db)
    tuple val(genus), path(blast_tsv)

    output:
    tuple val(genus), path("*_vs_clustered_nr_lineages.tsv"), emit: tsv

    script:
    def prefix = task.ext.prefix ?: "${genus}"
    """
    blastp_add_taxonomy_info.R ${input_sqlite_db} ${blast_tsv} ${prefix}_vs_clustered_nr_lineages.tsv 
    """
}
