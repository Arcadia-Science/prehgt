process blastp_to_hgt_candidates_kingdom {
    tag "$genus"
    label 'process_low'

    conda "conda-forge::r-tidyverse=2.0.0"
    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/prehgt-tidyverse:2.0.0':
        '' }"
    // keep this in case tidyverse=2.0.0 ever becomes available as a container by default
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/diamond:2.1.6--h5b5514e_1':
    //    'quay.io/biocontainers/diamond:2.1.6--h5b5514e_1' }"

    input:
    tuple val(genus), path(blast_lineages_tsv)

    output:
    tuple val(genus), path("*_blastp_kingdom_scores.tsv")  , emit: blast_scores
    tuple val(genus), path("*_blastp_kingdom_gene_lst.txt"), emit: gene_lst

    script:
    def prefix = task.ext.prefix ?: "${genus}"
    """
    blastp_to_hgt_candidates_kingdom.R ${blast_lineages_tsv} ${prefix}_blastp_kingdom_scores.tsv ${prefix}_blastp_kingdom_gene_lst.txt
    """
}
