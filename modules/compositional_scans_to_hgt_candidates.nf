process compositional_scans_to_hgt_candidates {
    tag "$genus"
    label 'process_high_memory'

    conda "conda-forge::r-tidyverse=2.0.0 conda-forge::r-fastcluster=1.2.3"
    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/tidy-prehgt:2.0.0':
        '' }"

    input:
    tuple val(genus), path(input_pepstats_txt)

    output:
    tuple val(genus), path("*_clusters.tsv")         , emit: tsv
    tuple val(genus), path("*_pepstats_gene_lst.txt"), emit: gene_lst

    script:
    def prefix = task.ext.prefix ?: "${genus}"
    """
    compositional_scans_to_hgt_candidates.R ${input_pepstats_txt} ${prefix}_clusters.tsv ${prefix}_pepstats_gene_lst.txt
    """
}
