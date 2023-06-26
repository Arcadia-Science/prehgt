process combine_and_parse_gff_per_genus {
    tag "$genus"
    label 'process_low'

    conda "conda-forge::r-tidyverse=2.0.0"
    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/tidy-prehgt:2.0.0':
        '' }"

    input:
    tuple val(genus), path(input_gffs)

    output:
    tuple val(genus), path("*_gff_info.tsv"), emit: tsv

    script:
    def prefix = task.ext.prefix ?: "${genus}"
    """
    combine_and_parse_gff_per_genus.R ${prefix}_gff_info.tsv ${input_gffs}
    """
}
