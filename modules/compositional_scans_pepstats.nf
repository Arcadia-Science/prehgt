process compositional_scans_pepstats {
    tag "$genus"
    label 'process_low'

    conda "bioconda::emboss=6.6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/emboss:6.6.0--hf657eab_5':
        'quay.io/biocontainers/emboss:6.6.0--h440b012_4' }"

    input:
    tuple val(genus), path(input_aa_rep_seq)

    output:
    tuple val(genus), path("*_pepstats.txt"), emit: txt

    script:
    def prefix = task.ext.prefix ?: "${genus}"
    """
    pepstats -sequence ${input_aa_rep_seq} -outfile ${prefix}_pepstats.txt
    """
}
