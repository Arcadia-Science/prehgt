process translate_pangenome {
    tag "$genus"
    label 'process_single'

    conda "bioconda::emboss=6.6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/emboss:6.6.0--h86d058a_5':
        'quay.io/biocontainers/emboss:6.6.0--h86d058a_5' }"

    input:
    tuple val(genus), path(input_rep_seq)

    output:
    tuple val(genus), path("*aa_rep_seq.fasta"), emit: aa_rep_seq

    script:
    def prefix = task.ext.prefix ?: "${genus}"
    """
    transeq -sequence ${input_rep_seq} -outseq ${prefix}_aa_rep_seq.fasta
    """
}
