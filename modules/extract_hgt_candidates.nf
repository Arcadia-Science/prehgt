process extract_hgt_candidates {
    tag "$genus"
    label 'process_single'

    conda "bioconda::seqtk=1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.3--h7132678_4':
        'quay.io/biocontainers/seqtk:1.3--h7132678_4' }"

    input:
    tuple val(genus), path(input_aa_fasta), path(input_gene_lst)
    
    output:
    tuple val(genus), path("*_aa.fasta")    , emit: fasta

    script:
    def prefix = task.ext.prefix ?: "${genus}"
    """
    seqtk subseq ${input_aa_fasta} ${input_gene_lst} > ${prefix}_aa.fasta
    """
}
