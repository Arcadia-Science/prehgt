process combine_hgt_candidates {
    tag "$genus"
    label 'process_single'

    conda "bioconda::csvtk=0.25.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.25.0--h9ee0642_0':
        'quay.io/biocontainers/csvtk:0.25.0--h9ee0642_0' }"

    input:
    tuple val(genus), path(input_pepstats_gene_lst)
    tuple val(genus), path(input_blastp_kingdom_gene_lst)
    tuple val(genus), path(input_blastp_subkingdom_gene_lst)
    
    output:
    tuple val(genus), path("*_combined_gene_lst.txt"), emit: combined_gene_lst

    script:
    def prefix = task.ext.prefix ?: "${genus}"
    """
    cat ${input_pepstats_gene_lst} ${input_blastp_kingdom_gene_lst} ${input_blastp_subkingdom_gene_lst} | csvtk freq -H -f 1 | csvtk cut -f 1 -o ${prefix}_combined_gene_lst.txt
    """
}
