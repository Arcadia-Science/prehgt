process hmmscan_hgt_candidates {
    tag "$genus"
    label 'process_blast'

    conda "bioconda::hmmer=3.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.3.2--h87f3376_2':
        'quay.io/biocontainers/hmmer:3.3.2--h1b792b2_1' }"

    input:
    path(db)
    tuple val(genus), path(input_aa_fasta)
    
    output:
    tuple val(genus), path("*.tblout"), emit: tblout
    tuple val(genus), path("*.out")   , emit: out

    script:
    def prefix = task.ext.prefix ?: "${genus}"
    """
    hmmpress ${db}
    hmmscan --cpu $task.cpus --tblout ${prefix}.tblout -o ${prefix}.out ${db} ${input_aa_fasta}
    """
}
