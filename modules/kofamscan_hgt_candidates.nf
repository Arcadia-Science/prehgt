process kofamscan_hgt_candidates {
    tag "$genus"
    label 'process_blast'

    conda "bioconda::kofamscan=1.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kofamscan:1.3.0--hdfd78af_2':
        'quay.io/biocontainers/kofamscan:1.3.0--hdfd78af_2' }"

    input:
    path(ko_list)
    path(ko_profiles)
    tuple val(genus), path(input_aa_fasta)
    
    output:
    tuple val(genus), path("*_kofamscan.tsv"), emit: kofamscan_tsv

    script:
    def prefix = task.ext.prefix ?: "${genus}"
    """
    mkdir -p tmp
    gunzip -c ${ko_list} > ko_list
    tar xf ${ko_profiles}
    exec_annotation --format detail-tsv --ko-list ko_list --profile profiles --cpu $task.cpus -o ${prefix}_kofamscan.tsv ${input_aa_fasta}
    """
}
