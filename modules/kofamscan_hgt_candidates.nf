process kofamscan_hgt_candidates {
    tag "$genus"
    label 'process_blast'

    conda "bioconda::kofamscan=1.3.0 conda-forge::wget=1.20.3"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/kofamscan:1.3.0--hdfd78af_2':
    //    'quay.io/biocontainers/kofamscan:1.3.0--hdfd78af_2' }"

    input:
    tuple val(genus), path(input_aa_fasta)
    
    output:
    tuple val(genus), path("*_kofamscan.tsv"), emit: tsv

    script:
    def prefix = task.ext.prefix ?: "${genus}"
    """
    mkdir -p tmp
    # I decided to download the databases in place instead of supplying them as inputs.
    # ko_list.gz is 810Kb
    # profiles.tar.gz is 1.3Gb.
    wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz && tar xf profiles.tar.gz
    wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz && gunzip ko_list.gz
    exec_annotation --ko-list ko_list --profile profiles --cpu $task.cpus --format mapper -o ${prefix}_kofamscan.tsv ${input_aa_fasta}
    """
}
