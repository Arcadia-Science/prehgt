process build_genus_pangenome {
    tag "$genus"
    label 'process_medium'

    conda "bioconda::mmseqs2=14.7e284"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321hf1761c0_0':
        'quay.io/biocontainers/mmseqs2:14.7e284--pl5321hf1761c0_0' }"

    input:
    tuple val(genus), path(input_cds)

    output:
    tuple val(genus), path("*_rep_seq.fasta"), emit: rep_seq
    tuple val(genus), path("*_cluster.tsv")  , emit: cluster

    script:
    def prefix = task.ext.prefix ?: "${genus}"
    """
    mkdir -p tmp_mmseqs2 # create tmp dir for mmseqs if it doesn't exist already
    mmseqs easy-cluster ${input_cds} ${prefix} tmp_mmseqs2 --min-seq-id 0.9 --cov-mode 1 --cluster-mode 2
    """
}
