process blastp_against_clustered_nr {
    tag "$genus"
    label 'process_blast'

    conda "bioconda::diamond=2.1.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/diamond:2.1.6--h5b5514e_1':
        'quay.io/biocontainers/diamond:2.1.6--h5b5514e_1' }"

    input:
    path(input_db)
    tuple val(genus), path(input_aa_rep_seq)

    output:
    tuple val(genus), path("*_vs_clustered_nr.tsv"), emit: tsv

    script:
    def prefix = task.ext.prefix ?: "${genus}"
    """
    diamond blastp --db ${input_db} --query ${input_aa_rep_seq} --out ${prefix}_vs_clustered_nr.tsv \
        --outfmt 6 qseqid qtitle sseqid stitle pident approx_pident length mismatch gapopen qstart qend qlen qcovhsp sstart send slen scovhsp evalue bitscore score corrected_bitscore \
        --max-target-seqs 100 --threads $task.cpus --faster
    """
}
