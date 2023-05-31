process eggnog_hgt_candidates {
    tag "$genus"
    label 'process_blast'

    conda "bioconda::eggnog-mapper=2.1.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.10--pyhdfd78af_0':
        'quay.io/biocontainers/eggnog-mapper:2.1.10--pyhdfd78af_0' }"

    input:
    path(eggnog_db)
    path(eggnog_dmnd)
    tuple val(genus), path(input_aa_fasta)
    
    output:
    tuple val(genus), path("*.emapper.annotations"), emit: annotations

    script:
    def prefix = task.ext.prefix ?: "${genus}"
    """
    mkdir -p tmpdir
    # the data_dir must contain the EGGNOG db files. Since nextflow stages the input database files in the run dir, the data_dir becomes the working dir. 
    emapper.py --cpu $task.cpus -i ${input_aa_fasta} --output ${prefix} -m diamond --tax_scope none --seed_ortholog_score 60 --override --temp_dir tmpdir --data_dir .
    """
}
