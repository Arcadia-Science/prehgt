process download_reference_genomes {
    tag "$genus"
    label 'process_single'

    conda "bioconda::ncbi-genome-download=0.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-genome-download:0.3.1--pyh5e36f6f_0' :
        'quay.io/biocontainers/ncbi-genome-download:0.3.1--pyh5e36f6f_0' }" 

    input:
    val genus

    output:
    tuple val(genus), path("*_cds_from_genomic.fna.gz"), emit: cds
    tuple val(genus), path("*_genomic.gff.gz")         , emit: gff
    tuple val(genus), path("*_cds.fna")                , emit: cds_all
    tuple val(genus), path("*_genomes.csv")            , emit: csv

    script:
    def prefix = task.ext.prefix ?: "${genus}"
    """
    # NOTE this needs to change to accomodate bacteria or archaea inputs.
    # I'm nervous about genus collisions in the NCBI taxonomy (esp happens with viruses, i don't know if it happens between bacteria <-> euks). 
    # I could change this to go from taxid, which will be unique, harder for a user to interpret at a glance.
    # I'll fix this later, for now it's fine.

    # genbank is separate from refseq, so run twice
    for section in genbank refseq
    do
        ncbi-genome-download vertebrate_mammalian,vertebrate_other,invertebrate,plant,fungi,protozoa --output-folder ./ --flat-output --genera ${prefix} -F gff,cds-fasta -s \$section --retries 3
    done

    # when downloading from RefSeq and GenBank sections, there might be duplicate
    # accessions that vary only by GCA_* or GCF_*, where GCF_* files are from RefSeq.
    # delete_gca_files.sh deletes GCA_* files only when there is a corresponding GCF_ file.
    delete_gca_files.sh
   
    # combine & unzip all genus-level gene sequences
    cat *_cds_from_genomic.fna.gz | gunzip > ${prefix}_cds.fna

    # record the genomes were downloaded in a csv file
    echo "genus,genome" > ${prefix}_genomes.csv
    ls *_cds_from_genomic.fna.gz | while read -r line; do
       echo "${prefix},\${line}" >> ${prefix}_genomes.csv
    done
    """
}