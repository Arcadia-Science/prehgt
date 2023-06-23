process download_reference_genomes {
    tag "$genus"
    label 'process_single'
    // setting errorStrategy to ignore when exit 1 is seen makes it so that 
    // genera that don't have any matches on GenBank or RefSeq don't progress 
    // through the pipeline, but don't stop everything else from running
    errorStrategy { task.exitStatus == 1 ? 'ignore' : 'terminate' } 

    conda "$baseDir/envs/ncbi-genome-download.yml"
    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/ncbi-genome-download-patch:4c5c24e':
        '' }"
    //keep this here bc these will be the download instructions once the software is patched
    //conda "bioconda::ncbi-genome-download=0.3.1"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/ncbi-genome-download:0.3.1--pyh5e36f6f_0' :
    //    'quay.io/biocontainers/ncbi-genome-download:0.3.1--pyh5e36f6f_0' }" 

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
    # don't exit immediately if we get an exit 1, keep running the script.
    set +e
    # NOTE this needs to change to accomodate bacteria or archaea inputs.
    # I'm nervous about genus collisions in the NCBI taxonomy (esp happens with viruses, i don't know if it happens between bacteria <-> euks). 
    # I could change this to go from taxid, which will be unique, harder for a user to interpret at a glance.
    # I'll fix this later, for now it's fine.

    # genbank is separate from refseq, so run twice
    for section in genbank refseq
    do
        ncbi-genome-download vertebrate_mammalian,vertebrate_other,invertebrate,plant,fungi,protozoa --output-folder ./ --flat-output --genera ${prefix} -F gff,cds-fasta -s \$section --dry-run
        exit_code=\$?
        # if this gives an exit code 1 (ERROR: No downloads matched your filter. Please check your options.), skip that section as there are no outputs
        if [ \$exit_code -eq 1 ]; then
            echo "Skipping the next command..."
            (exit 0) # reset exit code. failing the dry run is not problematic
        # otherwise, do the download and get the genomes!
        else
            ncbi-genome-download vertebrate_mammalian,vertebrate_other,invertebrate,plant,fungi,protozoa --output-folder ./ --flat-output --genera ${prefix} -F gff,cds-fasta -s \$section --retries 3
        fi
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
