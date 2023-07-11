import pandas as pd
import csv
import os

#metadata = pd.read_csv("inputs/venoms.tsv", header = 0, sep = "\t")
metadata = pd.read_csv("inputs/candidate_fungi_for_bio_test_data_set.tsv", header = 0, sep = "\t")
source = ["genome"]
metadata = metadata.loc[metadata['source'].isin(source)] 
GENUS = metadata['genus'].unique().tolist()

# remove Trametes String from genus list
while("Trametes" in GENUS):
    GENUS.remove("Trametes")

# explanation of wildcards:
# genus (defined by GENUS): All of the genera that the pipeline will be executed on. This is defined from an input metadata file. 
# accession (inferred from checkpoint_download_reference_genomes): Finds all accessions for genomes associated with a given genus.

rule all:
    input: "outputs/hgt_candidates_final/all_results.tsv"
        
###################################################
## download references & build pangenome
###################################################

checkpoint download_reference_genomes:
    '''
    Using user-provided genera names, this rule uses the ncbi-genome-download tool to download the GFF annotation file and coding domain sequence (CDS) file for all member of that genus in GenBank and RefSeq.
    It then deletes any genomes that are in RefSeq and GenBank, keeping only the RefSeq version.
    Next, it combines the raw genome files for that genus into a single FASTA file.
    Lastly, it reports the genomes that were downloaded into a CSV file.
    While this rule could be broken out into many smaller rules, keeping it in this format keeps it consistent with the nextflow implementation and simplifies wildcards throughout the workflow.
    '''
    output: 
        outdir=directory("inputs/genbank/{genus}"),
        genome_csv="inputs/genbank/{genus}_genomes.csv",
        cds_combined="inputs/genbank/{genus}_cds.fna"
    conda: "envs/ncbi-genome-download.yml"
    params: outdir= lambda wildcards: "inputs/genbank/" + wildcards.genus
    benchmark: "benchmarks/download_reference_genomes/{genus}.tsv"
    shell:'''
    # NOTE this needs to change to accomodate bacteria or archaea inputs.
    # I'm nervous about genus collisions in the NCBI taxonomy (esp happens with viruses, i don't know if it happens between bacteria <-> euks). 
    # I could change this to go from taxid, which will be unique, harder for a user to interpret at a glance.
    # I'll fix this later, for now it's fine.

    # genbank is separate from refseq, so run twice
    for section in genbank refseq
    do
        ncbi-genome-download vertebrate_mammalian,vertebrate_other,invertebrate,plant,fungi,protozoa --output-folder {params.outdir} --flat-output --genera {wildcards.genus} -F gff,cds-fasta -s \$section --retries 3
    done

    # when downloading from RefSeq and GenBank sections, there might be duplicate
    # accessions that vary only by GCA_* or GCF_*, where GCF_* files are from RefSeq.
    # delete_gca_files.sh deletes GCA_* files only when there is a corresponding GCF_ file.
    cd ./inputs/genbank/{wildcards.genus}
    ./bin/delete_gca_files.sh
    cd ../../../

    # combine & unzip all genus-level gene sequences
    cat {params.outdir}/*_cds_from_genomic.fna.gz | gunzip > inputs/genbank/{wildcards.genus}_cds.fna

    # record the genomes were downloaded in a csv file
    echo "genus,genome" > {params.outdir}/{wildcards.genus}_genomes.csv
    ls {params.outdir}/*_cds_from_genomic.fna.gz | while read -r line; do
       bn=$(basename ${{line}})
       echo "{wildcards.genus},\${{bn}}" >> inputs/genbank/{wildcards.genus}_genomes.csv
    done
    '''

def checkpoint_download_reference_genomes(wildcards):
    # solve for gff file paths
    checkpoint_output = checkpoints.download_reference_genomes.get(**wildcards).output[0] # grab output directory name 
    file_names = expand("inputs/genbank/{{genus}}/{accession}_genomic.gff.gz",
                        accession = glob_wildcards(os.path.join(checkpoint_output, "{accession}_genomic.gff.gz")).accession)
    return file_names


rule combine_and_parse_gff_per_genus:
    '''
    Using the class checkpoint_accessions_to_genus, this rule combines all GFF annotation files from all accessions that belong to a given genus into one file.
    After pangenome clustering, this file will be used to retrieve information (genomic coords, etc) for rep sequences.
    '''
    input: gff = checkpoint_download_reference_genomes,
    output: gff = "outputs/genus_pangenome_raw/{genus}_gff_info.tsv"
    benchmark: "benchmarks/combine_gff_per_genus/{genus}.tsv"
    conda: "envs/tidy-prehgt.yml"
    shell:'''
    bin/combine_and_parse_gff_per_genus.R {output} {input}
    '''

rule build_genus_pangenome:
    '''
    This rule clusters all CDS sequences from all accessions in a given genus into a pangenome.
    This reduces redundancy for subsequent searches, and this information can be leveraged to interpret the HGT candidate genes that will eventually be predicted by the pipeline.
    selecting clustering threshold:
    0.9 https://www.science.org/doi/full/10.1126/sciadv.aba0111
    0.98, 0.99 https://www.sciencedirect.com/science/article/pii/S0960982220314263
    0.9 https://www.pnas.org/doi/abs/10.1073/pnas.2009974118
    '''
    input: "inputs/genbank/{genus}_cds.fna"
    output: 
        "outputs/genus_pangenome_clustered/{genus}_cds_rep_seq.fasta",
        "outputs/genus_pangenome_clustered/{genus}_cds_cluster.tsv"
    conda: "envs/mmseqs2.yml"
    benchmark:"benchmarks/genus_pangenome/{genus}.tsv"
    params: outprefix = lambda wildcards: "outputs/genus_pangenome_clustered/" + wildcards.genus + "_cds"
    shell:'''
    mmseqs easy-cluster {input} {params.outprefix} tmp_mmseqs2 --min-seq-id 0.9
    '''

rule translate_pangenome:
    '''
    Translate the representative pangenome sequences from nucleotide to amino acid
    '''
    input: "outputs/genus_pangenome_clustered/{genus}_cds_rep_seq.fasta"
    output: "outputs/genus_pangenome_clustered/{genus}_aa_rep_seq.fasta"
    conda: "envs/emboss.yml"
    benchmark: "benchmarks/translate_pangenome/{genus}.tsv"
    shell:'''
    transeq -sequence {input} -outseq {output}
    '''

##################################################
## DNA compositional screen
###################################################

rule compositional_scans_pepstats:
    '''
    This rule uses emboss pepstats to calculate relative amino acid usage per pangenome seq
    '''
    input: "outputs/genus_pangenome_clustered/{genus}_aa_rep_seq.fasta"
    output: 'outputs/compositional_scans_pepstats/{genus}_pepstats.txt'
    conda: "envs/emboss.yml"
    benchmark: "benchmarks/compositional_scans_pepstats/{genus}.tsv"
    shell:'''
    pepstats -sequence {input} -outfile {output}
    '''

rule compositional_scans_to_hgt_candidates:
    '''
    This script performs hierarchical clustering on genes based on their relative amino acid usage and identifies clusters with fewer than 150 genes. 
    The output is a list of genes that belong to these small clusters, which are considered HGT candidates.
    '''
    input: raau='outputs/compositional_scans_pepstats/{genus}_pepstats.txt',
    output: 
        tsv="outputs/compositional_scans_hgt_candidates/{genus}_clusters.tsv",
        gene_lst="outputs/compositional_scans_hgt_candidates/{genus}_pepstats_gene_lst.txt"
    benchmark: "benchmarks/compositional_scans_to_hgt_candidates/{genus}.tsv"
    conda: "envs/tidy-prehgt.yml"
    shell:'''
    bin/compositional_scans_to_hgt_candidates.R {input.raau} {output.tsv} {output.gene_lst}
    '''

###################################################
## BLASTP (diamond) against clustered nr
###################################################

rule blast_against_clustered_nr:
    '''
    This rule uses the diamond implementation of BLASTP to compare each CDS in the genus-level pangenome to a clustered version of NR.
    For more information on the database, see this repository: https://github.com/Arcadia-Science/2023-nr-clustering
    Using the diamond implementation and the clustered database decreases the time it takes to run this step. 
    '''
    input:
        db="inputs/nr_rep_seq.fasta.gz", # downloaded from S3...TBD on how to make available, its 60GB
        query="outputs/genus_pangenome_clustered/{genus}_aa_rep_seq.fasta"
    output: "outputs/blast_diamond/{genus}_vs_clustered_nr.tsv"
    conda: "envs/diamond.yml"
    benchmark: "benchmarks/blast_against_clustered_nr/{genus}.tsv"
    threads: 16
    shell:'''
    diamond blastp --db {input.db} --query {input.query} --out {output} \
        --outfmt 6 qseqid qtitle sseqid stitle pident approx_pident length mismatch gapopen qstart qend qlen qcovhsp sstart send slen scovhsp evalue bitscore score corrected_bitscore \
        --max-target-seqs 100 --threads {threads} --faster
    '''

rule blast_add_taxonomy_info:
    '''
    This rule adds taxonomic lineages to each BLAST match.
    The taxonomy sheet records the lowest common ancestor for a BLAST match, given that all BLAST matches represent a cluster of proteins.
    Because the taxonomy sheet is so large (>60GB), the script uses an sql query executed via dplyr and dbplyr to decrease search times.
    '''
    input: 
        tsv="outputs/blast_diamond/{genus}_vs_clustered_nr.tsv",
        sqldb="inputs/nr_cluster_taxid_formatted_final.sqlite" # downloaded from S3...TBD on how to make available, it's 63 GB
    output: tsv="outputs/blast_diamond/{genus}_vs_clustered_nr_lineages.tsv"
    conda: "envs/tidy-prehgt.yml"
    benchmark: "benchmarks/blast_add_taxonomy_info/{genus}.tsv"
    shell:'''
    bin/blastp_add_taxonomy_info.R {input.sqldb} {input.tsv} {output.tsv}
    '''

rule blast_to_hgt_candidates_kingdom:
    '''
    This script processes BLAST matches and their taxonomic lineages to identify HGT candidates using alien index, horizontal gene transfer index, donor distribution index, and acceptor lowest common acnestor calculations.
    It scores all candidates and writes the scores and other relevant information to a TSV file and outputs a list of candidate gene IDs.
    '''
    input: tsv="outputs/blast_diamond/{genus}_vs_clustered_nr_lineages.tsv"
    output: 
        gene_lst="outputs/blast_hgt_candidates/{genus}_blastp_kingdom_gene_lst.txt",
        tsv="outputs/blast_hgt_candidates/{genus}_blastp_kingdom_scores.tsv"
    conda: "envs/tidy-prehgt.yml"
    benchmark: "benchmarks/blast_to_hgt_candidates_kingdom/{genus}.tsv"
    shell:'''
    bin/blastp_to_hgt_candidates_kingdom.R {input.tsv} {output.tsv} {output.gene_lst}
    '''

rule blast_to_hgt_candidates_subkingdom:
    '''
    This script processes BLAST matches and their taxonomic lineages to identify HGT candidates using transfer index, a sub-kingdom identification algorithm
    It scores all candidates and writes the scores and other relevant information to a TSV file and outputs a list of candidate gene IDs.
    '''
    input: tsv="outputs/blast_diamond/{genus}_vs_clustered_nr_lineages.tsv"
    output: 
        gene_lst="outputs/blast_hgt_candidates/{genus}_blastp_subkingdom_gene_lst.txt",
        tsv="outputs/blast_hgt_candidates/{genus}_blastp_subkingdom_scores.tsv"
    conda: "envs/tidy-prehgt.yml"
    benchmark: "benchmarks/blast_to_hgt_candidates_subkingdom/{genus}.tsv"
    shell:'''
    bin/blastp_to_hgt_candidates_subkingdom.R {input.tsv} 0.01 {output.tsv} {output.gene_lst}
    '''

###################################################
## candidate characterization
###################################################

rule combine_hgt_candidates:
    '''
    This rule combines the lists of HGT candidate genes identified through two different methods - BLAST and compositional scans. 
    The output is a single list containing the unique genes from both input lists.
    * cat {input}: stream the input file to stdin
    * csvtk freq -H -f 1: for the headerless CSV that is being streamed in (it's a one column txt file so delimiter isn't important), count the frequency of each term in the first field.
      The frequency is reported as a new column, and duplicated terms are collapsed in the first column.
    * csvtk cut -f 1 -o {output}: cut the first column (the no de-duplicated hgt candidate names) and output it to a file.
    '''
    input: 
       "outputs/blast_hgt_candidates/{genus}_blastp_kingdom_gene_lst.txt",
       "outputs/blast_hgt_candidates/{genus}_blastp_subkingdom_gene_lst.txt",
       "outputs/compositional_scans_hgt_candidates/{genus}_pepstats_gene_lst.txt"
    output: "outputs/hgt_candidates/{genus}_blastp_gene_lst.txt"
    conda: "envs/csvtk.yml"
    shell:'''
    cat {input} | csvtk freq -H -f 1 | csvtk cut -f 1 -o {output}
    '''

rule extract_hgt_candidates:
    '''
    This rule extracts the DNA sequences of the HGT candidate genes from the pangenome clusters. 
    The output is a FASTA file containing the sequences of the identified HGT candidate genes.
    '''
    input:
        fa = "outputs/genus_pangenome_clustered/{genus}_aa_rep_seq.fasta", 
        gene_lst = "outputs/hgt_candidates/{genus}_gene_lst.txt"
    output: "outputs/hgt_candidates/{genus}_aa.fasta"
    benchmark: "benchmarks/extract_hgt_candidates/{genus}.tsv"
    conda: "envs/seqtk.yml"
    shell:'''
    seqtk subseq {input.fa} {input.gene_lst} > {output}
    '''

rule download_kofamscan_ko_list:
    """
    This rule downloads the kofamscan KEGG list file.
    """
    output: "inputs/kofamscandb/ko_list"
    params: outdir = "inputs/kofamscandb/"
    shell:'''
    wget -O {output}.gz ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz && gunzip {output}.gz -C {params.outdir}
    '''

rule download_kofamscan_profiles:
    """
    This rule downloads the kofamscan KEGG hmm profiles.
    """
    output: "inputs/kofamscandb/profiles/prokaryote.hal"
    params: outdir = "inputs/kofamscandb/"
    shell:'''
    wget -O {params.outdir}/profiles.tar.gz ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz && tar xf {params.outdir}/profiles.tar.gz -C {params.outdir}
    '''

rule kofamscan_hgt_candidates:
    """
    This rule uses the kofamscan to perform KEGG ortholog annotation on the HGT candidate genes. 
    """
    input:
        fa="outputs/hgt_candidates/{genus}_aa.fasta",
        kolist="inputs/kofamscandb/ko_list",
        profiles="inputs/kofamscandb/profiles/prokaryote.hal"
    output: "outputs/hgt_candidates_annotation/kofamscan/{genus}_kofamscan.tsv"
    conda: "envs/kofamscan.yml"
    params: profilesdir = "inputs/kofamscandb/profiles"
    threads: 8
    benchmark: "benchmarks/kofamscan_hgt_candidates/{genus}.tsv"
    shell:'''
    mkdir -p tmp
    exec_annotation --format detail-tsv --ko-list {input.kolist} --profile {params.profilesdir} --cpu {threads} -o {output} {input.fa}
    '''

rule hmmscan_hgt_candidates:
    '''
    Using the hmm file built in make_hmm_db.snakefile, this rule uses hidden markov models to annotate specific protein classes of interest.
    At the moment it targets viruses and biosynthetic gene clusters, but the HMM file can be expanded in the aforementioned snakefile as desired. 
    '''
    input:
        hmmdb="inputs/hmms/all_hmms.hmm",
        fa="outputs/hgt_candidates/{genus}_aa.fasta"
    output:
        tblout="outputs/hgt_candidates_annotation/hmmscan/{genus}.tblout",
        out="outputs/hgt_candidates_annotation/hmmscan/{genus}.out"
    conda: "envs/hmmer.yml"
    threads: 8
    benchmark: "benchmarks/hmmscan_hgt_candidates/{genus}.tsv"
    shell:'''
    hmmpress {input.hmmdb}
    hmmscan --cpu {threads} --tblout {output.tblout} -o {output.out} {input.hmmdb} {input.fa}
    '''

###################################################
## Combine results
###################################################

rule combine_results_genus:
    '''
    Combine all of the results into a single mega TSV file.
    The results are joined either on the genus or on the HGT candidate gene name, derived from the pangenome FASTA file.
    '''
    input: 
        compositional = "outputs/compositional_scans_hgt_candidates/{genus}_clusters.tsv",
        blast_kingdom = "outputs/blast_hgt_candidates/{genus}_blastp_kingdom_scores.tsv",
        blast_subkingdom = "outputs/blast_hgt_candidates/{genus}_blastp_subkingdom_scores.tsv",
        genome_csv = "inputs/genbank/{genus}_genomes.csv",
        pangenome_cluster = "outputs/genus_pangenome_clustered/{genus}_cds_cluster.tsv",
        gff = "outputs/genus_pangenome_raw/{genus}_gff_info.tsv",
        kofamscan = "outputs/hgt_candidates_annotation/kofamscan/{genus}_kofamscan.tsv",
        hmmscan = "outputs/hgt_candidates_annotation/hmmscan/{genus}.tblout",
    output: 
        all_results = "outputs/hgt_candidates_final/{genus}_results.tsv",
        method_tally = "outputs/hgt_candidates_final/{genus}_method_tally.tsv"
    conda: "envs/tidy-prehgt.yml"
    shell:'''
    bin/combine_results_genus.R {input.compositional} {input.blast_kingdom} {input.blast_subkingdom} {input.genome_csv} {input.pangenome_cluster} {input.gff} {input.kofamscan} {input.hmmscan}
    '''

rule combine_results:
    input: expand("outputs/hgt_candidates_final/{genus}_results.tsv", genus = GENUS)
    output: "outputs/hgt_candidates_final/all_results.tsv"
    conda: "envs/tidy-prehgt.yml"
    shell:'''
    bin/combine_results.R {input}
    '''
