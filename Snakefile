import pandas as pd
import csv
import os

class checkpoint_accessions_to_genus:
    """
    Define a class to bind accession to genus for pangenome calculation. 
    This approach is documented at this url:
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern):
        self.pattern = pattern

    def get_genome_accessions(self, genus):
        genus_csv = f'outputs/accessions_to_genus/{genus}.csv'
        assert os.path.exists(genus_csv)

        genome_accessions = []
        with open(genus_csv, 'rt') as fp:
           r = csv.DictReader(fp)
           for row in r:
               accession = row['accession'].split(' ')[0]
               genome_accessions.append(accession)

        return genome_accessions

    def __call__(self, w):
        global checkpoints
        checkpoints.accessions_to_genus.get(**w) # run when rule accessions_to_genus is finished
        genome_accessions = self.get_genome_accessions(w.genus) # parse accessions in gather output file
        p = expand(self.pattern, accession=genome_accessions, **w)
        return p


metadata = pd.read_csv("inputs/spomb_pos_control.tsv", header = 0, sep = "\t")
source = ["genome"]
metadata = metadata.loc[metadata['source'].isin(source)] 
GENUS = metadata['genus'].unique().tolist()

rule all:
    input:
        expand("outputs/genus_pangenome_clustered/{genus}_cds_rep_seq.fasta", genus = GENUS),
        expand("outputs/compositional_scans_codonw/{genus}_indices.txt", genus = GENUS),
        expand("outputs/compositional_scans_bbmap/{genus}_tetramerfreq.tsv", genus = GENUS)

###################################################
## download references
###################################################

rule download_reference_genomes:
    output: 
        cds="inputs/genbank/{accession}_cds_from_genomic.fna.gz",
        gff="inputs/genbank/{accession}_genomic.gff.gz",
        metadata="inputs/genbank/{accession}.tsv"
    conda: "envs/ncbi-genome-download.yml"
    params: outdir="inputs/genbank/"
    benchmark: "benchmarks/download_reference_genomes/{accession}.tsv"
    shell:'''
    ncbi-genome-download all -s genbank -o {params.outdir} --flat-output -m {output.metadata} --assembly-accessions {wildcards.accession} -F gff,cds-fasta
    mv {params.outdir}{wildcards.accession}*_cds_from_genomic.fna.gz {output.cds}
    mv {params.outdir}{wildcards.accession}*_genomic.gff.gz {output.gff}
    '''

rule shorten_gene_names_for_codonw:
    '''
    codonw truncates multifasta names (output at 25 characters, bulk output at 20 characters).
    The protein names should be unique identifiers.
    The following code parses the FASTA headers so they end up as the protein ids.
    This also matches annotations in the GFF file.
    * awk '{print $1;next}1' removes everything after the first space
    * sed 's/lcl|//g' removes the prefix lcl|
    * sed 's/_cds//g' removes the string _cds if it exists in the header
    * cut -d_ -f1,2 removes the second underscore and everything after it
    double curly braces escape snakemake and are parsed as single curly braces in the shell command
    '''
    input: "inputs/genbank/{accession}_cds_from_genomic.fna.gz"
    output: "outputs/genbank/{accession}_cds_from_genomic.fna"
    benchmark: "benchmarks/shorten_gene_names/{accession}.tsv"
    shell:'''
    gunzip -c {input} | awk '{{print $1;next}}1' | sed 's/lcl|//g' | sed 's/_cds//g' | cut -d_ -f1,2 > {output}
    '''

###################################################
## generate pangenome
###################################################

checkpoint accessions_to_genus:
    input: metadata="inputs/spomb_pos_control.tsv"
    output: genus="outputs/accessions_to_genus/{genus}.csv",
    conda: "envs/tidyverse.yml"
    benchmark: "benchmarks/accessions_to_genus/{genus}.tsv"
    script: "scripts/accessions_to_genus.R"

rule combine_cds_per_genus:
    input: checkpoint_accessions_to_genus('outputs/genbank/{accession}_cds_from_genomic.fna')
    output: "outputs/genus_pangenome_raw/{genus}_cds.fna"
    benchmark: "benchmarks/combine_cds_per_genus/{genus}.tsv"
    shell:'''
    cat {input} > {output}
    '''

rule build_genus_pangenome:
    '''
    selecting clustering threshold:
    0.9 https://www.science.org/doi/full/10.1126/sciadv.aba0111
    0.98, 0.99 https://www.sciencedirect.com/science/article/pii/S0960982220314263
    0.9 https://www.pnas.org/doi/abs/10.1073/pnas.2009974118
    '''
    input: "outputs/genus_pangenome_raw/{genus}_cds.fna"
    output: "outputs/genus_pangenome_clustered/{genus}_cds_rep_seq.fasta"
    conda: "envs/mmseqs2.yml"
    benchmark:"benchmarks/genus_pangenome/{genus}.tsv"
    params: outprefix = lambda wildcards: "outputs/genus_pangenome_clustered/" + wildcards.genus + "_cds"
    shell:'''
    mmseqs easy-cluster {input} {params.outprefix} tmp_mmseqs2 --min-seq-id 0.9
    '''

###################################################
## DNA compositional screen
###################################################

rule compositional_scans_codonw:
    input: "outputs/genus_pangenome_clustered/{genus}_cds_rep_seq.fasta"
    output: 
        indices="outputs/compositional_scans_codonw/{genus}_indices.txt",
        raau='outputs/compositional_scans_codonw/{genus}_raau.txt'
    conda: "envs/codonw.yml"
    benchmark: "benchmarks/compositional_scans_codonw/{genus}.tsv"
    shell:'''
    codonw {input} {output.indices} {output.raau} -all_indices -nomenu -machine -silent -raau
    '''

rule compositional_scans_to_hgt_candidates:
    input:
        raau='outputs/compositional_scans_codonw/{genus}_raau.txt',
    output: 
        tsv="outputs/compositional_scans_hgt_candidates/{genus}_clusters.tsv"
        gene_lst="outputs/compositional_scans_hgt_candidates/{genus}_gene_lst.txt",
    benchmark: "benchmarks/compositional_scans_to_hgt_candidates/{genus}.tsv"
    conda: "envs/tidyverse.yml"
    script: "scripts/compositional_scans_to_hgt_candidates.R"

###################################################
## BLASTP (diamond) against clustered nr
###################################################

rule blast_against_clustered_nr:
    input:
        db="inputs/nr_rep_seq.fasta.gz", # downloaded from S3...TBD on how to make available, its 60GB
        query="outputs/genus_pangenome_clusters/{genus}_cds_rep_seq.fasta"
    output: "outputs/blast_diamond/{genus}_vs_clustered_nr.tsv"
    conda: "envs/diamond.yml"
    benchmark: "benchmarks/blast_against_clustered_nr/{genus}.tsv"
    threads: 16
    shell:'''
    diamond blastp --db {input.db} --query {input.query} --out {output} --outfmt 6 --max-target-seqs 100 --threads {threads} --faster
    '''

rule blast_add_taxonomy_info:
    input: 
        tsv="outputs/blast_diamond/{genus}_vs_clustered_nr.tsv",
        sqldb="inputs/nr_cluster_taxid_formatted_final.sqlite" # downloaded from S3...TBD on how to make available, it's 63 GB
    output: tsv="outputs/blast_diamond/{genus}_vs_clustered_nr_lineages.tsv"
    params: sqldb_tbl="nr_cluster_taxid_table"
    conda: "envs/r-sql.yml"
    benchmark: "benchmarks/blast_add_taxonomy_info/{genus}.tsv"
    script: "scripts/blast_add_taxonomy_info.R"

rule blast_to_hgt_candidates:
    input: tsv="outputs/blast_diamond/{genus}_vs_clustered_nr_lineages.tsv"
    output: 
        gene_lst="outputs/blast_hgt_candidates/{genus}_gene_lst.txt",
        tsv="outputs/blast_hgt_candidates/{genus}_blast_scores.tsv"
    conda: "envs/tidyverse.yml"
    benchmark: "outputs/blast_to_hgt_candidates/{genus}.tsv"
    script: "scripts/blast_to_hgt_candidates.R"

###################################################
## candidate characterization
###################################################

rule combine_hgt_candidates:
    input: 
       "outputs/blast_hgt_candidates/{genus}_gene_lst.txt",
       "outputs/compositional_scans_hgt_candidates/{genus}_gene_lst.txt"
    output: "outputs/hgt_candidates/{genus}_gene_lst.txt"
    conda: "envs/csvtk.yml"
    shell:'''
    cat {input} | csvtk summary -H -f 1:uniq -o {output}
    '''

rule extract_hgt_candidates:
    input:
        fa = "outputs/genus_pangenome_clusters/{genus}_cds_rep_seq.fasta", 
        gene_lst = "outputs/hgt_candidates/{genus}_gene_lst.txt"
    output: "outputs/hgt_candidates/{genus}_cds.fasta"
    benchmark: "benchmarks/extract_hgt_candidates/{genus}.tsv"
    conda: "envs/seqtk.yml"
    shell:'''
    seqtk subseq {input.fa} {input.gene_lst} > {output}
    '''

rule download_eggnog_db:
    output: "inputs/eggnog_db/eggnog.db"
    conda: "envs/eggnog.yml"
    shell:'''
    download_eggnog_data.py -H -d 2 -y --data_dir inputs/eggnog_db
    '''

rule eggnog_hgt_candidates:
    input:
        db="inputs/eggnog_db/eggnog.db",
        fa="outputs/hgt_candidates/{genus}_cds.fasta"
    output: "outputs/hgt_candidates_annotation/eggnog/{genus}.emapper.annotations",
    conda: "envs/eggnog.yml"
    params:
        outdir="outputs/hgt_candidates_annotation/eggnog/",
        dbdir = "inputs/eggnog_db" 
    threads: 4
    benchmark: "benchmarks/eggnog_hgt_candidates/{genus}.tsv"
    shell:'''
    emapper.py --cpu {threads} -i {input.fa} --output {genus} \
       --output_dir {params.outdir} -m hmmer -d none --tax_scope none \
       --seed_ortholog_score 60 --override --temp_dir tmp/ \
       --data_dir {params.dbdir} --itype CDS --translate True
    '''

rule download_antismash_hmms:
    input:
    output:
    shell:'''
    
    '''
