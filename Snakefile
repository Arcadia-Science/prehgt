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


metadata = pd.read_csv("inputs/candidate_fungi_for_bio_test_data_set.tsv", header = 0, sep = "\t")
source = ["genome"]
metadata = metadata.loc[metadata['source'].isin(source)] 
GENUS = metadata['genus'].unique().tolist()
#ACCESSION = metadata['accession'].unique().tolist()

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
    * sed 's|\(.*\)_.*|\1|' removes the last underscore and everything after it
    * sed 's/^\([^_]*\(_[^_]*\)\{1\}\)_/>/' removes everything before and including the second underscore and replaces with ">"
    double curly braces escape snakemake and are parsed as single curly braces in the shell command
    '''
    input: "inputs/genbank/{accession}_cds_from_genomic.fna.gz"
    output: "outputs/genbank/{accession}_cds_from_genomic.fna"
    benchmark: "benchmarks/shorten_gene_names/{accession}.tsv"
    shell:'''
    gunzip {input} | awk '{{print $1;next}}1' | sed 's|\(.*\)_.*|\1|' | sed 's/^\([^_]*\(_[^_]*\)\{{1\}}\)_/>/' > {output}
    '''

###################################################
## generate pangenome
###################################################

checkpoint accessions_to_genus:
    input: metadata="inputs/candidate_fungi_for_bio_test_data_set.tsv"
    output: genus="outputs/accessions_to_genus/{genus}.csv",
    conda: "envs/tidyverse.yml"
    benchmark: "benchmarks/accessions_to_genus/{genus}.tsv"
    script: "scripts/accessions_to_genus.R"

rule combine_cds_per_genus:
    input: checkpoint_accessions_to_genus('outputs/genbank/{accession}_cds_from_genomic.fna')
    output: "outputs/genus_pangenome_raw/{genus}_cds.fna"
    benchmark: "benchmarks/combine_cds_per_genus/{genus}.tsv"
    shell:'''
    cat {input} | gunzip > {output}
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
        bulk='outputs/compositional_scans_codonw/{genus}_rauu.txt'
    conda: "envs/codonw.yml"
    benchmark: "benchmarks/compositional_scans_codonw/{genus}.tsv"
    shell:'''
    codonw {input} {output.indices} {output.bulk} -all_indices -nomenu -machine -silent -raau
    '''

rule compositional_scans_bbmap:
    input: "outputs/genus_pangenome_clustered/{genus}_cds_rep_seq.fasta"
    output: "outputs/compositional_scans_bbmap/{genus}_tetramerfreq.tsv"
    conda: "envs/bbmap.yml"
    benchmark: "benchmarks/compositional_scans_bbmap/{genus}.tsv"
    shell:'''
    tetramerfreq.sh in={input} out={output} w=0 short=T k=4
    '''


###################################################
## BLASTP (diamond) against clustered nr
###################################################

###################################################
## candidate characterization
###################################################
