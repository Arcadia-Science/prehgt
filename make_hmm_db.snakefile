"""
This snakefile downloads HMMs from various sources and creates a single file that can be searched with hmmscan (hmmer3).
The goal of this file is to create an extensible framework for compiling HMMs from a variety of sources into one file that can be used for one annotation run.
This snakefile currently targets viruses and biosynthetic gene clusters, but can be extended to include new classes of genes as those use cases crop up.
This annotation approach is a complement to the eggnog annotation that also occurs in the main Snakefile.
"""

import pandas as pd
import shutil

hmm_urls = pd.read_csv("inputs/hmms/hmm_urls.csv", header = 0)
HMMS = hmm_urls['hmm'].tolist()

rule all:
    input: "inputs/hmms/all_hmms.hmm.h3i"

rule download_vog_hmms:
    '''
    Download viral HMMs from VOGDB (Virus Orthologous Groups)
    '''
    output: "inputs/hmms/vog.hmm.tar.gz"
    shell:'''
    curl -JLo {output} http://fileshare.csb.univie.ac.at/vog/latest/vog.hmm.tar.gz
    '''

checkpoint decompress_vog_hmms:
    '''
    Decompress the viral HMMs. 
    Use a checkpoint to re-evaluate the DAG after decompression, automatically inferring the VOG wildcard values.
    '''
    input: "inputs/hmms/vog.hmm.tar.gz"
    output: directory("inputs/hmms/vogs/")
    params: outdir = "inputs/hmms/vogs/"
    shell:'''
    mkdir -p {params.outdir}
    tar xf {input} -C {params.outdir}
    '''

def checkpoint_decompress_vog_hmms(wildcards):
    # expand checkpoint to get VOG hmm prefixes, and place them in the final file name that uses that wildcard
    # checkpoint_output encodes the output dir from the checkpoint rule. 
    checkpoint_output = checkpoints.decompress_vog_hmms.get(**wildcards).output[0]    
    file_names = expand("inputs/hmms/vogs/{vog}.hmm",
                        vog = glob_wildcards(os.path.join(checkpoint_output, "{vog}.hmm")).vog)
    return file_names    

rule download_single_hmms:
    '''
    From the CSV metadata file, download single HMMs.
    Sources include antismash and arcadia-created HMMs.
    '''
    output: "inputs/hmms/single_hmms/raw/{hmm}.hmm"
    run:
        hmm = wildcards.hmm
        hmm_df = hmm_urls.loc[(hmm_urls['hmm'] == wildcards.hmm)]
        if hmm_df is None:
            raise TypeError("'None' value provided for hmm_df. Are you sure the hmm df was not empty?")

        prefix_url = hmm_df['prefix_url'].values[0]
        shell("curl -JLo {output} {prefix_url}/{hmm}.hmm")


rule convert_single_hmms:
   '''
   Convert HMMs to the format for the most recent version of hmmer.
   This is required for all hmms to be combined in a single file, which enables all vs all annotation via hmmscan.
   '''
   input: "inputs/hmms/single_hmms/raw/{hmm}.hmm"
   output: "inputs/hmms/single_hmms/hmmconvert/{hmm}.hmm"
   conda: "envs/hmmer.yml"
   shell:'''
   hmmconvert {input} > {output}
   '''


rule combine_hmms:
    '''
    Combine all hmms into a single file.
    '''
    input:
        vogs = checkpoint_decompress_vog_hmms,
        single_hmms = expand("inputs/hmms/single_hmms/hmmconvert/{hmm}.hmm", hmm = HMMS),
    output: "inputs/hmms/all_hmms.hmm"
    run:
        with open(str(output[0]), 'wb') as wfd:
            for f in input:
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
 
rule hmmpress:
    '''
    Press the hmm file to prepare it for search via hmmscan.
    '''
    input: "inputs/hmms/all_hmms.hmm"
    output: "inputs/hmms/all_hmms.hmm.h3i"
    conda: "envs/hmmer.yml"
    shell:'''
    hmmpress {input}
    '''
