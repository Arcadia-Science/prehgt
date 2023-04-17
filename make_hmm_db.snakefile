import pandas as pd
import shutil

hmm_urls = pd.read_csv("inputs/hmms/hmm_urls.csv", header = 0)
HMMS = hmm_urls['hmm'].tolist()

rule all:
    input: "inputs/hmms/all_hmms.hmm.h3i"

rule download_vog_hmms:
    output: "inputs/hmms/vog.hmm.tar.gz"
    shell:'''
    curl -JLo {output} http://fileshare.csb.univie.ac.at/vog/latest/vog.hmm.tar.gz
    '''

checkpoint decompress_vog_hmms:
    input: "inputs/hmms/vog.hmm.tar.gz"
    output: directory("inputs/hmms/vogs/")
    params: outdir = "inputs/hmms/vogs/"
    shell:'''
    mkdir -p {params.outdir}
    tar xf {input} -C {params.outdir}
    '''

rule download_single_hmms:
    output: "inputs/hmms/single_hmms/raw/{hmm}.hmm"
    run:
        hmm = wildcards.hmm
        hmm_df = hmm_urls.loc[(hmm_urls['hmm'] == wildcards.hmm)]
        if hmm_df is None:
            raise TypeError("'None' value provided for hmm_df. Are you sure the hmm df was not empty?")

        prefix_url = hmm_df['prefix_url'].values[0]
        shell("curl -JLo {output} {prefix_url}/{hmm}.hmm")


rule convert_single_hmms:
   input: "inputs/hmms/single_hmms/raw/{hmm}.hmm"
   output: "inputs/hmms/single_hmms/hmmconvert/{hmm}.hmm"
   conda: "envs/hmmer.yml"
   shell:'''
   hmmconvert {input} > {output}
   '''

def checkpoint_decompress_vog_hmms(wildcards):
    # expand checkpoint to get VOG hmm prefixes, and place them in the final file name that uses that wildcard
    # checkpoint_output encodes the output dir from the checkpoint rule. 
    checkpoint_output = checkpoints.decompress_vog_hmms.get(**wildcards).output[0]    
    file_names = expand("inputs/hmms/vogs/{vog}.hmm",
                        vog = glob_wildcards(os.path.join(checkpoint_output, "{vog}.hmm")).vog)
    return file_names    


rule combine_hmms:
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
    input: "inputs/hmms/all_hmms.hmm"
    output: "inputs/hmms/all_hmms.hmm.h3i"
    conda: "envs/hmmer.yml"
    shell:'''
    hmmpress {input}
    '''
