process combine_results {
    //tag "$genus"
    label 'process_low'

    conda "conda-forge::r-tidyverse=2.0.0 conda-forge::r-janitor=2.2.0"
    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/tidy-prehgt:2.0.0':
        '' }"

    input:
    path(tsv_files)

    output:
    path("all_results.tsv"), emit: all_results

    script:
    """
    combine_results.R all_results.tsv *_results.tsv
    """
}
