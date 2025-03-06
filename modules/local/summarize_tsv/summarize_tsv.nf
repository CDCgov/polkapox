process SUMMARIZE_TSV {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bionumpy:0.2.17--pyha8f3691_0  ' :
        'quay.io/biocontainers/bionumpy:0.2.17--pyha8f3691_0' }"
        
    input:
    tuple val(meta), path(ivar_tsv)
    tuple val(coord1_start), val(coord1_end), val(coord2_start), val(coord2_end)
    val save_vars

    output:
    tuple val(meta), path("*_ivar_summary.txt"), emit: vars
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vars = save_vars ? "| tee ${prefix}_ivar_summary.txt" : ""
    
    """
    filter_tsv.py \\
        -t $ivar_tsv \\
        -s $meta.id \\
        -fs $coord1_start\\
        -fe $coord1_end\\
        -es $coord2_start\\
        -ee $coord2_end\\
        -af ${params.af_cutoff} \\
        -DP ${params.DP_cutoff} \\
        -AD ${params.altDP_cutoff} \\
        -aq ${params.altQ_cutoff} \\
        -r ${params.outdir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
