process SRATOOLS_PREFETCH {
    tag "$id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-tools:3.1.0--h9f5acd7_0' :
        'quay.io/biocontainers/sra-tools:3.1.0--h9f5acd7_0' }"

    input:
    tuple val(meta), val(id)

    output:
    tuple val(meta), path(id), emit: sra
    //path 'versions.yml'      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    """
    prefetch $id -O .

    """
}