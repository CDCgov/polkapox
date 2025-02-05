process SRATOOLS_PREFETCH {
    tag "$id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5f89fe0cd045cb1d615630b9261a1d17943a9b6a:6a9ff0e76ec016c3d0d27e0c0d362339f2d787e6-0' :
        'ncbi/sra-tools:2.11.3' }"

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
