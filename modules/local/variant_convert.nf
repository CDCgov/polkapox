process VARIANT_CONVERT {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext. singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(tsv)
    val(af_cutoff)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-po '
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    ivar2vcf.py \\
        $args \\
        -af $af_cutoff \\
        $tsv \\
        ${prefix}.ivar.vcf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
