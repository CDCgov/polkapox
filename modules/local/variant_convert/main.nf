process VARIANT_CONVERT {
    tag "$meta.id"
    label 'process_medium'

        conda "${moduleDir}/environment.yml"
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bionumpy:0.2.17--pyha8f3691_0  ' :
        'quay.io/biocontainers/bionumpy:0.2.17--pyha8f3691_0' }"


    input:
    tuple val(meta), path(tsv)
    val(af_cutoff)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //g"'), topic: versions, emit: versions_python

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
    """
}
