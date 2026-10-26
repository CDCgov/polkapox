process AGGREGATE_TSVS {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bionumpy:0.2.17--pyha8f3691_0  ' :
        'quay.io/biocontainers/bionumpy:0.2.17--pyha8f3691_0' }"
                
    input:
    path(individual_tsvs)
    val save_aggr

    output:
    path("*.txt"), emit: txt
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //g"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def tee_output = save_aggr ? "| tee all_samples.vcf.summary.txt" : ""
    """
    aggregate_tsvs.py
    """
}
