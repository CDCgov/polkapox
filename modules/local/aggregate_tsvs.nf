process AGGREGATE_TSVS {
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bionumpy:0.2.17--pyha8f3691_0  ' :
        'quay.io/biocontainers/bionumpy' }"
                
    input:
    path(individual_tsvs)
    val save_aggr

    output:
    path("*.txt"), emit: txt
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def save_aggr = save_aggr ? "| tee all_samples.vcf.summary.txt" : ""
    """
    aggregate_tsvs.py
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
