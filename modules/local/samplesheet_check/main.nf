process SAMPLESHEET_CHECK {
    tag "$samplesheet"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path samplesheet

    output:
    path '*.csv'       , emit: csv
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //g"'), topic: versions, emit: versions_python

    script: // This script is bundled with the pipeline, in polkapox/bin/
    """
    check_samplesheet.py \\
        $samplesheet \\
        samplesheet.valid.csv
    """
}
