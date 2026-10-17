process CREATE_SAMPLESHEET {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/python:3.8.3'

input:
path indir

output:
path("*.csv"),       emit: samplesheet
tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //g"'), topic: versions, emit: versions_python

script: // This script is bundled with the pipeline, in polkapox/bin
single = params.paired ? "" : "--single"  // $single is empty string if paired=true, '--single' if paired=false
// Check if params.project_name is not empty and not null, then add --project_name flag
project_name = params.project_name ? "--project_name ${params.project_name}" : ""


"""
    create_samplesheet.py \\
        --indir $indir \\
        --outdir . \\
        --file_levels ${params.file_levels} \\
        $single \\
        $project_name
    """
}
