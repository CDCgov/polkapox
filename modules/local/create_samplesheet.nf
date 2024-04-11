process CREATE_SAMPLESHEET {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container 'quay.io/biocontainers/python:3.8.3'

input:
path indir

output:
path("*.csv"),       emit: samplesheet
path "versions.yml", emit: versions

script: // This script is bundled with the pipeline, in polkapox/bin
single = params.paired ? "" : "--single"  // $single is empty string if paired=true, '--single' if paired=false
outfile = "${params.project_name}_samplesheet.csv" // Samplesheet name = basename of the project dir containing samples

"""
    create_samplesheet.py \\
        --indir $indir \\
        --outdir . \\
        --file_levels ${params.file_levels} \\
        $single \\
        --project_name ${params.project_name} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
