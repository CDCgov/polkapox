process SUMMARIZE_QC {

    conda (params.enable_conda ? "conda-forge::python=3.8.3 conda-forge::pandas" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas' }"

    input:
    path summarize_qc_files
    path samplesheet

    output:
    path '*.tsv'       , emit: tsv
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in polkapox/bin/
    def args = task.ext.args   ?: ''
    """
    summarize_qc.py \\
        --analysis_dir . \\
        --samplesheet $samplesheet \\
        --reference_genome ${params.fasta} \\
        --project_outdir ${params.outdir} \\
        --kraken_db ${params.kraken_db} \\
        --kraken_tax_ids ${params.kraken2_tax_ids} \\
        --filter ${params.filter} \\
        --workflow ${params.workflow} \\
        --coords ${params.coords} \\
        --locus1 ${params.locus1_name} \\
        --locus2 ${params.locus2_name}
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
