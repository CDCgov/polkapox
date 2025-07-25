process SUMMARIZE_QC {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    path summarize_qc_files
    path samplesheet

    output:
    path '*.tsv'       , emit: tsv
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in polkapox/bin/
    def args = task.ext.args   ?: ''
    // Convert relative path to absolute path
    def absolute_outdir = file(params.outdir)
    """
    summarize_qc.py \\
        --analysis_dir . \\
        --samplesheet $samplesheet \\
        --reference_genome ${params.fasta} \\
        --project_outdir ${absolute_outdir} \\
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
