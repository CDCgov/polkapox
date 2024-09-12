process VCFTOOLS_DIFF {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::vcftools=0.1.16" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcftools:0.1.16--he513fc3_4' :
        'quay.io/biocontainers/vcftools:0.1.16--he513fc3_4' }"

    input:
    tuple val(meta), path(variant1)
    tuple val(meta), path(variant2)

    output:
    tuple val(meta), path("*.sites_in_files"), emit: diff
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}." + task.ext.caller
    def args_list = args.tokenize()

    """
    vcftools \\
        --vcf $variant1 \\
        --diff $variant2 \\
        --diff-site \\
        --out $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcftools: \$(echo \$(vcftools --version 2>&1) | sed 's/^.*VCFtools (//;s/).*//')
    END_VERSIONS
    """
}


