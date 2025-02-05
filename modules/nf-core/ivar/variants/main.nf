process IVAR_VARIANTS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ivar=1.3.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ivar:1.3.1--h089eab3_0' :
        'quay.io/biocontainers/ivar:1.3.1--h089eab3_0' }"

    input:
    tuple val(meta), path(bam)
    val save_tsv

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}." + task.ext.caller
    def tsv = save_tsv ? "| tee ${prefix}.tsv" : ""

    """
    samtools \\
        mpileup \\
        --reference ${params.fasta} \\
        $args \\
        $bam \\
        | ivar \\
        variants \\
        $args2 \\
        -m ${params.DP_cutoff} \\
        -t ${params.af_cutoff} \\
        -r ${params.fasta} \\
        -g ${params.gff} \\
        -p $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(echo \$(ivar version 2>&1) | sed 's/^.*iVar version //; s/ .*\$//')
    END_VERSIONS
    """

}
