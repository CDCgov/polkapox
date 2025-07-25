process UNICYCLER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/unicycler:0.4.8--py38h8162308_3' :
        'quay.io/biocontainers/unicycler:0.4.8--py38h8162308_3' }"

    input:
    tuple val(meta), path(shortreads)

    output:
    tuple val(meta), path('*.scaffolds.fa.gz'), emit: scaffolds optional true
    tuple val(meta), path('*bridges_applied.gfa'), emit: gfa optional true
    tuple val(meta), path('*.log')            , emit: log
    path  "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def short_reads = shortreads ? ( meta.single_end ? "-s $shortreads" : "-1 ${shortreads[0]} -2 ${shortreads[1]}" ) : ""
    """
    unicycler \\
        --threads $task.cpus \\
        $args \\
        $short_reads \\
        --out ./ \\
        --keep 2

    mv unicycler.log ${prefix}.unicycler.log
    mv assembly.fasta ${prefix}.scaffolds.fa
    gzip -n ${prefix}.scaffolds.fa
    mv *_bridges_applied.gfa ${prefix}.bridges_applied.gfa
    mv assembly.gfa ${prefix}.assembly.gfa
    gzip -n ${prefix}.assembly.gfa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        unicycler: \$(echo \$(unicycler --version 2>&1) | sed 's/^.*Unicycler v//; s/ .*\$//')
    END_VERSIONS
    """
}

