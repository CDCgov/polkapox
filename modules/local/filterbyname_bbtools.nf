process BBMAP_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bbmap=38.90" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:38.90--he522d1c_1' :
        'quay.io/biocontainers/bbmap:38.90--he522d1c_1' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(mpxreads)

    output:
    tuple val(meta), path('*.fq.gz')         , emit: reads
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = meta.single_end ? "in=${reads[0]}" : "in=${reads[0]} in2=${reads[1]}"
    def filtered  = meta.single_end ? "out=${prefix}.fastq.gz" : "out=${prefix}_1.fq.gz out2=${prefix}_2.fq.gz"
    """
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    filterbyname.sh \
        $raw \
        $filtered \
        names=${mpxreads} \
        include=t

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh)
    END_VERSIONS
    """
}
