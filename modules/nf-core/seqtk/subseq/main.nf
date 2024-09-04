process SEQTK_SUBSEQ {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::seqtk=1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3' :
        'quay.io/biocontainers/seqtk:1.3--h5bf99c6_3' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(classifiedreads)

    output:
    tuple val(meta), path('*.fq.gz')         , emit: reads
    path "versions.yml"                      , emit: versions
    tuple val(meta), path('*.opxreads.txt')  , emit: opxv_reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ext = "fa"
    if ("$reads" ==~ /.+\.fq|.+\.fq.gz|.+\.fastq|.+\.fastq.gz/) {
        ext = "fq"
    }
    """
    awk 'NR==FNR { tax_ids[\$1]++; next} \$3 in tax_ids{print \$2}' ${params.kraken2_tax_ids} ${classifiedreads} > ${prefix}.opxreads.txt
    
    seqtk \
        subseq \
        $args \
        ${reads[0]} \
        ${prefix}.opxreads.txt | \
        gzip --no-name > ${prefix}_1.${ext}.gz
        
    if [ "${meta.single_end}" == "false" ] 
    then
        seqtk \
        subseq \
        $args \
        ${reads[1]} \
        ${prefix}.opxreads.txt | \
        gzip --no-name > ${prefix}_2.${ext}.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
