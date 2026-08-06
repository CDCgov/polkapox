
process SRA_TO_SAMPLESHEET {
    executor 'local'
    memory 100.MB

    input:
    tuple val(meta), path(fastq1)
    tuple val(meta), path(fastq2)

    output:
    tuple val(meta), path("*samplesheet.csv"), emit: samplesheet

    script:
    """
    echo "sample_id,fastq_1,fastq_2" > samplesheet.csv
    echo "${meta},${fastq1},${fastq2}" >> samplesheet.csv
    """

}
