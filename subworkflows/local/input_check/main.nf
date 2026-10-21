//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../../modules/local/samplesheet_check/main'
include { FQ_LINT           } from '../../../modules/nf-core/fq/lint/main'

workflow INPUT_CHECK {
    take:
    //file samplesheet from samplesheet_to_check
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    FQ_LINT(reads)

    // check that the fastq files are not empty
    FQ_LINT.out.lint.map { meta, lintReport ->
        checkFqLintReport(meta, lintReport)
    }

    validated_reads = reads
        .join(FQ_LINT.out.lint)
        .map { meta, fastq, lint -> [meta, fastq] }

    emit:
    reads = validated_reads                   // channel: [ val(meta), [ reads ] ]
    lint = FQ_LINT.out.lint                    // channel: [ val(meta), lint_report ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}

/**
 * Checks an fq lint report for the number of records
 *
 * Errors if the report does not contain a recognizable count line,
 * or if the reported record count is zero.
 *
 * @param meta       Sample metadata map, expected to include an `id` field
 * @param lintReport Path to the fq lint report file
 * @return           No value; this function only validates and logs
 */
def checkFqLintReport(meta, lintReport) {
    // fq writes colored log lines; strip ANSI codes so the regex can match
    def logText = lintReport.text.replaceAll(/\u001B\[[0-9;]*m/, '')

    // fq logs "read N * 2 records" (pass 1, paired only) and "read N records" (pass 2); keep the latter
    def counts = (logText =~ /(?m)read\s+(\d+)\s+records\s*$/).collect { it[1] as long }

    if (counts.isEmpty()) {
        error("Could not find a record count in the fq lint report for sample '${meta.id}': ${lintReport}")
    }

    if (counts.last() == 0) {
        error("fq lint found zero records for sample '${meta.id}'. The FastQ file(s) appear to be empty: ${lintReport}")
    }
}
