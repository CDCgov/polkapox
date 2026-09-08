/*
 * Thin compatibility layer for an originally modified nf-core SAMTOOLS_FLAGSTAT
 * module.
 *
 * The old de novo-specific process wrote `${meta.id}.denovo.flagstat` and
 * emitted legacy `versions.yml`. This wrapper keeps the de novo workflow-facing
 * name and outputs while delegating flagstat execution to the upstream nf-core
 * SAMTOOLS_FLAGSTAT process. The `.denovo` output prefix is configured in
 * conf/modules.config.
 */

include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_DENOVO_RUN } from '../../modules/nf-core/samtools/flagstat/main'

workflow SAMTOOLS_FLAGSTAT_DENOVO {
    take:
    bambai

    main:
    SAMTOOLS_FLAGSTAT_DENOVO_RUN (
        bambai
    )

    emit:
    flagstat = SAMTOOLS_FLAGSTAT_DENOVO_RUN.out.flagstat
    versions = SAMTOOLS_FLAGSTAT_DENOVO_RUN.out.versions
}
