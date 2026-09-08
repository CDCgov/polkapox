/*
 * Thin compatibility layer for an originally modified nf-core VCFTOOLS module.
 *
 * The old custom process was named VCFTOOLS_DIFF, accepted two VCF channels
 * tuple(meta, variant1) and tuple(meta, variant2), ran `vcftools --diff` with
 * `--diff-site`, and emitted `*.sites_in_files` as `diff`. This wrapper keeps
 * that workflow-facing interface while delegating the vcftools execution to the
 * upstream nf-core VCFTOOLS process. The diff arguments are configured for the
 * wrapped process in conf/modules.config.
 */

include { VCFTOOLS as VCFTOOLS_DIFF_RUN } from '../../modules/nf-core/vcftools/main'

workflow VCFTOOLS_DIFF {
    take:
    variant1
    variant2

    main:
    ch_empty_bed = Channel.value([])
    ch_diff_variant = variant2.map { meta, diff_variant -> diff_variant }

    VCFTOOLS_DIFF_RUN (
        variant1,
        ch_empty_bed,
        ch_diff_variant
    )

    emit:
    diff     = VCFTOOLS_DIFF_RUN.out.diff_sites_in_files
    versions = VCFTOOLS_DIFF_RUN.out.versions
}
