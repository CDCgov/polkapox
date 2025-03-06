include { SRATOOLS_PREFETCH           } from '../../../modules/nf-core/sratools/prefetch/main'
include { SRATOOLS_FASTERQDUMP        } from '../../../modules/nf-core/sratools/fasterqdump/main'
include { SRA_TO_SAMPLESHEET          } from '../../../modules/local/sra_to_samplesheet/sra_to_samplesheet'
//
// Download FASTQ sequencing reads from the NCBI's Sequence Read Archive (SRA).
//
workflow SRA_TOOLS {
    take:
    ch_sra_ids   // channel: [ val(meta), val(id) ]

    main:

    ch_versions = Channel.empty()

    //
    // Prefetch sequencing reads in SRA format.
    //
    SRATOOLS_PREFETCH ( ch_sra_ids )
    //ch_versions = ch_versions.mix(SRATOOLS_PREFETCH.out.versions.first())

    //
    // Convert the SRA format into one or more compressed FASTQ files.
    //
    SRATOOLS_FASTERQDUMP ( SRATOOLS_PREFETCH.out.sra )
    ch_versions = ch_versions.mix(SRATOOLS_FASTERQDUMP.out.versions.first())

    //
    // Convert reads ch to samplesheet
    //
    // SRA_TO_SAMPLESHEET (
    //     SRATOOLS_FASTERQDUMP.out.forward, 
    //     SRATOOLS_FASTERQDUMP.out.reverse 
    //     )

    emit:
    //samplesheet = SRA_TO_SAMPLESHEET.out.samplesheet 
    forward = SRATOOLS_FASTERQDUMP.out.forward
    reverse = SRATOOLS_FASTERQDUMP.out.reverse
    versions = ch_versions                    // channel: [ versions.yml ]
}
