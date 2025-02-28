
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { KRAKEN2                                       } from '../../modules/nf-core/kraken2/main'
include { FASTP                                         } from '../../modules/nf-core/fastp/main'
include { SEQTK_SUBSEQ                                  } from '../../modules/nf-core/seqtk/subseq/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow READ_FILTER {
    take:
    input_reads

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Run Kraken to keep only orthopox reads 
    //
    ch_kraken2_db = file(params.kraken_db, checkIfExists: true)

    KRAKEN2 (
        input_reads,
        ch_kraken2_db,
        true,
        false,
        true
    )
    ch_kreads = KRAKEN2.out.classified_reads_fastq
    ch_orthoreads = KRAKEN2.out.classified_reads_assignment
    ch_versions = ch_versions.mix(KRAKEN2.out.versions.first().ifEmpty(null))

    SEQTK_SUBSEQ (
        ch_kreads,
        ch_orthoreads
    )
    ch_filt_fastq = SEQTK_SUBSEQ.out.reads
    ch_versions = ch_versions.mix(SEQTK_SUBSEQ.out.versions)
    
    //
    // MODULE: Run Fastp
    //

    FASTP (
        ch_filt_fastq,
        false,
        false
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    emit:
    trimmed_fastq = FASTP.out.reads 
    json = FASTP.out.json
    kraken2_report = KRAKEN2.out.report
    classified_reads_assignment = KRAKEN2.out.classified_reads_assignment
    seqtk_reads = SEQTK_SUBSEQ.out.opxv_reads
    versions      = ch_versions // channel: [ versions.yml ]

}