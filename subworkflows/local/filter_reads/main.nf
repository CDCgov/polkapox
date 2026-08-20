
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
include { KRAKEN2_KRAKEN2                               } from '../../../modules/nf-core/kraken2/kraken2/main'
include { FASTP                                         } from '../../../modules/nf-core/fastp/main'
include { SEQTK_SUBSEQ                                  } from '../../../modules/nf-core/seqtk/subseq/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow READ_FILTER {
    take:
    input_reads

    main:

    ch_versions = Channel.topic('versions')

    //
    // MODULE: Run Kraken to keep only orthopox reads 
    //
    ch_kraken2_db = file(params.kraken_db, checkIfExists: true)

    //TODO replace with KRAKEN2_KRAKEN2
    KRAKEN2_KRAKEN2 (
        input_reads,
        ch_kraken2_db,
        true,
        true
    )
    ch_kreads = KRAKEN2_KRAKEN2.out.classified_reads_fastq
    ch_orthoreads = KRAKEN2_KRAKEN2.out.classified_reads_assignment

    SEQTK_SUBSEQ (
        ch_kreads,
        ch_orthoreads
    )
    ch_filt_fastq = SEQTK_SUBSEQ.out.sequences
    
    //
    // MODULE: Run Fastp
    //

    FASTP (
        ch_filt_fastq,
        false, // writes reads that pass trimming
        false,
        false
    )

    emit:
    trimmed_fastq = FASTP.out.reads 
    json = FASTP.out.json
    kraken2_report = KRAKEN2_KRAKEN2.out.report
    classified_reads_assignment = KRAKEN2_KRAKEN2.out.classified_reads_assignment
    seqtk_reads = SEQTK_SUBSEQ.out.sequences
    versions      = ch_versions // channel: [ versions.yml ]

}