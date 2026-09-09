include { KRAKEN2_KRAKEN2                               } from '../../../modules/nf-core/kraken2/kraken2/main'
include { FASTP                                         } from '../../../modules/nf-core/fastp/main'
include { KRAKENTOOLS_EXTRACTKRAKENREADS                } from '../../../modules/nf-core/krakentools/extractkrakenreads/main'

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

    KRAKEN2_KRAKEN2 (
        input_reads,
        ch_kraken2_db,
        true,
        true
    )
  
    //using native extract reads instead of subseq
    KRAKENTOOLS_EXTRACTKRAKENREADS (
        "10242 10244",                                          // orthopox taxon IDs (space-separated) TODO: have as user input
        KRAKEN2_KRAKEN2.out.classified_reads_assignment,        
        KRAKEN2_KRAKEN2.out.classified_reads_fastq,             
        KRAKEN2_KRAKEN2.out.report                              
    )

    FASTP (
        KRAKENTOOLS_EXTRACTKRAKENREADS.out.extracted_kraken2_reads.map { meta, files -> [meta, files, []] },
        false, // writes reads that pass trimming
        false,
        false
    )

    emit:
    trimmed_fastq = FASTP.out.reads 
    json = FASTP.out.json
    kraken2_report = KRAKEN2_KRAKEN2.out.report
    classified_reads_assignment = KRAKEN2_KRAKEN2.out.classified_reads_assignment
    orthopox_reads = KRAKENTOOLS_EXTRACTKRAKENREADS.out.extracted_kraken2_reads
    versions      = ch_versions // channel: [ versions.yml ]

}