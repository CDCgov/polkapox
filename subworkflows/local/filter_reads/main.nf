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
    input_reads.view { it -> "input_reads: ${it}" }
    //TODO replace with KRAKEN2_KRAKEN2
    KRAKEN2_KRAKEN2 (
        input_reads,
        ch_kraken2_db,
        true,
        true
    )
  
    //No need for subseq sibce tge 

    FASTP (
        SEQTK_SUBSEQ.out.sequences.map { meta, files -> [meta, files, []] },
        false, // writes reads that pass trimming
        false,
        false
    )
    FASTP.out.reads.view { it -> "FASTP.out.reads: ${it}" }

    emit:
    trimmed_fastq = FASTP.out.reads 
    json = FASTP.out.json
    kraken2_report = KRAKEN2_KRAKEN2.out.report
    classified_reads_assignment = KRAKEN2_KRAKEN2.out.classified_reads_assignment
    seqtk_reads = SEQTK_SUBSEQ.out.sequences
    versions      = ch_versions // channel: [ versions.yml ]

}