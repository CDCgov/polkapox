//
// Generate reference genome related files for analysis
//

include { BWA_INDEX     } from '../../modules/nf-core/modules/bwa/index/main'

workflow PREPARE_GENOME {
    main:
    
    //
    // generate index for bwa
    // 

    ch_versions = Channel.empty()
    ch_fasta = file(params.fasta)

    BWA_INDEX ( ch_fasta )
    
    ch_bwa_index = BWA_INDEX.out.index
    ch_versions      = ch_versions.mix(BWA_INDEX.out.versions)

    emit:
    
    fasta                = ch_fasta            // path: genome.fasta
    bwa_index            = ch_bwa_index        // path: bwa/index/
    versions             = ch_versions         // channel: [ versions.yml ]
}