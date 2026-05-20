/*
 * Thin compatibility layer for an originally modified nf-core BWA_MEM module.
 *
 * The local de novo path used to run a custom BWA_MEM process that accepted
 * tuple(meta, fasta, reads), built a BWA index from the reconstructed FASTA in
 * the same task, ran `bwa mem`, and emitted BAM plus BAI outputs. This wrapper
 * preserves that de novo-facing interface while delegating the tool execution
 * to upstream nf-core modules: BWA_INDEX and BWA_MEM.
 */

include { BWA_INDEX as BWA_INDEX_DENOVO     } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM   as BWA_MEM_DENOVO_ALIGN } from '../../modules/nf-core/bwa/mem/main'

workflow BWA_MEM_DENOVO {
    take:
    fasta_reads
    sort_bam

    main:
    ch_fasta = fasta_reads.map { meta, fasta, reads -> fasta }
    ch_reads = fasta_reads.map { meta, fasta, reads -> [ meta, reads ] }

    BWA_INDEX_DENOVO (
        ch_fasta
    )

    BWA_MEM_DENOVO_ALIGN (
        ch_reads,
        BWA_INDEX_DENOVO.out.index,
        sort_bam
    )

    ch_versions = BWA_INDEX_DENOVO.out.versions
        .mix(BWA_MEM_DENOVO_ALIGN.out.versions)

    emit:
    bam      = BWA_MEM_DENOVO_ALIGN.out.bam
    bambai   = BWA_MEM_DENOVO_ALIGN.out.bambai
    bai      = BWA_MEM_DENOVO_ALIGN.out.bai
    versions = ch_versions
}
