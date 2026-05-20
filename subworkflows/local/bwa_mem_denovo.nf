/*
 * Thin compatibility layer for an originally modified nf-core BWA_MEM module.
 *
 * The local de novo path used to run a custom BWA_MEM process that accepted
 * tuple(meta, fasta, reads), built a BWA index from the reconstructed FASTA in
 * the same task, ran `bwa mem`, and emitted BAM plus BAI outputs. This wrapper
 * preserves that de novo-facing interface while delegating the tool execution
 * to upstream nf-core modules: BWA_INDEX, BWA_MEM, and SAMTOOLS_INDEX.
 */

include { BWA_INDEX     as BWA_INDEX_DENOVO       } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM       as BWA_MEM_DENOVO_ALIGN   } from '../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DENOVO } from '../../modules/nf-core/samtools/index/main'

workflow BWA_MEM_DENOVO {
    take:
    fasta_reads
    sort_bam

    main:
    ch_fasta = fasta_reads.map { meta, fasta, reads -> [ meta, fasta ] }
    ch_reads = fasta_reads.map { meta, fasta, reads -> [ meta, reads ] }
    ch_fasta_for_mem = fasta_reads.map { meta, fasta, reads -> [ meta, fasta ] }

    BWA_INDEX_DENOVO (
        ch_fasta
    )

    BWA_MEM_DENOVO_ALIGN (
        ch_reads,
        BWA_INDEX_DENOVO.out.index,
        ch_fasta_for_mem,
        sort_bam
    )

    SAMTOOLS_INDEX_DENOVO (
        BWA_MEM_DENOVO_ALIGN.out.bam
    )

    ch_bambai = BWA_MEM_DENOVO_ALIGN.out.bam.join(SAMTOOLS_INDEX_DENOVO.out.index)
    ch_versions = BWA_INDEX_DENOVO.out.versions_bwa
        .mix(BWA_MEM_DENOVO_ALIGN.out.versions_bwa)
        .mix(BWA_MEM_DENOVO_ALIGN.out.versions_samtools)
        .mix(SAMTOOLS_INDEX_DENOVO.out.versions_samtools)

    emit:
    bam      = BWA_MEM_DENOVO_ALIGN.out.bam
    bambai   = ch_bambai
    bai      = SAMTOOLS_INDEX_DENOVO.out.index
    versions = ch_versions
}
