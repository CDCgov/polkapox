/*
 * Thin compatibility layer for an originally modified nf-core IVAR_CONSENSUS
 * module.
 *
 * The old polishing process accepted tuple(meta, bam, fasta), ran iVar with a
 * `.final` prefix, optionally emitted mpileup, and cleaned the consensus FASTA
 * by removing the `Consensus_` header prefix and trimming terminal Ns. This
 * wrapper preserves that de novo polishing interface while delegating consensus
 * generation to the upstream nf-core IVAR_CONSENSUS process; only the FASTA
 * cleanup remains local.
 */

include { IVAR_CONSENSUS as IVAR_CONSENSUS_POLISH_RUN } from '../../../modules/nf-core/ivar/consensus/main'

workflow IVAR_CONSENSUS_POLISH {
    take:
    bam_fasta
    save_mpileup

    main:
    ch_versions = Channel.topic('versions')

    ch_bam = bam_fasta.map { meta, bam, fasta -> [ meta, bam ] }
    ch_fasta = bam_fasta.map { meta, bam, fasta -> fasta }

    IVAR_CONSENSUS_POLISH_RUN (
        ch_bam,
        ch_fasta,
        save_mpileup
    )

    IVAR_CONSENSUS_POLISH_CLEANUP (
        IVAR_CONSENSUS_POLISH_RUN.out.fasta
    )

    emit:
    fasta    = IVAR_CONSENSUS_POLISH_CLEANUP.out.fasta
    mpileup  = IVAR_CONSENSUS_POLISH_RUN.out.mpileup
    versions = ch_versions
}

process IVAR_CONSENSUS_POLISH_CLEANUP {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(fasta, stageAs: 'ivar_consensus_input.fa')

    output:
    tuple val(meta), path("*.fa"), emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.final"
    """
    cp $fasta ${prefix}.fa
    sed -i 's/Consensus_//;s/^\\(N\\)\\{1,\\}//g;s/\\(N\\)\\{1,\\}\$//g' ${prefix}.fa
    """
}
