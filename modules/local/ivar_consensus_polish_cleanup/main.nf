//Small mofule that cleans up the final consensus FASTA file by removing leading and trailing Ns and the 'Consensus_' prefix
process IVAR_CONSENSUS_POLISH_CLEANUP {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(fasta, stageAs: 'ivar_consensus_input.fa')

    output:
    tuple val(meta), path("*.fa"), emit: fasta
    tuple val("${task.process}"), val('ivar_consensus_polish_cleanup'), val('1.0.0'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.final"
    """
    cp $fasta ${prefix}.fa
    sed -i 's/Consensus_//;s/^\\(N\\)\\{1,\\}//g;s/\\(N\\)\\{1,\\}\$//g' ${prefix}.fa
    """
}