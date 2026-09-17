//Small module that cleans up the final consensus FASTA file by removing leading and trailing Ns and the 'Consensus_' prefix
process IVAR_CONSENSUS_POLISH_CLEANUP {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(fasta, stageAs: 'ivar_consensus_input.fa')
    val correct_ns

    output:
    tuple val(meta), path("*.fa"), emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.final"

    if (correct_ns) {
        """
        cp $fasta ${prefix}.fa
        sed -i 's/Consensus_//;s/polished/final/;s/^\\(N\\)\\{1,\\}//g;s/\\(N\\)\\{1,\\}\$//g' ${prefix}.fa
        """
    } else {
        """
        cp $fasta ${prefix}.fa
        sed -i 's/Consensus_//g;s/${meta.id}/${prefix}/g' ${prefix}.fa
        """
    }

}