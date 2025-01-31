process PUBLISH_CONTIGS {
    tag "$meta.id"
    label 'process_low'

   conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bionumpy:0.2.17--pyha8f3691_0  ' :
        'quay.io/biocontainers/bionumpy:0.2.17--pyha8f3691_0' }"
       
    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*draft.fa"), optional: true
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def draft_asm = "${meta.id}.draft.fa"
    def contigs = "${meta.id}.contigs.fasta"
    def assembly = "${meta.id}.final.fa"
    """
    if [ ! -f $assembly ]
    then cp $contigs $draft_asm
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}

