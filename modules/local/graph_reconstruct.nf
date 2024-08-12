process GRAPH_RECON {

    conda (params.enable_conda ? "conda-forge::python=3.10.6 conda-forge::biopython conda-forge::numpy conda-forge::mkl conda-forge::networkx bioconda::gfapy bioconda::blast" : null)
    container "staphb/polkapox"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path('*asm.fasta')       , emit: gfa_assembly optional true
    path '*longest.fasta'                     , emit: gfa_longest optional true
    tuple val(meta), path('*.contigs.fasta')  , emit: unicycler_contigs optional true
    path '*.log'                              , emit: log
    path '*.summary'                          , emit: summary optional true
    path "versions.yml"                       , emit: versions

    shell: // This script is bundled with the pipeline, in polkapox/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gfa_unzip = "$gfa".replaceAll(/.gz/, "")
    """
    #gunzip -f $gfa
    #cat $gfa_unzip | awk 'BEGIN { FS="\\t" } /^S/{ if( length(\$3) >= 1) print ">Contig"\$2"_len"substr(\$4,6)"_cov"substr(\$5,6,5)"\\n"\$3}' | fold > ${prefix}.contigs.fasta
    
    AssemblyGraph_gfaPy.py \\
        -i $gfa \\
        -r "$projectDir/assets/MPXV-UK_P2.noN_39086_40204.fasta" \\
        -o .

    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
