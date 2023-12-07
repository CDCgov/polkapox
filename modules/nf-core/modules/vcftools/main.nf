process VCFTOOLS {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::vcftools=0.1.16" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcftools:0.1.16--he513fc3_4' :
        'quay.io/biocontainers/vcftools:0.1.16--he513fc3_4' }"

    input:
    tuple val(meta), path(variant_file)
    
    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}." + task.ext.caller
    def args_list = args.tokenize()

    def input_file = ("$variant_file".endsWith(".vcf")) ? "--vcf ${variant_file}" :
        ("$variant_file".endsWith(".vcf.gz")) ? "--gzvcf ${variant_file}" :
        ("$variant_file".endsWith(".bcf")) ? "--bcf ${variant_file}" : ''
    def bed_arg  = (args.contains('--bed')) ? "--bed ${bed}" :
        (args.contains('--exclude-bed')) ? "--exclude-bed ${bed}" :
        (args.contains('--hapcount')) ? "--hapcount ${bed}" : ''
    args_list.removeIf { it.contains('--bed') }
    args_list.removeIf { it.contains('--exclude-bed') }
    args_list.removeIf { it.contains('--hapcount') }

    """
    vcftools \\
        $input_file \\
        --bed "$projectDir/assets/MPXV-UK_P2.noN_drugres.bed" \\
        --recode \\
        --out $prefix \\
        --recode-INFO-all

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcftools: \$(echo \$(vcftools --version 2>&1) | sed 's/^.*VCFtools (//;s/).*//')
    END_VERSIONS
    """
}


