include { BWA_MEM                                       } from '../../modules/nf-core/modules/bwa/mem/main'
include { IVAR_CONSENSUS as IVAR_CONSENSUS_BWA          } from '../../modules/nf-core/modules/ivar/consensus/main'
include { IVAR_VARIANTS                                 } from '../../modules/nf-core/modules/ivar/variants/main'
include { VARIANT_CONVERT                               } from '../../modules/local/variant_convert.nf'
include { SAMTOOLS_SORT                                 } from '../../modules/nf-core/modules/samtools/sort/main'
include { SAMTOOLS_FLAGSTAT                             } from '../../modules/nf-core/modules/samtools/flagstat/main'
include { SAMTOOLS_DEPTH                                } from '../../modules/nf-core/modules/samtools/depth/main'
include { VCFTOOLS as VCFTOOLS_IVAR                     } from '../../modules/nf-core/modules/vcftools/main'
include { SUMMARIZE_TSV                                 } from '../../modules/local/summarize_tsv.nf'
include { AGGREGATE_TSVS                                } from '../../modules/local/aggregate_tsvs.nf'

workflow REFBASED {
    take: 
        ch_trimmed_fastq_bwa
        ch_bwa_index

    main: 
        ch_versions = Channel.empty()

        //
        // Module: run BWA MEM alignment
        //
        BWA_MEM (
            ch_trimmed_fastq_bwa,
            ch_bwa_index,
            true
        )
        ch_bwa_aln = BWA_MEM.out.bam
        ch_bwa_ivar = BWA_MEM.out.bam
        ch_bwa_depth = BWA_MEM.out.bam
        ch_bwa_flagstat = BWA_MEM.out.bambai
        ch_bai_lofreq = BWA_MEM.out.bai
        ch_versions = ch_versions.mix(BWA_MEM.out.versions)

        SAMTOOLS_FLAGSTAT (
            ch_bwa_flagstat
        )

        SAMTOOLS_DEPTH (
            ch_bwa_depth
        )

        //
        // Module: run ivar
        //
    
        IVAR_CONSENSUS_BWA (
            ch_bwa_aln,
            params.fasta,
            true
        )
        ch_versions = ch_versions.mix(IVAR_CONSENSUS_BWA.out.versions)

        IVAR_VARIANTS (
            ch_bwa_ivar,
            true
        )
        ch_ivar_out = IVAR_VARIANTS.out.tsv
        ch_versions = ch_versions.mix(IVAR_VARIANTS.out.versions)
        
        VARIANT_CONVERT (
            ch_ivar_out,
            params.af_cutoff
        )
        ch_ivar_vcf = VARIANT_CONVERT.out.vcf
        ch_versions = ch_versions.mix(VARIANT_CONVERT.out.versions)
    
        if ( params.filter ) {
            //
            // Module: run summarize ivar tsv
            //
            params.coords_list=params.coords?.split(',') as List

            SUMMARIZE_TSV (
                ch_ivar_out,
                params.coords_list,
                true
            )
            ch_tsv_vars = SUMMARIZE_TSV.out.vars
            ch_versions = ch_versions.mix(SUMMARIZE_TSV.out.versions)

            //
            // Module: aggregate ivar tsv summary files
            //
            ch_aggregate_tsvs = Channel.empty()
            ch_aggregate_tsvs = ch_aggregate_tsvs.mix(SUMMARIZE_TSV.out.vars.collect{it[1]}.ifEmpty([]))

            AGGREGATE_TSVS (
                ch_aggregate_tsvs,
                true
            )
            ch_versions = ch_versions.mix(AGGREGATE_TSVS.out.versions)
        }
        else {
            ch_tsv_vars = Channel.empty()
        }

        emit:
            flagstat = SAMTOOLS_FLAGSTAT.out.flagstat
            depth_tsv = SAMTOOLS_DEPTH.out.tsv
            tsv_vars = ch_tsv_vars
            ivar_tsv = IVAR_VARIANTS.out.tsv
            versions      = ch_versions // channel: [ versions.yml ]
}