include { BWA_MEM                                       } from '../../../modules/nf-core/bwa/mem/main'
include { IVAR_CONSENSUS                                } from '../../../modules/nf-core/ivar/consensus/main'
include { IVAR_VARIANTS                                 } from '../../../modules/nf-core/ivar/variants/main'
include { VARIANT_CONVERT                               } from '../../../modules/local/variant_convert/main'
include { SAMTOOLS_SORT                                 } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX                                } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT                             } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_DEPTH                                } from '../../../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_FAIDX                                } from '../../../modules/nf-core/samtools/faidx/main'
include { VCFTOOLS as VCFTOOLS_IVAR                     } from '../../../modules/nf-core/vcftools/main'
include { SUMMARIZE_TSV                                 } from '../../../modules/local/summarize_tsv/main'
include { AGGREGATE_TSVS                                } from '../../../modules/local/aggregate_tsvs/main'

workflow REFBASED {
    take: 
        ch_trimmed_fastq_bwa
        ch_bwa_index

    main: 
        ch_versions = Channel.topic('versions')

        //
        // Module: run BWA MEM alignment
        //
        BWA_MEM (
            ch_trimmed_fastq_bwa,
            ch_bwa_index,
            [[:], []], //fasta only required for cram output
            true //sort the bam file
        )

        //new version of BWA_MEM does not support index
        SAMTOOLS_INDEX(BWA_MEM.out.bam)

        // Join BAM with its index for downstream tools
        ch_bam_bai = BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.index)

        SAMTOOLS_FLAGSTAT (
            ch_bam_bai
        )

        SAMTOOLS_DEPTH (
            ch_bam_bai.map { meta, bam, bai -> [ meta, bam, bai, [] ] } //ie. no interval files
        )

        //
        // Module: run ivar
        //
    
        IVAR_CONSENSUS (
            BWA_MEM.out.bam,
            params.fasta,
            true
        )

        // Generate .fai index if not provided by the user
        if ( params.fai ) {
            ch_fai = file(params.fai, checkIfExists: true)
        } else {
            SAMTOOLS_FAIDX (
                [ [id:'index_fasta'], file(params.fasta) , [] ],
                false
            )
            ch_fai = SAMTOOLS_FAIDX.out.fai.map { meta, fai -> fai }
        }

        IVAR_VARIANTS (
            BWA_MEM.out.bam,
            params.fasta,
            ch_fai,
            params.gff ?: [], //if the user passes gff file
            false
        )
        ch_ivar_out = IVAR_VARIANTS.out.tsv
        
        VARIANT_CONVERT (
            ch_ivar_out,
            params.af_cutoff
        )
        ch_ivar_vcf = VARIANT_CONVERT.out.vcf
    
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

            //
            // Module: aggregate ivar tsv summary files
            //
            ch_aggregate_tsvs = Channel.empty()
            ch_aggregate_tsvs = ch_aggregate_tsvs.mix(SUMMARIZE_TSV.out.vars.collect{it[1]}.ifEmpty([]))

            AGGREGATE_TSVS (
                ch_aggregate_tsvs,
                true
            )
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