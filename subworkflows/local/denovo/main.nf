include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_DENOVO } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_DENOVO } from '../../../modules/local/samtools_coverage/main'
include { UNICYCLER                                     } from '../../../modules/nf-core/unicycler/main'
include { BANDAGE_IMAGE                                       } from '../../../modules/nf-core/bandage/image/main'
include { GRAPH_RECON                                   } from '../../../modules/local/graph_reconstruct/main'
include { BWA_DENOVO                                    } from '../../../modules/local/bwa_denovo/main'
include { PUBLISH_CONTIGS                               } from '../../../modules/local/publish_contigs/main'
include { MUMMER                                        } from '../../../modules/nf-core/mummer/main'
include { QUAST                                         } from '../../../modules/nf-core/quast/main'
include { IVAR_CONSENSUS_POLISH                         } from '../ivar_consensus_polish/main'

workflow DENOVO {

    take:
    trimmed_fastq

    main:

    ch_versions = Channel.topic('versions')

    //
    // Module: run Unicycler
    //

    UNICYCLER (
        trimmed_fastq.map { meta, reads -> [meta, reads, []] }
    )
    ch_gfa = UNICYCLER.out.gfa

    //
    // Module: Publish Bandage PNG Plot
    //
    BANDAGE_IMAGE (
        ch_gfa,
    )

    //
    // Module: Genome Reconstruction from Unicycler GFA
    //
    GRAPH_RECON (
        ch_gfa,
    )
    ch_graph_fasta = GRAPH_RECON.out.gfa_assembly
    ch_gfaassm_compare = GRAPH_RECON.out.gfa_assembly
    ch_gfa_forpolishing = GRAPH_RECON.out.gfa_assembly
    ch_uni_contigs = GRAPH_RECON.out.unicycler_contigs
    ch_gfa_summary = GRAPH_RECON.out.summary

    //
    // Module: Align reads to reconstructed genome
    //
    ch_denovo_joined = ch_graph_fasta.join(trimmed_fastq, by: 0)

    BWA_DENOVO (
        ch_denovo_joined,
        true
    )
    ch_mapped_denovo = BWA_DENOVO.out.bam
    ch_mapped_denovo_flagstat = BWA_DENOVO.out.bambai
    ch_mapped_denovo_coverage = BWA_DENOVO.out.bambai

    //
    // Module: Calculate statistics for de novo bwa mapping
    //
    SAMTOOLS_FLAGSTAT_DENOVO (
        ch_mapped_denovo_flagstat
    )

    SAMTOOLS_COVERAGE_DENOVO (
        ch_mapped_denovo_coverage
    )    

    //
    // Module: Polish assembly with IVAR Consensus
    //
    ch_polishing_input = ch_mapped_denovo.join(ch_gfa_forpolishing, by: 0)
    IVAR_CONSENSUS_POLISH (
        ch_polishing_input,
        true
    )
    ch_tocompare = ch_gfaassm_compare.join(IVAR_CONSENSUS_POLISH.out.fasta, by: 0)
    
    //join unicycler contigs with the polished fasta, and only keep contigs if fasta doesn't exist
    ch_assemblies = ch_uni_contigs.join(IVAR_CONSENSUS_POLISH.out.fasta, remainder: true)
    | map { meta, contigs, fasta -> [meta + [ "final" : "draft"], fasta ?: contigs ]}

    //
    // Module: Publish contigs if IVAR Polished fasta doesn't exist
    //
    PUBLISH_CONTIGS (
       ch_assemblies
    )

    //
    // Module: run MUMMER for assembly stats
    //
    MUMMER (
        ch_tocompare
    )
    //ch_mummer = MUMMER.out.summary //No more summary in the new version - TODO: determine if this is important

    //
    // Module: run QUAST for assembly stats
    //
    QUAST (
        GRAPH_RECON.out.unicycler_contigs,
        [[:], []], //no reference fasta
        [[:], []] //no annotations
    )
    
    emit:
    quast_tsv       = QUAST.out.tsv
    flagstat        = SAMTOOLS_FLAGSTAT_DENOVO.out.flagstat
    coverage        = SAMTOOLS_COVERAGE_DENOVO.out.coverage
    graph_recon_log = GRAPH_RECON.out.log
    gfa_assembly    = GRAPH_RECON.out.gfa_assembly
    //mummer_summary  = MUMMER.out.summary
    fasta           = IVAR_CONSENSUS_POLISH.out.fasta
    mpileup         = IVAR_CONSENSUS_POLISH.out.mpileup
    versions        = ch_versions // channel: [ versions.yml ]
}
