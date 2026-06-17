//include { SAMTOOLS_FLAGSTAT_DENOVO                     } from './samtools_flagstat_denovo'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_DENOVO } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_DENOVO } from '../../modules/local/samtools_coverage/samtools_coverage'
include { UNICYCLER                                     } from '../../modules/nf-core/unicycler/main'
include { BANDAGE                                       } from '../../modules/nf-core/bandage/image/main'
include { GRAPH_RECON                                   } from '../../modules/local/graph_reconstruct/graph_reconstruct'
//include { BWA_MEM_DENOVO                                } from './bwa_mem_denovo'
include { BWA_DENOVO                                     } from '../../modules/local/bwa_denovo/bwa_denovo'
include { PUBLISH_CONTIGS                               } from '../../modules/local/publish_contigs/publish_contigs'
include { MUMMER                                        } from '../../modules/nf-core/mummer/main'
include { QUAST                                         } from '../../modules/nf-core/quast/main'
include { IVAR_CONSENSUS_POLISH                         } from './ivar_consensus_polish'

workflow DENOVO {

    take:
    trimmed_fastq

    main:

    ch_versions = Channel.empty()

    //
    // Module: run Unicycler
    //

    UNICYCLER (
        trimmed_fastq
    )
    ch_versions = ch_versions.mix(UNICYCLER.out.versions)
    ch_gfa = UNICYCLER.out.gfa

    //
    // Module: Publish Bandage PNG Plot
    //
    BANDAGE (
        ch_gfa,
    )

    //
    // Module: Genome Reconstruction from Unicycler GFA
    //
    GRAPH_RECON (
        ch_gfa,
    )
    ch_versions = ch_versions.mix(GRAPH_RECON.out.versions)
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
    ch_versions = ch_versions.mix(BWA_DENOVO.out.versions)
    
    //BWA_MEM_DENOVO (
    //    ch_denovo_joined,
    //    true
    //)
    //ch_mapped_denovo = BWA_MEM_DENOVO.out.bam
    //ch_mapped_denovo_flagstat = BWA_MEM_DENOVO.out.bambai
    //ch_versions = ch_versions.mix(BWA_MEM_DENOVO.out.versions)

    //
    // Module: Calculate statistics for de novo bwa mapping
    //
    SAMTOOLS_FLAGSTAT_DENOVO (
        ch_mapped_denovo_flagstat
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT_DENOVO.out.versions)

    SAMTOOLS_COVERAGE_DENOVO (
        ch_mapped_denovo_coverage
    )

    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE_DENOVO.out.versions)    
    

    //
    // Module: Polish assembly with IVAR Consensus
    //
    ch_polishing_input = ch_mapped_denovo.join(ch_gfa_forpolishing, by: 0)
    IVAR_CONSENSUS_POLISH (
        ch_polishing_input,
        true
    )
    ch_versions = ch_versions.mix(IVAR_CONSENSUS_POLISH.out.versions)
    //ch_gfapolish_compare = IVAR_CONSENSUS_POLISH.out.fasta
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
    ch_mummer = MUMMER.out.summary

    //
    // Module: run QUAST for assembly stats
    //
    QUAST (
        GRAPH_RECON.out.unicycler_contigs.collect{it[1]}.ifEmpty([]),
        true,
        true
    )
    ch_versions = ch_versions.mix(QUAST.out.versions)
    
    emit:
    quast_tsv = QUAST.out.tsv
    flagstat = SAMTOOLS_FLAGSTAT_DENOVO.out.flagstat
    coverage = SAMTOOLS_COVERAGE_DENOVO.out.coverage
    graph_recon_log = GRAPH_RECON.out.log
    gfa_assembly = GRAPH_RECON.out.gfa_assembly
    mummer_summary = MUMMER.out.summary
    fasta = IVAR_CONSENSUS_POLISH.out.fasta
    mpileup = IVAR_CONSENSUS_POLISH.out.mpileup
    versions      = ch_versions // channel: [ versions.yml ]
}
