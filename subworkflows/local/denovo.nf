include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_DENOVO } from '../../modules/nf-core/modules/samtools/flagstat/main.denovo'
include { UNICYCLER                                     } from '../../modules/nf-core/modules/unicycler/main'
include { GRAPH_RECON                                   } from '../../modules/local/graph_reconstruct'
include { BWA_MEM as BWA_MEM_DENOVO                     } from '../../modules/nf-core/modules/bwa/mem/main_denovo'
include { PUBLISH_CONTIGS                               } from '../../modules/local/publish_contigs.nf'
include { MUMMER                                        } from '../../modules/nf-core/modules/mummer/main'
include { QUAST                                         } from '../../modules/nf-core/modules/quast/main'
include { IVAR_CONSENSUS as IVAR_CONSENSUS_POLISH       } from '../../modules/nf-core/modules/ivar/consensus/main_denovo'

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

    BWA_MEM_DENOVO (
        ch_denovo_joined,
        true
    )
    ch_mapped_denovo = BWA_MEM_DENOVO.out.bam
    ch_mapped_denovo_flagstat = BWA_MEM_DENOVO.out.bambai
    ch_versions = ch_versions.mix(BWA_MEM_DENOVO.out.versions)

    //
    // Module: Calculate statistics for de novo bwa mapping
    //
    SAMTOOLS_FLAGSTAT_DENOVO (
        ch_mapped_denovo_flagstat
    )

    ch_polishing_input = ch_mapped_denovo.join(ch_gfa_forpolishing, by: 0)

    //
    // Module: Polish assembly with IVAR Consensus
    //
    IVAR_CONSENSUS_POLISH (
        ch_polishing_input,
        false
    )
    ch_gfapolish_compare = IVAR_CONSENSUS_POLISH.out.fasta
    ch_tocompare = ch_gfaassm_compare.join(ch_gfapolish_compare, by: 0)
    
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
    graph_recon_log = GRAPH_RECON.out.log
    gfa_assembly = GRAPH_RECON.out.gfa_assembly
    mummer_summary = MUMMER.out.summary
    fasta = IVAR_CONSENSUS_POLISH.out.fasta
    mpileup = IVAR_CONSENSUS_POLISH.out.mpileup
    versions      = ch_versions // channel: [ versions.yml ]
}